! Pcz_io: general 'pcz' format file I/O routines
! Version 4.0 - added support for dcd format output courtesy of
! Yiming Chen
! Copyright (C) Charlie Laughton 2011
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License along
!   with this program; if not, write to the Free Software Foundation, Inc.,
!   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! Charlie Laughton
! School of Pharmacy
! University of Nottingham
! NG7 2RD
! UK
! charles.laughton@nottingham.ac.uk
module pcz_io
implicit none
  type :: pczfile
    character(160):: filename
    character(80) :: title
    character(4)  :: version
    integer       :: unit
    integer       :: status
    integer       :: natoms
    integer       :: nvecs
    integer       :: nframes
    real          :: quality
    real, pointer :: avg(:)
    real, pointer :: evals(:)
    real, pointer :: evecs(:,:)
    real, pointer :: proj(:,:)
  end type pczfile
  contains
    function pczcreate(un,fname,title,natoms,nframes,nvecs,quality,avg) &
        result(pczf)
    implicit none
    character(160) :: fname     ! filename
    character(80)  :: title     ! title from traj file
    integer        :: un        ! unit number
    integer        :: natoms    ! number of atoms
    integer        :: nvecs     ! number of eigenvectors
    integer        :: nframes   ! number of snapshots
    real           :: quality   ! quality of pczfile (% of total variance)
    real           :: avg(3*natoms) ! average structure
    type(pczfile)  :: pczf      ! pczfile structure returned
! specimen record contents
    integer           :: iol, iol1, iol2
! initialise the pczfile structure
    pczf%filename=fname
    pczf%unit=un
    pczf%title=title
    pczf%version='PCZ4'
    pczf%status=0
    pczf%natoms=natoms
    pczf%nvecs=nvecs
    pczf%nframes=nframes
    pczf%quality=quality
    allocate(pczf%avg(3*natoms),pczf%evals(nvecs),pczf%evecs(3*natoms,nvecs))
    allocate(pczf%proj(nframes,nvecs))

    pczf%avg=avg
!
! all done.
!
    end function pczcreate

    function pczclose(pczf) result(status)
    implicit none
    type(pczfile)  :: pczf
    integer        :: i,status
    integer(2)     :: ip(pczf%nframes)
    real           :: evec(3*pczf%natoms),eval,proj(pczf%nframes)
    real           :: p0,pinc
!
! first we convert the quality measure,
! specified originally as a percentage, into a sum of the
! eigenvalues. This is so that if the file is subsequently truncated
! it will still report the correct quality measure when opened
!
    if(pczf%quality > 0.0) pczf%quality=100.0*sum(pczevals(pczf))/pczf%quality
    open(pczf%unit,file=pczf%filename,form='unformatted',access='stream')
    write(pczf%unit,iostat=pczf%status) pczf%version,pczf%title, &
      pczf%natoms, pczf%nframes, pczf%nvecs, pczf%quality, 0,0,0,0,pczf%avg
    do i=1,pczf%nvecs
!
! pointer components of structures not allowed in I/O lists...
!
      evec=pczf%evecs(:,i)
      eval=pczf%evals(i)
      proj=pczf%proj(:,i)
      if(pczf%version=='PCZ6') then
!
! new super-compressed format
!
        pinc=(maxval(proj)-minval(proj))/65534
        p0=(maxval(proj)+minval(proj))/2
        proj=proj-p0
        ip=proj/pinc
        write(pczf%unit,iostat=pczf%status) evec,eval,p0,pinc,ip
      else
        write(pczf%unit,iostat=pczf%status) evec,eval,proj
      endif
    end do
    close(pczf%unit)
!
! free memory
!
    deallocate(pczf%avg)
    deallocate(pczf%evals)
    deallocate(pczf%evecs)
    deallocate(pczf%proj)
    status=pczf%status
    end function pczclose

    function pczwrite(pczf,ivec,evec,eval,proj) result(status)
    implicit none
    type(pczfile)   :: pczf
    integer         :: ivec,status
    real            :: evec(3*pczf%natoms),eval,proj(pczf%nframes)

    status=0
    if(ivec>pczf%nvecs) then
      status=-1
      return
    end if
    pczf%evecs(:,ivec)=evec
    pczf%evals(ivec)=eval
    pczf%proj(:,ivec)=proj

    end function pczwrite

    function pczopen(un,fname) result(pczf)
    use rscode
    implicit none
    character(160) :: fname ! filename
    integer        :: un    ! unit number
    type(pczfile)  :: pczf  ! pczfile structure returned
! specimen record contents
    real, allocatable :: evec(:),proj(:)
    integer(2),allocatable  :: ip(:)
    character(5),allocatable :: cbuffer(:)
    real              :: eval,pczfq
    real              :: p0,pinc
    integer           :: i,p, pczfn, iol, iol1, iol2, hasnames
    integer           :: atomnumber,residuenumber
    character(4)      :: atomname
    character(3)      :: residuename
    character(1)      :: chain
    logical           :: exists
    integer           :: j

! initialise the pczfile structure
    pczf%filename=fname
    pczf%unit=un
    pczf%title=''
    pczf%version='????'
    pczf%status=0
    pczf%natoms=0; pczf%nvecs=0; pczf%nframes=0; pczf%quality=0.0
!
! does it exist?
!
    inquire(file=fname,exist=exists)
    if(.not.exists) then
      write(0,*) trim(fname),': no such file?'
      call exit(1)
    end if
!
! first peek at the file to check it out
!
    open(un,file=fname,status='old',form='unformatted', access='stream',iostat=pczf%status)
    read(un) pczf%version
    close(un)
    if(pczf%version/='PCZ4'.and.pczf%version/='PCZ5' &
       .and.pczf%version/='PCZ0'.and.pczf%version/='PCZ6') then
!  
! try opening it unformatted direct access (PCZ2-style)
!
      close(un)
      inquire(iolength=iol) pczf%version
      open(un,file=fname,status='old',iostat=pczf%status,access='direct', &
      form='unformatted', recl=iol)
! give up if there is a problem
      if(pczf%status /= 0) return
! read key information
      read(un,rec=1,iostat=pczf%status) pczf%version
! give up if there is a problem
      if(pczf%status /= 0) return
! give up if not a recognised format
      if(pczf%version/='PCZ2') then
        write(0,*) trim(fname),':  unrecognised pcz file format'
        call exit(1)
      endif
      close(un)
    endif
!
! OK - now we know we have a known PCZ file format - proceed accordingly
!
    if(pczf%version=="PCZ0".or.pczf%version=="PCZ5") then
      open(un,file=fname,status='old')
      read(un,'(a4)') pczf%version
      read(un,'(a80)') pczf%title
      read(un,'(3i8,f8.3)') pczf%natoms,pczf%nframes,pczf%nvecs,pczf%quality
!
! now read in everything else...
!
! Begin by allocating all arrays, including extra ones as cannot use pointers
! in i/o sttements...
!
      allocate(pczf%avg(pczf%natoms*3),evec(pczf%natoms*3), STAT=pczf%status)
      if(pczf%status/=0) call die('Memory allocation error')
      allocate(pczf%evals(pczf%nvecs),STAT=pczf%status)
      if(pczf%status/=0) call die('Memory allocation error')
      allocate(pczf%evecs(3*pczf%natoms,pczf%nvecs),STAT=pczf%status)
      if(pczf%status/=0) call die('Memory allocation error')
      allocate(pczf%proj(pczf%nframes,pczf%nvecs),proj(pczf%nframes),STAT=pczf%status)
      if(pczf%status/=0) call die('Memory allocation error')
      if(pczf%version=='PCZ0') then
!
! this is the vanilla formatted format
!
        read(un,'(10f8.3)') evec
        pczf%avg=evec
        do i=1,pczf%nvecs
          read(un,'(10f8.3)',iostat=pczf%status) evec
          if(pczf%status/=0) return
          pczf%evecs(:,i)=evec
          read(un,'(10f8.3)',iostat=pczf%status) eval
          if(pczf%status/=0) return
          pczf%evals(i)=eval
          read(un,'(10f8.3)',iostat=pczf%status) proj
          if(pczf%status/=0) return
          pczf%proj(:,i)=proj
        end do
        deallocate(evec,proj)
      else if(pczf%version=='PCZ5') then
!
! this is the rare base64 encoded format
!
        allocate(cbuffer(3*pczf%natoms+1+pczf%nframes),STAT=pczf%status)
        if(pczf%status/=0) call die('Memory allocation error')
        read(pczf%unit,'(16a5)',iostat=pczf%status) cbuffer(1:3*pczf%natoms)
        do i=1,3*pczf%natoms
          pczf%avg(i)=stor(cbuffer(i))
        end do
        do p=1,pczf%nvecs
          read(pczf%unit,'(16a5)',iostat=pczf%status) cbuffer
          if(pczf%status/=0) return
          do i=1,3*pczf%natoms
            pczf%evecs(i,p)=stor(cbuffer(i))
          end do
          pczf%evals(p)=stor(cbuffer(3*pczf%natoms+1))
          do i=1,pczf%nframes
            pczf%proj(i,p)=stor(cbuffer(i+1+3*pczf%natoms))
          end do
        end do
        deallocate(cbuffer)
        close(un)
      endif
!
! can now correct quality measure - is stored as
! the sum of all eigenvalues, rather than a percentage
!
      if(pczf%quality > 0.0) pczf%quality=100.0*sum(pczevals(pczf))/pczf%quality
      return
    else if(pczf%version=='PCZ2') then
!
! this is the fortran unformatted direct access format
!
      inquire(iolength=iol) pczf%version,pczf%title, pczf%natoms,pczf%nframes, &
        pczf%nvecs, pczf%quality
      open(un,file=fname,status='old',iostat=pczf%status,access='direct', &
        form='unformatted', recl=iol)
      read(un,rec=1,iostat=pczf%status) pczf%version,pczf%title, &
          pczf%natoms,pczf%nframes,pczf%nvecs, pczf%quality
      if(pczf%status /= 0) then
        return
      end if
!
! correct number of atoms in PCZ2 format
      pczf%natoms=pczf%natoms/3
!
! we now have the information required to reopen the file with
! the correct record length
!
      close(un)
      allocate(pczf%avg(pczf%natoms*3),evec(pczf%natoms*3), STAT=pczf%status)
      if(pczf%status/=0) call die('Memory allocation error')
      allocate(pczf%evals(pczf%nvecs),STAT=pczf%status)
      if(pczf%status/=0) call die('Memory allocation error')
      allocate(pczf%evecs(3*pczf%natoms,pczf%nvecs),STAT=pczf%status)
      if(pczf%status/=0) call die('Memory allocation error')
      allocate(pczf%proj(pczf%nframes,pczf%nvecs),proj(pczf%nframes),STAT=pczf%status)
      if(pczf%status/=0) call die('Memory allocation error')
! we need to check which is longer - the header structure or the data structure
      inquire(iolength=iol1) evec,eval,proj
      inquire(iolength=iol2) pczf%version,pczf%title,pczf%natoms,pczf%nvecs, &
        pczf%nframes, pczf%quality, evec
      iol=max(iol1,iol2)
      open(un,file=fname,status='old',iostat=pczf%status, &
        access='direct',form='unformatted',recl=iol)
      read(pczf%unit,rec=1,iostat=pczf%status) pczf%version,pczf%title, &
        pczfn,pczf%nframes,pczf%nvecs, pczfq, evec
      if(pczf%status/=0) return
      pczf%avg=evec
      do p=1,pczf%nvecs
        read(un,rec=p+1,iostat=pczf%status) evec,eval,proj
        if(pczf%status/=0) return
        pczf%evecs(:,p)=evec
        pczf%evals(p)=eval
        pczf%proj(:,p)=proj
      end do
    else
!
! PCZ4/6 format - pure binary stream format
!
      open(un,file=fname,status='old',form='unformatted',access='stream')
      read(un) pczf%version,pczf%title,pczf%natoms,pczf%nframes
      read(un) pczf%nvecs, pczf%quality
!
! read the Barcelona-version stuff we currently ignore...
      read(un) iol,iol1,iol2,hasnames
      if(hasnames==1) then
        do i=1,pczf%natoms
          read(un) atomnumber,atomname,residuenumber,residuename,chain
        end do
      endif
      allocate(pczf%avg(pczf%natoms*3),evec(pczf%natoms*3), STAT=pczf%status)
      if(pczf%status/=0) call die('Memory allocation error')
      allocate(pczf%evals(pczf%nvecs),STAT=pczf%status)
      if(pczf%status/=0) call die('Memory allocation error')
      allocate(pczf%evecs(3*pczf%natoms,pczf%nvecs),STAT=pczf%status)
      if(pczf%status/=0) call die('Memory allocation error')
      allocate(pczf%proj(pczf%nframes,pczf%nvecs),proj(pczf%nframes),STAT=pczf%status)
      if(pczf%status/=0) call die('Memory allocation error')
      if(pczf%version=='PCZ6') then
        allocate(ip(pczf%nframes),STAT=pczf%status)
        if(pczf%status/=0) call die('Memory allocation error')
      endif
      read(un) evec
      pczf%avg=evec
      do p=1,pczf%nvecs
        if(pczf%version=='PCZ4') then
          read(un,iostat=pczf%status) evec,eval,proj
        else
          read(un,iostat=pczf%status) evec,eval,p0,pinc,ip
          proj=p0
          proj=proj+ip*pinc
        endif
        if(pczf%status/=0) return
        pczf%evecs(:,p)=evec
        pczf%evals(p)=eval
        pczf%proj(:,p)=proj
      end do
!      if(pczf%quality > 0.0) pczf%quality=100.0*sum(pczevals(pczf))/pczf%quality
      close(un)
    endif
!
! can now correct quality measure - is stored as
! the sum of all eigenvalues, rather than a percentage
!
    if(pczf%quality > 0.0) pczf%quality=100.0*sum(pczevals(pczf))/pczf%quality
!
!  all done.
!
    if(allocated(evec)) deallocate(evec,proj)
    if(allocated(ip)) deallocate(ip)
    end function pczopen

    function pczavg(pczf) result(avg)
    implicit none
    type(pczfile) :: pczf
! next 2 are to stop overwriting numbers corrected in pczopen
    real          :: avg(3*pczf%natoms)

    avg=pczf%avg
    end function pczavg

    function pczevec(pczf,ivec) result(evec)
    implicit none
    type(pczfile)  :: pczf
    integer        :: ivec
    real           :: evec(3*pczf%natoms)

    evec=pczf%evecs(:,ivec)
    end function pczevec

    function pczevals(pczf) result(evals)
    implicit none
    type(pczfile)    :: pczf
    real             :: evals(pczf%nvecs)

    evals=pczf%evals
    end function pczevals

    function pczproj(pczf,ivec) result(proj)
    implicit none
    type(pczfile)  :: pczf
    integer        :: ivec
    real           :: proj(pczf%nframes)

    proj=pczf%proj(:,ivec)
    end function pczproj

    function pczmaha(pczf,nvecs) result(maha)
    implicit none
    type(pczfile)  :: pczf
    integer        :: nvecs
    real           :: maha(pczf%nframes)

    real           :: scores(pczf%nvecs),evals(pczf%nvecs)
    integer        :: i

    evals=pczevals(pczf)
    do i=1,pczf%nframes
      scores=pczscores(pczf,i)
      scores=scores*scores/evals
      maha(i)=sqrt(sum(scores(1:nvecs)))
    end do
    end function pczmaha

    function pczscores(pczf,iframe) result(scores)
    implicit none
    type(pczfile)  :: pczf
    integer        :: iframe
    real           :: scores(pczf%nvecs)
    real           :: proj(pczf%nframes)
    integer        :: i

    do i=1,pczf%nvecs
      proj=pczproj(pczf,i)
      scores(i)=proj(iframe)
    end do
    end function pczscores

    function pczsnap(pczf,iframe,iv1,iv2) result(snap)
    implicit none
    type(pczfile)   :: pczf
    integer         :: iframe
    integer, optional :: iv1,iv2
    real            :: snap(3*pczf%natoms)
    real            :: scores(pczf%nvecs)
    integer         :: ivec,liv1,liv2

    snap=pczavg(pczf)
    if(pczf%status /= 0) return
    if(iframe > pczf%nframes) return
    if(iframe == 0) return
    if(pczf%status /= 0) return
    if(.not.present(iv1)) then
      liv1=1
    else
      liv1=iv1
    endif
    if(.not.present(iv2)) then
      liv2=pczf%nvecs
    else
      liv2=min(iv2,pczf%nvecs)
    endif
    do ivec = liv1,liv2
      snap=snap+pczf%evecs(:,ivec)*pczf%proj(iframe,ivec)
      if(pczf%status /= 0) exit
    end do
    end function pczsnap

    function pczfluc(pczf,iv) result(fluc)
    implicit none
    type(pczfile)    :: pczf
    real             :: fluc(pczf%natoms), evec(pczf%natoms*3)
    integer          :: iv, i,j

    evec=pczevec(pczf,iv)
    if(pczf%status /= 0) return
    do i=1,pczf%natoms
      j=3*i-2
      fluc(i)=evec(j)*evec(j)+evec(j+1)*evec(j+1)+evec(j+2)*evec(j+2)
      fluc(i)=sqrt(fluc(i))
    end do
    end function pczfluc

end module pcz_io

