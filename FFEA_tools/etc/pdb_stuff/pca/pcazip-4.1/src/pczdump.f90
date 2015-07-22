! Pczdump: extract information from an MD file compressed with pcazip
! version 3.2
! Copyright (C) Charlie Laughton 2008
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
program pczdump
use pcz_io
use mygetargs
implicit none

  type(pczfile)              :: pczf
  integer,parameter          :: STDERR=0, STDOUT=6
  integer                    :: i,ounit,nr,na,iv,ip,stat,ir
  character(160)             :: infile,outfile,pdbfile
  character(4)               :: testline
  real, allocatable          :: avg(:),evals(:),evec(:),proj(:)
  character(30), allocatable :: pdbline(:)
  logical                    :: dumpprj

  if(getflag('-help')) then
    call helptext(STDERR)
    call exit(0)
  end if

  dumpprj=.not.getflag('-noproj')

  if(.not.mygetarg('-i',infile)) then
    write(STDERR,*) 'Specify input file with -i'
    call exit(1)
  end if
  pczf=pczopen(1,infile)
  if(pczf%status /=0) then
    write(STDERR,*) 'Error with input file: ',trim(infile)
    write(STDERR,*) '      pczopen returns: ',pczf%status
!    call perror(trim(infile))
    call exit(1)
  end if

  if(mygetarg('-o',outfile)) then
    ounit=2
    open(ounit,file=outfile,status='new')
  else
    ounit=STDOUT
  end if

  if(mygetarg('-pdb',pdbfile)) then
    open(3,file=pdbfile,status='old')
    nr=0
    na=0
    stat=0
    do while(stat==0)
      read(3,'(a4)',iostat=stat) testline
      if(stat == 0) then
        nr=nr+1
        if(testline=='ATOM') na=na+1
      end if
    end do
    if(na /= pczf%natoms) then
      write(STDERR,*) 'Error: pdbfile contains wrong number of atoms'
      call exit(1)
    end if
    allocate(pdbline(nr))
    rewind(3)
    do i=1,nr
      read(3,'(a30)') pdbline(i)
    end do
    close(3)
  end if

  if(getflag('-avg')) then
    if(mygetarg('-pdb',pdbfile)) then
      call dumpavgpdb(pczf,ounit,pdbline,nr)
    else
      call dumpavgx(pczf,ounit)
    end if
  else if(getflag('-evals')) then
    call dumpevals(pczf,ounit)
  else if(getflag('-coll')) then
    call dumpcoll(pczf,ounit)
  else if(getflag('-info')) then
    call dumpinfo(pczf,ounit)
  else if(mygetarg('-evec',iv)) then
    call dumpevec(pczf,ounit,iv)
  else if(mygetarg('-maha',iv)) then
    call dumpmaha(pczf,ounit,iv)
  else if(mygetarg('-proj',ip)) then
    call dumpproj(pczf,ounit,ip)
  else if(mygetarg('-rms',ir)) then
    call dumprms(pczf,ounit,ir)
  else if(mygetarg('-fluc',iv)) then
    call dumpfluc(pczf,ounit,iv)
  else if(mygetarg('-anim',iv)) then
    if(mygetarg('-pdb',pdbfile)) then
      call dumpanimpdb(pczf,ounit,iv,pdbline,nr)
    else
      call dumpanimx(pczf,ounit,iv)
    end if
  else
    call dumpinfo(pczf,ounit)
    write(ounit,*) 'Average structure:'
    call dumpavgx(pczf,ounit)
    write(ounit,*) 'Eigenvalues:'
    call dumpevals(pczf,ounit)
    do iv=1,pczf%nvecs
      write(ounit,*) 'Eigenvector ',iv
      call dumpevec(pczf,ounit,iv)
      if(dumpprj) then
        write(ounit,*) 'Projections along eigenvector ',iv
        call dumpproj(pczf,ounit,iv)
      endif
    end do
  end if
  close(ounit)
  stat=pczclose(pczf)
end program pczdump

subroutine dumpavgpdb(pczf,iunit,pdbline,nr)
use pcz_io
implicit none

  type(pczfile)         :: pczf
  integer               :: iunit,nr
  character(30)         :: pdbline(nr)

  real                  :: avg(pczf%natoms*3)
  integer               :: na,i,stat

  avg=pczavg(pczf)
  na=0
  stat=0
  do i=1,nr
    if(pdbline(i)(1:4)=='ATOM') then
      na=na+1
      if(stat==0) write(iunit,'(a30,3f8.3)',iostat=stat) pdbline(i),avg(3*na-2:3*na)
    else
      if(stat==0) write(iunit,'(a)',iostat=stat) trim(pdbline(i))
    end if
  end do
end subroutine dumpavgpdb

subroutine dumpavgx(pczf,iunit)
use pcz_io
implicit none

  type(pczfile)        :: pczf
  integer              :: iunit,stat

  write(iunit,'(a80)',iostat=stat) pczf%title
  write(iunit,'(10f8.3)',iostat=stat) pczavg(pczf)
end subroutine dumpavgx

subroutine dumpevals(pczf,iunit)
use pcz_io
implicit none

  type(pczfile)       :: pczf
  integer             :: iunit

  write(iunit,'(f15.3)') pczevals(pczf)
end subroutine dumpevals

subroutine dumpcoll(pczf,iunit)
use pcz_io
implicit none

  type(pczfile)       :: pczf
  integer             :: iunit
  real                :: evec(pczf%natoms*3),coll(pczf%nvecs),r2
  integer             :: i,j

  do i=1,pczf%nvecs
    evec=pczevec(pczf,i)
    coll(i)=0.0
    do j=1,3*pczf%natoms,3
      r2=evec(j)*evec(j)+evec(j+1)*evec(j+1)+evec(j+2)*evec(j+2)
      coll(i)=coll(i)-r2*log(r2)
    end do
  end do
  coll=exp(coll)/pczf%natoms
  write(iunit,'(f15.3)') coll
end subroutine dumpcoll

subroutine dumpinfo(pczf,iunit)
use pcz_io
implicit none

  type(pczfile)       :: pczf
  integer             :: iunit

  write(iunit,*) 'Title:   ',trim(pczf%title)
  write(iunit,*) 'Atoms:   ',pczf%natoms
  write(iunit,*) 'Vectors: ',pczf%nvecs
  write(iunit,*) 'Frames:  ',pczf%nframes
  if(pczf%quality == 0.0) then
    write(iunit,*) 'Quality:  (unknown)'
  else
    write(iunit,'(" Quality: ",f5.1,"% (of ",f8.2,")")') pczf%quality,&
         100.0*sum(pczevals(pczf))/pczf%quality
  end if
end subroutine dumpinfo

subroutine dumpevec(pczf,iunit,iv)
use pcz_io
implicit none

  type(pczfile)       :: pczf
  integer             :: iunit,iv,stat

  write(iunit,'(10f8.4)',iostat=stat) pczevec(pczf,iv)
end subroutine dumpevec

subroutine dumpmaha(pczf,iunit,iv)
use pcz_io
implicit none

  type(pczfile)       :: pczf
  integer             :: iunit,iv,stat

  write(iunit,'(f8.4)',iostat=stat) pczmaha(pczf,iv)
end subroutine dumpmaha

subroutine dumpproj(pczf,iunit,ip)
use pcz_io
implicit none

  type(pczfile)      :: pczf
  integer            :: iunit,ip

  write(iunit,'(f8.3)') pczproj(pczf,ip)
end subroutine dumpproj

subroutine dumprms(pczf,iunit,ir)
use pcz_io
implicit none

  type(pczfile)        :: pczf
  integer              :: iunit,ir

  real                 :: s1(pczf%nvecs),s2(pczf%nvecs),rmsd
  integer              :: i

  if(ir /= 0) s2=pczscores(pczf,ir)
  do i=1,pczf%nframes
    rmsd=0.0
    s1=pczscores(pczf,i)
    if(ir /= 0) s1=s1-s2
    rmsd=sum(s1*s1)/pczf%natoms
    write(iunit,'(f8.3)') sqrt(rmsd)
  end do
end subroutine dumprms

subroutine dumpfluc(pczf,iunit,iv)
use pcz_io
implicit none

  type(pczfile)       :: pczf
  integer             :: iunit,iv

  write(iunit,'(f8.3)') pczfluc(pczf,iv)
end subroutine dumpfluc

subroutine dumpanimpdb(pczf,iunit,iv,pdbline,nr)
use pcz_io
implicit none

  type(pczfile)      :: pczf
  integer            :: iunit,iv,nr
  character(30)      :: pdbline(nr)

  real               :: avg(pczf%natoms*3),evec(pczf%natoms*3),x(pczf%natoms*3)
  real               :: proj(max(21,pczf%nframes)),rmin,rmax,rinc
  integer            :: na,i,j,stat

  avg=pczavg(pczf)
  evec=pczevec(pczf,iv)
  proj=pczproj(pczf,iv)

  rmin=minval(proj)
  rmax=maxval(proj)
  rinc=(rmax-rmin)*0.1
  proj(1)=0.0
  do i=2,6
   proj(i)=proj(i-1)+rinc
  end do
  do i=7,16
   proj(i)=proj(i-1)-rinc
  end do
  do i=17,20
   proj(i)=proj(i-1)+rinc
  end do

  write(iunit,'(a34,i4)',iostat=stat) 'REMARK   Animation of eigenvector ',iv
  do j=1,20
    x=avg+evec*proj(j)
    na=0
    write(iunit,'(a7,i5)') 'MODEL  ',j
    do i=1,nr
      if(pdbline(i)(1:4)=='ATOM') then
        na=na+1
        if(stat==0) write(iunit,'(a30,3f8.3)',iostat=stat) pdbline(i),x(3*na-2:3*na)
      else
        if(stat==0) write(iunit,'(a)',iostat=stat) trim(pdbline(i))
      end if
    end do
    if(stat==0) write(iunit,'(a6)') 'ENDMDL'
  end do
end subroutine dumpanimpdb

subroutine dumpanimx(pczf,iunit,iv)
use pcz_io
implicit none

  type(pczfile)      :: pczf
  integer            :: iunit,iv

  real               :: avg(pczf%natoms*3),evec(pczf%natoms*3),x(pczf%natoms*3)
  real               :: proj(max(21,pczf%nframes)),rmin,rmax,rinc
  integer            :: na,i,j,stat

  avg=pczavg(pczf)
  evec=pczevec(pczf,iv)
  proj=pczproj(pczf,iv)

  rmin=minval(proj)
  rmax=maxval(proj)
  rinc=(rmax-rmin)*0.1
  proj(1)=0.0
  do i=2,6
   proj(i)=proj(i-1)+rinc
  end do
  do i=7,16
   proj(i)=proj(i-1)-rinc
  end do
  do i=17,20
   proj(i)=proj(i-1)+rinc
  end do

  write(iunit,'(a26,i4)') 'Animation of eigenvector ',iv
  stat=0
  do j=1,20
    x=avg+evec*proj(j)
    if(stat==0) write(iunit,'(10f8.3)',iostat=stat) x
  end do
end subroutine dumpanimx

subroutine helptext(iunit)
implicit none
  integer   :: iunit

  write(iunit,"('*** pczdump version 3.2 ***')")
  write(iunit,"('Usage:'/'pcadump -i infile -o outfile [-info] [-avg]')")
  write(iunit,"('  [-evec iv] [-proj iv] [-rms iref] [-fluc iv] [-evals]')")
  write(iunit,"('  [-anim iv [-pdb pdbref]] [-coll]')")
end subroutine helptext

subroutine die(unit,msg)
implicit none
  integer            :: unit
  character(len=*)   :: msg

  write(unit,*) msg
  call exit(1)
end subroutine die
