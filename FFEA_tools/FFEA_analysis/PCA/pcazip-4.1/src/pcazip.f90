! Pcazip: compress MD trajectory files using a PCA-based method
! Version 4.1
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
program pcazip
use mygetargs
use x_io
use pcz_io
use trajutils
use eigenutils
implicit none

  integer,parameter             :: STDERR=0, IDQ=90
  real, allocatable             :: x(:,:),xt(:),pn(:,:),evec(:),proj(:),xtmp(:)
  real, allocatable             :: mass(:)
  logical, allocatable          :: mask(:)
  real            , allocatable :: cm(:,:),w(:),z(:,:)
  integer, allocatable          :: freezeidx(:)
  character(160)                :: mdfile,outfile,maskfile,massfile
  character(80)                 :: title
  logical                       :: hasbox,ok,verbose,binary,album
  logical                       :: exists,dofast,lowmem
  integer                       :: trajframes,i,j,k,i2,nf2
  integer                       :: nvecs,iqual,natom,natom3,nframes,ierr
  integer                       :: nmask,nfreeze,ip1,ip2,stat
  real                          :: rframes
  real                          :: vsum, vtot,eval,inc,ri2
  type(xfile)                   :: xf
  type(pczfile)                 :: pczf
  type(xalbum)                  :: xalbf

  call f77_molfile_init
  ok=.true.
  nvecs=0
  iqual=0
  if(getflag('-help')) then
    call helptext(STDERR)
    call exit(0)
  endif
  if(.not.(mygetarg('-i',mdfile).or.mygetarg('-a',mdfile))) then
    write(STDERR,*) 'Specify md file with -i or album with -a'
    ok=.false.
  endif
  if(.not.mygetarg('-o',outfile)) then
    write(STDERR,*) 'Specify compressed file with -o'
    ok=.false.
  endif
  if(.not.mygetarg('-n',natom)) then
    write(STDERR,*) 'Specify number of atoms with -n'
    ok=.false.
  endif
  if(.not.ok) call exit(1)
  verbose=getflag('-v')
  lowmem=getflag('-lowmem')
  natom3=3*natom
!
! check output file doesn't already exist...
!
  inquire(file=outfile,exist=exists)
  if(exists) then
    write(STDERR,*) trim(outfile),': already exists'
    call exit(1)
  endif
!
!  are we chosing set no. of eigenvectors or quality measure?
!
  ok=mygetarg('-e',nvecs)
  ok=mygetarg('-q',iqual)
  if(iqual.eq.0.and.nvecs.eq.0) iqual=IDQ
  if(iqual.ne.0) nvecs=0
  if(iqual.ge.100) then
    write(STDERR,*) 'Error: -q must be in range 1-99'
    call exit(1)
  endif
  album=getflag('-a')
!
! open input file, read all in, masking down if asked
!
  allocate(xt(natom3),mask(natom3),STAT=stat)
  if(stat /= 0) then
    write(STDERR,*) 'Memory allocation error'
    call exit(1)
  endif
  if(album) then
    xalbf=xalbopen(8,mdfile,natom)
    if(xalbf%status /= 0) then
      write(STDERR,*) 'Error opening xalbumfile, xalbopen returns ',xalbf%status
      call exit(1)
    end if
    nframes=xalbf%nframes
    title=xalbf%title
    if(verbose) then
      write(STDERR,*) 'Album file contains:'
      i=1
      j=1
      do while(xalbf%lastframe(j) > 0)
        write(STDERR,*) trim(xalbf%track(j)%filename),': frames ',i, &
          '- ',xalbf%lastframe(j)
        i=xalbf%lastframe(j)+1
        j=j+1
      end do
    endif
  else
    xf=xopen(mdfile,natom)
    if(xf%status /= 0) then
      write(STDERR,*) 'Error opening mdfile, xopen returns ',xf%status
      call exit(1)
    end if
    nframes=xf%nframes
    title=xf%title
  endif
  if(mygetarg('-mask', maskfile).or.mygetarg('-freeze',maskfile)) then
    inquire(file=maskfile,exist=exists)
    if(.not.exists) then
      write(STDERR,*) trim(maskfile),': no such file?'
      call exit(1)
    endif
    open(9,file=maskfile,status='old')
    call maskread(9,natom,nmask,mask)
    if(verbose) write(STDERR,*) 'Masking down to ',nmask,' atoms'
    natom3=nmask*3
    natom=nmask
    if(getflag('-freeze')) then
      rewind(9)
      allocate(xtmp(3*natom),freezeidx(natom))
      call freezeread(9,xtmp,freezeidx,natom,nfreeze)
      if(verbose) write(STDERR,*) 'Freezing down to ',nfreeze,' groups'
    endif
    close(9)
  else
    mask=.true.
  endif
  if(album) then
    allocate(x(natom3,xalbf%nframes),STAT=stat)
  else
    allocate(x(natom3,xf%nframes),STAT=stat)
  endif
  if(stat /= 0) then
    write(STDERR,*) 'Memory allocation error'
    call exit(1)
  endif
  do i=1,nframes
    if(album) then
      xt=xalbsnap(xalbf,i)
      if(xalbf%status /= 0) then
        write(STDERR,*) 'Error in album file at snapshot ',i
        call exit(1)
      endif
    else
      xt=xsnap(xf,i)
    endif
    x(:,i)=pack(xt,mask)
    if(getflag('-freeze')) call xfreeze(x(:,i),xtmp,freezeidx,natom,nfreeze)
  end do
  if(getflag('-freeze')) deallocate(xtmp)
  if(verbose) print*,nframes,' frames read'
!
! sanity check on number of eigenvectors requested
!
  if(nvecs.gt.nframes) then
    write(STDERR,*) 'Warning - will only use ',nframes-1,' eigenvectors'
    nvecs=nframes-1
  end if
!
! bail out if no vectors possible...
!
  if(nframes<2) then
    write(STDERR,*) 'Error - insufficient frames to proceed'
    call exit(1)
  end if
!
!  fitting process: 1. fit to 1st frame, 2: calc avg, 3: fit to this

  if(getcarg('-mass',massfile)) then
    allocate(mass(natom*3))
    open(9,file=massfile,status='old')
    mass=1.0
    read(9,*) (mass(i),i=1,natom3,3)
    close(9)
  endif
!
  if(.not.getflag('-nofit')) then
    if(getflag('-mass')) then
      call trajfit(x(:,2),x(:,1),natom,nframes-1,ierr,mass(::3))
      call trajavg(x,xt,natom,nframes,ierr)
      call trajfit(x,xt,natom,nframes,ierr,mass(::3))
!
! now prepare the masses for CV matrix calculation...
!
      do i=1,natom3,3
        mass(i)=sqrt(mass(i))
        mass(i+1)=mass(i)
        mass(i+2)=mass(i)
      end do
    else
      call trajfit(x(:,2),x(:,1),natom,nframes-1,ierr)
      call trajavg(x,xt,natom,nframes,ierr)
      call trajfit(x,xt,natom,nframes,ierr)
    endif
  else
    call trajavg(x,xt,natom,nframes,ierr)
  endif
!
! remove average structure
!
  do i=1,nframes
    x(:,i)=x(:,i)-xt
  end do
!
! see if we are to use Ian's fastpca method...
!
  dofast=.not.(getflag('-nofast').or.getflag('-mass'))
!
! the value '5' here is pretty arbitrary - should be optimised...
!
  dofast=(dofast.and.nframes*5<natom3)
  if(dofast) then
    if(verbose) write(*,*) 'Using fastpca method'
    nf2=nframes
    allocate(w(natom3),z(natom3,nf2),STAT=stat)
    if(stat /= 0) then
      write(STDERR,*) 'Memory allocation error'
      call exit(1)
    endif
    call fastpca(x,w,z,natom3,nf2)
  else
!
!  calculate covariance matrix
!
    if(verbose) write(*,*) 'Calculating covariance matrix'
    rframes=1.0/nframes
    allocate(cm(natom3,natom3),STAT=stat)
    if(stat /= 0) then
      write(STDERR,*) 'Memory allocation error'
      call exit(1)
    endif
    call &
! matmul is SLOOOOW - use BLAS instead...
!    cm=matmul(x,transpose(x))*rframes
    sgemm("N","T",natom3,natom3,nframes,1.0,x,natom3,x,natom3,0.0,cm,natom3)
    cm=cm*rframes
!
! do mass weighting if asked...
!
    if(getflag('-mass')) then
      do i=1,natom3
        do j=1,natom3
          cm(i,j)=cm(i,j)*mass(i)*mass(j)
        end do
      end do
    endif
  end if
!
! if operating in 'lowmem' mode, save x to a scratch file and
! deallocate the array
!
  if(lowmem) then
    open(31,form='unformatted',status='scratch')
    do i=1,nframes
      write(31) x(:,i)
    end do
    rewind(31)
    deallocate(x)
  endif
!
!  if nvecs = 0, we need to calc all eigenvectors to work out how
!  many are needed to reach specified quality limit
!
  vtot=0.0
  if(nvecs.eq.0) then
    if(verbose) print*, 'Calculating vectors required'
    if(.not.dofast) then
      allocate(w(natom3),STAT=stat)
      if(stat /= 0) then
        write(STDERR,*) 'Memory allocation error'
        call exit(1)
      endif
      w=evals(cm)
!
!  calc now many eigenvectors will be needed
!
      do i=1,natom3
        vtot=vtot+w(i)
      end do
      vsum=0.0
      do while (int(vsum*100/vtot).lt.iqual)
        nvecs=nvecs+1
        vsum=vsum+w(nvecs)
      end do
    else
      do i=1,nf2
        vtot=vtot+w(i)
      end do
      vsum=0.0
      do while(int(vsum*100/vtot).lt.iqual)
        nvecs=nvecs+1
        vsum=vsum+w(nvecs)
      end do
    endif
    if(verbose) print*,nvecs,' eigenvectors will be used'
    if(verbose) then
      print*,int(vsum*100/vtot),'% of variance captured'
    endif
  else
!
!  number of eigenvectors specified at start
!
    if(verbose) print*,nvecs,' eigenvectors will be used'
    allocate(w(natom3))
    w=evals(cm)
  endif
  if(.not.dofast) then
    nf2=nvecs
    allocate(z(natom3,nvecs),STAT=stat)
    if(stat /= 0) then
      write(STDERR,*) 'Memory allocation error'
      call exit(1)
    endif
!
! now the real thing...
!
    if(verbose) print*,'Calculating eigenvectors'
    z=evecs(cm,nvecs)
!
! don't need that big covariance matrix any more...
!
    deallocate(cm)
  endif
!
! now we can do the projections
!
  if(verbose) print*,'Calculating projections'
  allocate(pn(nf2,nframes),STAT=stat)
  if(stat /= 0) then
    write(STDERR,*) 'Memory allocation error'
    call exit(1)
  endif
  if(lowmem) then
    allocate(xtmp(natom3))
    do i=1,nframes
      read(31) xtmp
      pn(:,i)=matmul(transpose(z),xtmp)
    end do
    close(31)
    deallocate(xtmp)
  else
!    pn=matmul(transpose(z),x)
     call sgemm("T","N",nf2,nframes,natom3,1.0,z,natom3,x,natom3,0.0,pn,nf2)
  endif
!
!  now write the compressed output file
!
  if(verbose) print*,'Writing compressed file'
  binary=.not.getflag('-formatted')
  if(.not.binary) then
    open(9,file=outfile,status='new')
    write(9,'(a4)') "PCZ0"
    write(9,'(a80)') title
    write(9,'(3i8,f8.3)') natom,nframes,nvecs,vtot
    write(9,'(10f8.3)') xt
    do k=1,nvecs
      write(9,'(10f8.3)') (z(j,k),j=1,natom3)
      write(9,'(10f8.3)') w(k)
      write(9,'(10f8.3)') (pn(k,j),j=1,nframes) 
    end do
  else
    allocate(evec(natom3),proj(nframes),STAT=stat)
    if(stat /= 0) then
      write(STDERR,*) 'Memory allocation error'
      call exit(1)
    endif
    pczf=pczcreate(9,outfile,title,natom,nframes,nvecs,real(iqual),xt)
    if(pczf%status /= 0) then
      write(0,*) 'Error creating output file: pczcreate returns ',pczf%status
      call exit(1)
    end if
    do k=1,nvecs
      evec=z(:,k)
      eval=w(k)
      proj=pn(k,:)
      stat=pczwrite(pczf,k,evec,eval,proj)
    end do
!
! default to the new PCZ6 format
!
    pczf%version='PCZ6'
    stat=pczclose(pczf)
  endif
  close(9)
end program pcazip

subroutine maskread(iunit,natom,nmask,mask)
implicit none
!
!  read a mask file, in pdb format
!
  logical        :: mask(natom*3)
  character(80)  :: line
  integer        :: iunit,natom,nmask,stat,i

  nmask=0
  mask=.false.
  do
    read(iunit,'(a80)',iostat=stat) line
    if(stat /= 0) exit
    if(line(1:4).eq.'ATOM') then
      nmask=nmask+1
      read(line,'(6x,i5)') i
      mask(3*i-2)=.true.
      mask(3*i-1)=.true.
      mask(3*i)=.true.
    endif
  end do
end subroutine maskread

subroutine freezeread(iunit,x,idx,nx,nidx)
use geom
implicit none

  integer          :: iunit,nx,nidx
  real             :: x(3*nx)
  integer          :: idx(nx)
  character(80)    :: line
  integer          :: i,j,stat

  i=-2
  nidx=0
  do
    read(iunit,'(a80)',iostat=stat) line
    if(stat /= 0) exit
    if(line(1:4).eq.'ATOM') then
      i=i+3
      read(line,'(30x,3f8.3)') x(i:i+2)
    else if(line(1:4).eq.'TER ') then
      nidx=nidx+1
      idx(nidx)=i+2
    endif
  end do
  if(idx(nidx)<i+2) then
    nidx=nidx+1
    idx(nidx)=i+2
  endif

end subroutine freezeread


subroutine xfreeze(x,xref,idx,natom,nfreeze)
use geom
implicit none

  integer   :: natom,nfreeze
  real      :: x(3*natom),xref(3*natom)
  integer   :: idx(nfreeze)

  integer   :: i,j,k,l
  real      :: xtmp(3*natom)

  do i=1,nfreeze
    if(i>1) then
      j=idx(i-1)+1
    else
      j=1
    endif
    k=idx(i)
    l=(k-j+1)/3
    xtmp(j:k)=lsqfit(x(j:k),xref(j:k),l)
  end do
  x=xtmp

end subroutine xfreeze

subroutine fastpca(x,w,z,p,n)
use eigenutils
implicit none

  integer          :: p,n
  real             :: x(p,n)
  real             :: w(n),z(p,n)
  real             :: cm(n,n),lz(n,n),lw(n),vv(n),r
  integer          :: od(n),ot(1),i
  logical          :: m(n)

  r=1.0/n
  cm=matmul(transpose(x),x)*r
  lw=evals(cm)
  lz=evecs(cm,n)

  z=matmul(x,lz)

  do i=1,n
    vv(i)=sqrt(sum(z(:,i)*z(:,i)))
    z(:,i)=z(:,i)/vv(i)
  end do

  w=sqrt(abs(lw*r))*vv

  m=.true.
  do i=1,n
    ot=maxloc(w,m)
    od(i)=ot(1)
    m(ot(1))=.false.
  end do

  w=w(od)
  z=z(:,od)

end subroutine fastpca
    
subroutine helptext(iunit)
implicit none
!
!  print help text
!
  integer   :: iunit
  write(iunit,"('*** pcazip version 4.1 ***')")
  write(iunit,"('Usage:'/'pcazip {-i|-a} infile -o outfile -n natoms')")
  write(iunit,"('     [-v] [-mask maskfile] [-e nev] [-q qual] [-formatted]')")
  write(iunit,"('Details:'/'-i infile    | any VMD-supported trajectory file')")
  write(iunit,"('OR -a infile | album file containing list of traj. files')")
  write(iunit,"('-o outfile   | compressed output file')")
  write(iunit,"('-n natoms    | number of atoms in a snapshot')")
  write(iunit,"('-v           | verbose diagnostics')")
  write(iunit,"('-nofit       | bypass least-squares fitting step')")
  write(iunit,"('-nofast      | do not use fastpca method, even if the ',/, &
    &           '             | number of snapshots is less than 1/15',/, &
    &           '             | the number of atoms')")
  write(iunit,"('-formatted   | human readable format (for checking)')") 
  write(iunit,"('-mask mfile  | pdb format mfile - only atoms in this')")
  write(iunit,"('             | will be included in the compression. Atom')")
  write(iunit,"('             | numbers (2nd field) must be correct!')")
  write(iunit,"('-e nev       | Include nev eigenvectors.', &
    &' Overrides -q option')")
  write(iunit,"('-q qual      | Include enough eigenvectors', &
    &' to capture qual%')")
  write(iunit,"('             | (1<=qual<=99) of the total variance')")
  write(iunit,"('Note: in the absence of any -q or -e option,', &
     &' default is to ')")
  write(iunit,"('      include enough eigenevectors to capture', &
    &' 90% of variance')")
end subroutine helptext

subroutine die(unit,msg)
implicit none
  integer            :: unit
  character(len=*)   :: msg

  write(unit,*) msg
  call exit(1)
end subroutine die
