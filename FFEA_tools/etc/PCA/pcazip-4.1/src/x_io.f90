! x_io: General routines for I/O on MD trajectory files
! Copyright (C) Charlie Laughton 2011
! This version uses the VMD molfile plugins
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
module x_io
implicit none

  type :: xfile
    character(160):: filename
    character(80) :: title
    integer       :: handle
    integer       :: status
    integer       :: natoms
    integer       :: nframes
    integer       :: firstframe
    integer       :: lastframe
    integer       :: stepsize
    integer       :: nextframe
    logical       :: hasbox
    real          :: box(6)
    character*10  :: type
  end type xfile

  type :: xalbum
    character*160 :: filename
    character*80  :: title
    integer       :: unit
    integer       :: status
    integer       :: natoms
    integer       :: nframes
    integer       :: fidx
    type(xfile)   :: track(500)
    integer       :: lastframe(500)
    type(xfile)   :: xf
  end type xalbum
    
  contains
    function xalbopen(un,fname,natoms) result (xalb)
    implicit none
    character*160  :: fname  ! filename
    integer        :: un     ! unit number
    integer        :: natoms ! atoms
    type(xalbum)   :: xalb   ! xalbum structure returned

    character*160  :: line
    integer        :: ios,i
    logical        :: exists

!
! initialise xalbum structure
!
    xalb%filename=fname
    xalb%unit=un
    xalb%title=''
    xalb%status=0
!
! NB: natoms may be zero if we are hoping that the trajectory files
! themselves will provide this information later
!
    xalb%natoms=natoms
    xalb%nframes=0
    xalb%fidx=0
    xalb%lastframe=0
!
! open and extract filenames from album file
!
    inquire(file=fname,exist=exists)
    if(.not.exists) then
      write(0,*) trim(fname),': no such file?'
      call exit(1)
    endif
    open(un,file=fname,status='old',iostat=ios)
    xalb%status=ios
    if(ios /= 0) return
    i=0
    read(un,'(a160)',iostat=ios) line
    do while (ios==0)
      i=i+1
      xalb%fidx=i
!
! this next bit is hard-coded - improve one day...
!
      if(i>500) then
        write(0,*) 'Error - maximum number of 500 files in an album'
        call exit(1)
      endif
      xalb%track(i)=xopen(line,natoms)
      if(xalb%track(i)%status/=0) then
        write(0,*) 'Error with file ',trim(xalb%track(i)%filename)
        xalb%status=xalb%track(i)%status
        return
      end if
      if(i==1) then
        xalb%lastframe(i)=xalb%track(i)%nframes
      else
        xalb%lastframe(i)=xalb%track(i)%nframes+xalb%lastframe(i-1)
      end if
      xalb%nframes=xalb%nframes+xalb%track(i)%nframes
      read(un,'(a160)',iostat=ios) line
    end do
    close(un)
!
! reset everything to start of first track
!
    xalb%fidx=1
    xalb%xf=xalb%track(1)
    xalb%natoms=natoms
    xalb%status=xalb%xf%status
    end function xalbopen


    function xopen(fname,natoms) result(xf)
    implicit none
    character(160) :: fname  ! filename
    integer        :: natoms ! atoms - may be zero  on entry if unknown
    type(xfile)    :: xf     ! xfile structure returned

    character(80)  :: line
    real,allocatable  :: x(:)
    real           :: d12,d22,dtol
    integer        :: i,j,ierr,nat
    logical        :: exists
    character*10   :: intype

! initialise the xfile structure
    if(scan(fname,"(")>0) then
      i=scan(fname,"(")
      j=scan(fname,")")
      call xparse(fname(i+1:j-1),xf%firstframe,xf%lastframe,xf%stepsize,ierr)
      if(ierr /= 0) then
        write(0,*) trim(fname),': syntax error'
        call exit(1)
      endif
      xf%filename=fname(1:i-1)
    else
      xf%firstframe=1
      xf%lastframe=0
      xf%stepsize=1
      xf%filename=fname
    endif
    xf%type="auto"
    xf%title=''
    xf%status=0
    xf%natoms=natoms; xf%nframes=0; xf%nextframe=0
    xf%hasbox=.false.
    xf%box=0.0
!
! first peek at the file to check it out
!
    inquire(file=xf%filename,exist=exists)
    if(.not.exists) then
      write(0,*) trim(xf%filename),': no such file?'
      call exit(1)
    endif

    call f77_molfile_open_read(xf%handle,natoms,xf%filename,xf%type)
    if(xf%natoms==0.and.natoms>0) xf%natoms=natoms
!
!  if the file type can't be recognised automatically, try AMBER
!  - but only if we already know how many atoms per frame
!
    if(xf%handle<0.and.xf%natoms>0) then
!      call f77_molfile_close_read(xf%handle,xf%status)
      xf%type="crd"
      write(0,*) 'Trying AMBER...'
      call f77_molfile_open_read(xf%handle,natoms,xf%filename,xf%type)
!
! (beware that natoms will return as zero for an AMBER file)
!
!
! ok - lets try to find out if there is a box or not. To do this
! read in two snapshots and calculate the apparent separation between
! the first two atoms. If it is not more or less conserved we have
! probably read box dimensions in snapshot 2 instead of coordinates
!
      dtol=0.2
      allocate(x(3*xf%natoms))
      xf%status=1
      call f77_molfile_read_next(xf%handle,xf%natoms,x,xf%box,xf%status)
      if(xf%status==0) then
        write(0,*) 'Problem with trajectory file ',trim(xf%filename)
        call exit(1)
      endif
      d12=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
      xf%status=1
      call f77_molfile_read_next(xf%handle,xf%natoms,x,xf%box,xf%status)
      d22=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
      if(abs(d12-d22).gt.dtol) then
!
! hmm - distance looks dodgy - try crdbox instead
!
        call f77_molfile_close_read(xf%handle,xf%status)
        xf%type="crdbox"
        call f77_molfile_open_read(xf%handle,natoms,xf%filename,xf%type)
        xf%status=1
        call f77_molfile_read_next(xf%handle,xf%natoms,x,xf%box,xf%status)
        if(xf%status==0) then
          write(0,*) 'Problem with trajectory file ',trim(xf%filename)
          call exit(1)
        endif
        d12=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
        xf%status=1
        call f77_molfile_read_next(xf%handle,xf%natoms,x,xf%box,xf%status)
        d22=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
        if(abs(d12-d22).gt.dtol) then
          write(0,*) 'Error - unrecognised format file ',trim(xf%filename)
          !call exit(1)
        endif
      endif
    else if(xf%handle<0) then
      write(0,*) 'Error - need to know number of atoms per snapshot'
      call exit(1)
    endif
!
! if at this stage natoms is not zero, it's the correct value as read
! from some non-AMBER format trajectory file...
!
    xf%natoms=max(xf%natoms,natoms)
!
! ok, find number of frames
!
    xf%nframes=0
    xf%status=1
    if(.not.allocated(x)) allocate(x(3*xf%natoms))
    do while(xf%status>0)
      xf%status=0
      call f77_molfile_read_next(xf%handle,xf%natoms,x,xf%box,xf%status)
      xf%nframes=xf%nframes+1
    end do
    xf%nframes=xf%nframes-1
!
! remember that if it was an amber crd or crdbox file, we read a couple
! of frames to check for a box
!
    if(xf%type/="auto") xf%nframes=xf%nframes+2
!
! do some adjustments if firstframe,lastframe,stepsize are set
!
    if(xf%lastframe==0) xf%lastframe=xf%nframes
    if(xf%lastframe>xf%nframes) then
      xf%lastframe=xf%nframes
      write(0,*) "Warning: ",trim(xf%filename),' only contains ',xf%nframes,' frames'
    endif
    if(xf%firstframe>xf%nframes) then
      write(0,*) "Error: ",trim(xf%filename),' only contains ',xf%nframes,' frames'
      call exit(1)
    endif
    xf%nframes=1+(xf%lastframe-xf%firstframe)/xf%stepsize
!
! note we are now positioned at the end of the file - future reads
! will need to start with a rewind...
!
    xf%nextframe=xf%nframes+1
!
! make sure we return the right number of atoms in natoms
!
    natoms=max(natoms,xf%natoms)
    end function xopen

    function xalbsnap(xalb,is) result (x)
    implicit none

    type(xalbum) :: xalb
    integer      :: is,i,is2
    real         :: x(3*xalb%natoms)

    if(is > xalb%nframes) then
      xalb%status=-2
      return
    end if

    i=1
    do while (is > xalb%lastframe(i))
      i=i+1
    end do

    if(i /= xalb%fidx) then
!
! the snapshot we want is not in the currently associated track
!
      xalb%fidx=i
      xalb%xf=xalb%track(i)
      xalb%status=xalb%xf%status
    end if
    if(xalb%status /= 0) return
!
! convert snapshot number into right number for this file
!
    if(xalb%fidx==1) then
      is2=is
    else
      is2=is-xalb%lastframe(xalb%fidx-1)
    endif
    x=xsnap(xalb%xf,is2)
    xalb%status=xalb%xf%status
    end function xalbsnap

    function xsnap(xf,is) result(x)
    implicit none

    type(xfile) :: xf
    integer     :: is,i,isl,nat,j
    real        :: x(3*xf%natoms)
    character*4 :: magic
    character*80 :: line
    character*90 :: lline

    nat=xf%natoms
    if(is > xf%nframes) then
      xf%status=-2
      return
    end if

    isl=xf%firstframe+(is-1)*xf%stepsize
    if(isl < xf%nextframe) then
!
! close and reopen to do a rewind
!
      call f77_molfile_close_read(xf%handle,xf%status)
      call f77_molfile_open_read(xf%handle,nat,xf%filename,xf%type)
      xf%nextframe=1
    end if

    do i=xf%nextframe,isl-1
      xf%status=1
      call f77_molfile_read_next(xf%handle,xf%natoms,x,xf%box,xf%status)
    enddo
    xf%status=1
    call f77_molfile_read_next(xf%handle,xf%natoms,x,xf%box,xf%status)
    xf%nextframe=isl+1
    xf%status=0
  end function xsnap

  subroutine xparse(str,i,j,k,err)
  implicit none
  character(len=*) :: str
  integer          :: i,j,k,err
  integer          :: l1,l2,l3
  
    err=0
    i=1
    j=0
    k=1
    l1=index(str,":")
    l3=len_trim(str)
    l2=index(str(l1+1:l3),":")
    l2=l1+l2
  
!   check str for invalid characters
    err=verify(str,"0123456789: ")
    if (err /= 0) return
!   check str for length
    if(l3==0) return
!   check str for excessive number of colons
    if(l2>0.and.l2<l3.and.index(str(l2+1:l3),":")>0) then
      err=-2
      return
    endif
!   looks like a valid str, so begin to parse...
    if(l1==0) then
!   no colons - just a number, matches 1 frame
      read(str(1:l3),*) i
      j=i
      return
    else if(l1>1) then
!   a number before the 1st colon - set i
      read(str(1:l1-1),*) i
    endif
!   now move onto j...
!   nothing beyond the 1st colon - j and k take default values
    if(l3==l1) return
    if(l2==l1) then
!   only 1 colon - rest of str is j
      read(str(l1+1:l3),*) j
      return
    endif
!   two colons  - j is what's between them (if anything)
    if(l2>l1+1) then
      read(str(l1+1:l2-1),*) j
    end if
!   now move onto k...
!   nothing beyond second colon
    if(l3==l2) return
!   there is a k...
    read(str(l2+1:l3),*) k
    return
  end subroutine xparse

  subroutine xclose(xf)
    implicit none

    type(xfile) :: xf

    call f77_molfile_close_read(xf%handle,xf%status)
    return
  end subroutine xclose

end module x_io
