! molfile_substitute: Basic substitute for VMD molfile plugins
! Copyright (C) Charlie Laughton 2011
! A very basic emergency substitute where no proper molfile plugins
! are available. Just supports AMBER crd/crdbox and binpos formats
! plus (some?) dcd formats.
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
subroutine f77_molfile_init
implicit none

  integer, parameter :: MAXHANDLES=50
  integer            :: handles(MAXHANDLES),types(MAXHANDLES)

  common/handcom/handles,types

  handles=0
  types=0

  return
end subroutine f77_molfile_init

subroutine f77_molfile_open_read(handle,natoms,fname,type)
implicit none

  integer      :: handle,natoms
  character(*) :: fname, type
  
  integer            :: i,nat
  real               :: dum
  character(4)       :: magic
  character(80)      :: line
  logical            :: exists

  integer, parameter :: MAXHANDLES=50
  integer            :: handles(MAXHANDLES),types(MAXHANDLES)
  common/handcom/handles,types

  i=1
  do while(i<=MAXHANDLES.and.handles(i)>0)
    i=i+1
  end do
  if(i>MAXHANDLES) then
    write(0,*) 'Error - too many file handles needed'
    call exit(1)
  endif
!
! assume unit numbers above 50 are safe...
!
  handles(i)=i+50
  handle=handles(i)

  inquire(file=fname,exist=exists)
  if(.not.exists) then
    write(0,*) trim(fname),': no such file?'
    call exit(1)
  endif

  if(type=="auto") then
!
! look at the file extension
!
    i=index(fname,".")
    if(i>0) then
      select case(fname(i+1:len_trim(fname)))
      case("crd")
        types(handle-50)=1
        write(0,*) 'file type is crd'
      case("crdbox")
        types(handle-50)=2
        write(0,*) 'file type is crdbox'
      case("binpos")
        types(handle-50)=3
        write(0,*) 'file type is binpos'
!
! note type=4 is reserved for binpos with box...
!
      case("dcd")
        types(handle-50)=5
        write(0,*) 'file type is dcd'
      case default
        types(handle-50)=0
      end select
    endif
  else if(type=="crd") then
    types(handle-50)=1
  else if(type=="crdbox") then
    types(handle-50)=2
  else if(type=="binpos") then
    types(handle-50)=3
  else if(type=="dcd") then
    types(handle-50)=5
  endif

  if(types(handle-50)==0) then
    handles(handle-50)=0
    handle=-1
    write(0,*) 'Unrecognised file type'
    return
  else if(types(handle-50)==3) then
    write(0,*) 'Attempting to open binpos file'
    open(handle,file=fname,status='old',form='unformatted',access='stream')
    read(handle) magic
    if(magic/='fxyz') then
      write(0,*) 'Error - bad binpos file'
      close(handle)
      call exit(1)
    endif
    read(handle) nat
    if(nat/=natoms) then
      write(0,*) 'Warning: binpos file resetting natoms to ',nat
      natoms=nat
    endif
    do i=1,3*natoms
      read(handle) dum
    end do
    read(handle) nat
    if(nat/=natoms) then
!
! assume box is present
!
      types(handle-50)=4
      write(0,*) 'binpos file has box'
    endif
    rewind(handle)
    read(handle) magic
  else if(types(handle-50)==5) then
    write(0,*) 'Attempting to open dcd file'
    call dcdopen(handle,fname,nat)
    if(nat/=natoms) then
      write(0,*) 'Warning: dcd file resetting natoms to ',nat
      natoms=nat
    endif
  else
    write(0,*) 'Attempting to open crd or crdbox file'
    open(handle,file=fname,status='old')
    read(handle, '(a80)') line
    write(0,*) 'Hello'
    natoms=0
  end if
end subroutine f77_molfile_open_read

subroutine f77_molfile_read_next(handle,natoms,x,box,status)
implicit none

  integer      :: handle,natoms,status
  real         :: x(3*natoms),box(6)
  
  integer            :: i

  integer, parameter :: MAXHANDLES=50
  integer            :: handles(MAXHANDLES),types(MAXHANDLES)

  common/handcom/handles,types

  if(types(handle-50)==1) then
    read(handle,'(10f8.3)',iostat=status) x
  else if(types(handle-50)==2) then
    read(handle,'(10f8.3)',iostat=status) x
    read(handle,'(6f8.3)',iostat=status) box
  else if(types(handle-50)==3) then
    read(handle,iostat=status) i,x
  else if(types(handle-50)==4) then
    read(handle,iostat=status) i,x,box
  else if(types(handle-50)>4) then
    call dcdreadnext(handle,x,natoms,box,status)
  endif

  if(status==0) then
    status=1
  else
    status=0
  endif
end subroutine f77_molfile_read_next

subroutine f77_molfile_close_read(handle,status)
implicit none

  integer    :: handle, status

  integer, parameter :: MAXHANDLES=50
  integer            :: handles(MAXHANDLES),types(MAXHANDLES)

  common/handcom/handles,types

  if(handle>50) then
    close(handle)
    handles(handle-50)=0
    types(handle-50)=0
  endif
  status=0
end subroutine f77_molfile_close_read

  subroutine  dcdopen(unit,filename,natoms)
  implicit none

  integer      :: unit,natoms
  character(*) :: filename

  integer      :: status
!
! try opening the file...
!
  open(unit,file=filename,status='old',form='unformatted', access='stream', &
       iostat=status)

  if(status /= 0) then
    write(0,*) 'Error - can''t open dcd file'
    close(unit)
    unit=-1
    return
  endif
  call readdcdheader(unit,natoms)
  return
  end subroutine dcdopen

  subroutine readdcdheader(unit,natoms)
  implicit none
  integer      :: unit,natoms
  character*80 :: title,title2

  integer      :: irecl,idumm(10),nfixed,extrablock,i,nframes,status
  integer*8    :: lrecl
  real         :: delta
  double precision :: ddelta
  character(4) :: magic
  logical      :: ischarmm,longrecl,hasbox

  integer, parameter :: MAXHANDLES=50
  integer            :: handles(MAXHANDLES),types(MAXHANDLES)

  common/handcom/handles,types

  status=0
!
! see if we have an integer*4 or *8 record length indicator...
!
  read(unit) irecl
  if(irecl==84) then
    longrecl=.false.
  else
    rewind(unit)
    read(unit) lrecl
    if(lrecl==84) then
      longrecl=.true.
    else
!
! if we get to here then this is most likely a wrong-endian format file,
! for now, we can't deal with this - maybe next time...
!
      status=-2
      write(0,*) 'Problem with dcd file - wrong-endian?'
    endif
  endif
  if(status/=0) then
    close(unit)
    unit=-1
    return
  endif
!
! check the magic...
!
  read(unit) magic
  if(magic /= "CORD") then
    write(0,*) 'Error - no magic word in dcd file?'
    close(unit)
    unit=-1
    return
  endif
 if(longrecl) types(unit-50)=6
!
! read the rest of the first line and see if this is a charmm format file
!
  read(unit) idumm(1:9),delta,idumm
  ischarmm=(idumm(10) /= 0)
!
! OK, back to the start, now can read the whole first record correctly
!
  rewind(unit)
  if(longrecl) then
    read(unit) lrecl
  else
    read(unit) irecl
  endif
  read(unit) magic,nframes,idumm(1:8)
!
! if we have fixed atoms, give up...
!
  if(idumm(8)/=0) then
    write(0,*) 'Error - fixed atoms in dcd file not supported'
    close(unit)
    unit=-1
    return
  endif
!
! see if there is going to be box info
!
  extrablock=0
  if(ischarmm) then
    read(unit) delta,extrablock
  else
    read(unit) ddelta
  endif
  hasbox=(extrablock/=0)
  if(hasbox) types(unit-50)=types(unit-50)+2
!
! now read the rest of the line
!
  read(unit) idumm(1:9)
  if(idumm(1)/=0) then
!
! extra dimension in this dcd file - no thanks
!
    write(0,*) 'Error - extra dimension in dcd file not supported'
    close(unit)
    unit=-1
    return
  endif
!
! read the trailing recordlength indicator
!
  if(longrecl) then
    read(unit) lrecl
  else
    read(unit) irecl
  endif
!
! OK, onto record two - the title(s)
!
  if(longrecl) then
    read(unit) lrecl
  else
    read(unit) irecl
  endif
  
  read(unit) idumm(1)
  read(unit) title
!
! we discard all but the first title
!
  do i=1,idumm(1)-1
    read(unit) title2
  end do
 
  if(longrecl) then
    read(unit) lrecl
  else
    read(unit) irecl
  endif
!
! next line: the number of atoms
!
  if(longrecl) then
    read(unit) lrecl
  else
    read(unit) irecl
  endif

  read(unit) natoms

  if(longrecl) then
    read(unit) lrecl
  else
    read(unit) irecl
  endif
!
! OK - from now on in its frames...
!
  return
  end subroutine readdcdheader

  subroutine dcdreadnext(handle,x,natoms,box,status)
  implicit none
  integer       :: handle,natoms,status
  real          :: x(3*natoms),box(6)

  integer, parameter :: MAXHANDLES=50
  integer            :: handles(MAXHANDLES),types(MAXHANDLES)

  double precision :: dbox(6)
  integer          :: irecl,i
  integer*8        :: lrecl
  logical          :: longrecl,hasbox

  common/handcom/handles,types

  status=0
  select case(types(handle-50))
    case(5)
      longrecl=.false.
      hasbox=.false.
    case(6)
      longrecl=.true.
      hasbox=.false.
    case(7)
      longrecl=.false.
      hasbox=.true.
    case(8)
      longrecl=.true.
      hasbox=.true.
    case default
      write(0,*) 'Error - unrecognised dcd file type'
      status=-2
      return
  end select

  if(hasbox) then
    if(longrecl) then
      read(handle,iostat=status) lrecl
    else
      read(handle,iostat=status) irecl
    endif
    if(status/=0) return
    read(handle) dbox
    box=dbox
    if(longrecl) then
      read(handle) lrecl
    else
      read(handle) irecl
    endif
  endif
  
  if(longrecl) then
    read(handle,iostat=status) lrecl
  else
    read(handle,iostat=status) irecl
  endif
  if(status/=0) return

  read(handle) (x(i),i=1,3*natoms-2,3)
  if(longrecl) then
    read(handle) lrecl,lrecl
  else
    read(handle) irecl,irecl
  endif
  read(handle) (x(i),i=2,3*natoms-1,3)
  if(longrecl) then
    read(handle) lrecl,lrecl
  else
    read(handle) irecl,irecl
  endif
  read(handle) (x(i),i=3,3*natoms,3)
  if(longrecl) then
    read(handle) lrecl
  else
    read(handle) irecl
  endif
  status=0
  return
  end subroutine dcdreadnext

