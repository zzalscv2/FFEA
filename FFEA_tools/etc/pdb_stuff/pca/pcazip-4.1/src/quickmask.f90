! Quickmask: utility to concatenate and/or strip down trajectory files.
!            Writes Scripps .binpos format output files.
! Version 4.0
! Copyright (C) Charlie Laughton 2010
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
program quickmask
use mygetargs
use x_io
implicit none

  integer,parameter             :: STDERR=0
  real, allocatable             :: x(:,:),xt(:)
  logical, allocatable          :: mask(:)
  character(160)                :: mdfile,outfile,maskfile
  character(80)                 :: title
  logical                       :: hasbox,ok,verbose,binary,album
  logical                       :: exists
  integer                       :: trajframes,i,j,k,i2,nf2
  integer                       :: natom,natom3,nframes,ierr
  integer                       :: nmask,ip1,ip2,stat
  real                          :: rframes
  type(xfile)                   :: xf
  type(xalbum)                  :: xalbf

  call f77_molfile_init
  ok=.true.
  if(getflag('-help')) then
    call helptext(STDERR)
    call exit(0)
  endif
  if(.not.(mygetarg('-i',mdfile).or.mygetarg('-a',mdfile))) then
    write(STDERR,*) 'Specify md file with -i or album with -a'
    ok=.false.
  endif
  if(.not.mygetarg('-o',outfile)) then
    write(STDERR,*) 'Specify new binpos file with -o'
    ok=.false.
  endif
  if(.not.mygetarg('-n',natom)) then
    write(STDERR,*) 'Specify number of atoms with -n'
    ok=.false.
  endif
  if(.not.ok) call exit(1)
  verbose=getflag('-v')
  natom3=3*natom
!
! check output file doesn't already exist...
!
  inquire(file=outfile,exist=exists)
  if(exists) then
    write(STDERR,*) trim(outfile),': already exists'
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
  if(mygetarg('-mask', maskfile)) then
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
    close(9)
  else
    do i=1,natom3
      mask(i)=.true.
    end do
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
  end do
  if(verbose) print*,nframes,' frames read'
!
! write out the new binpos file
!
  open(9,file=outfile,status='new',form='unformatted', access='stream')
  write(9) 'fxyz'
  do i=1,nframes
    write(9) natom,x(:,i)
  end do
  close(9)
  call exit(0)
end program quickmask

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

subroutine helptext(iunit)
implicit none
!
!  print help text
!
  integer   :: iunit
  write(iunit,"('*** quickmask version 4.0 ***')")
  write(iunit,"('Usage:'/'quickmask {-i|-a} infile -o outfile -n natoms')")
  write(iunit,"('     [-v] [-mask maskfile]')")
  write(iunit,"('Details:'/'-i infile    | amber format trajectory file')")
  write(iunit,"('OR -a infile | album file containing list of traj. files')")
  write(iunit,"('-o outfile   | compressed output file')")
  write(iunit,"('-n natoms    | number of atoms in a snapshot')")
  write(iunit,"('-v           | verbose diagnostics')")
  write(iunit,"('-mask mfile  | pdb format mfile - only atoms in this')")
  write(iunit,"('             | will be included in the compression. Atom')")
  write(iunit,"('             | numbers (2nd field) must be correct!')")
end subroutine helptext
