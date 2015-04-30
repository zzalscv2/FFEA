! Pcaunzip: uncompress MD trajectory files compressed with pcazip
! version 4.1
!
! This version with enhanced decompression rate and DCD output support
! written by Yiming Chen.
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
program pcaunzip
use pcz_io
use mygetargs
implicit none

  character(160)    :: infile,outfile,fort
  integer,parameter :: STDERR=0, STDOUT=6
  integer           :: ounit,i,stat,iv1,iv2
  type(pczfile)     :: pczf

  
  integer           :: ICNTRL(20),NATOM
  integer           :: NTITLE
  character(4)      :: HDRC
  character(80)     :: TITLE
  parameter (HDRC='CORD')
  real(4)           :: DELTA
  real(4),allocatable     :: snapx(:),snapy(:),snapz(:)

  if(.not.getcarg('-i',infile)) then
    write(STDERR,*) 'Error - specify pczfile with -i'
    call exit(1)
  end if

  pczf=pczopen(1,infile)
  if(pczf%status /= 0) then
    write(STDERR,*) 'Error with pczfile: pczopen returns ',pczf%status
    call exit(1)
  end if

  if(getcarg('-o',outfile)) then
    ounit=2
  else
    ounit=STDOUT
  end if


  if(.not.getiarg('-iv1',iv1)) then
    iv1=1
    iv2=pczf%nvecs
  else
    if(.not.getiarg('-iv2',iv2)) iv2=iv1
  endif
  if(iv1>pczf%nvecs) then
    write(STDERR,*) &
      'Error - iv1 greater than number of eigenvectors available'
    call exit(1)
  endif
  if(iv2>pczf%nvecs) then
    write(STDERR,*) &
      'Warning - resetting iv2 to number of eigenvectors available'
    iv2=pczf%nvecs
  endif

! The new added option, which allows users to select the format of the
! output file (either amber ascii format or charmm dcd format)

  if(.not.getflag('-format')) then
    fort="amber"
  else
    if(.not.getcarg('-format',fort)) call exit(1)
  end if

  if (fort=='amber') then
    if(ounit/=STDOUT) open(ounit,file=outfile,status='new')
    write(ounit,'(a80)') pczf%title
    stat=0
    do i=1,pczf%nframes
      if(stat==0) write(ounit,'(10f8.3)',iostat=stat) pczsnap(pczf,i,iv1,iv2)
    end do
  else if (fort=='binpos') then
    if(ounit/=STDOUT) open(ounit,file=outfile,status='new',&
      form='unformatted',access='stream')
    write(ounit) 'fxyz'
    stat=0
    do i=1,pczf%nframes
      if(stat==0) write(ounit,iostat=stat) pczf%natoms,pczsnap(pczf,i,iv1,iv2)
    end do
  else if (fort=='charmm') then
    if(ounit/=STDOUT) open(ounit,file=outfile,status='new',form='unformatted')
    allocate (snapx(3*pczf%natoms))
!CHARMM DCD file header section setup    
    do I=1,20
      ICNTRL(I)=0
    end do
    DELTA=0.001
    ICNTRL(1)=pczf%nframes
    ICNTRL(2)=1
    ICNTRL(3)=1
    ICNTRL(4)=1
    ICNTRL(8)=pczf%natoms*3-6
    ICNTRL(20)=34
    NTITLE=1
    TITLE='         Created by PCAZIP             ' 
    NATOM=pczf%natoms
    call int2singl(ICNTRL(10),DELTA)
    write(ounit) HDRC,ICNTRL  
    write(ounit) NTITLE,TITLE  
    write(ounit) NATOM  
! writing frames    
    do i=1,pczf%nframes
      snapx=pczsnap(pczf,i,iv1,iv2)
      write(ounit,iostat=stat) snapx(1:3*pczf%natoms-2:3)
      write(ounit,iostat=stat) snapx(2:3*pczf%natoms-1:3)
      write(ounit,iostat=stat) snapx(3:3*pczf%natoms:3)
    end do
    deallocate (snapx)
  else
    write(STDERR,*) 'Error - specified format is not supported'
  end if  
  stat=pczclose(pczf)
end program pcaunzip
  
subroutine int2singl(dest,source)

      IMPLICIT NONE
      INTEGER DEST,SOURCE
      DEST=SOURCE
      RETURN
end subroutine int2singl

subroutine die(unit,msg)
implicit none
  integer            :: unit
  character(len=*)   :: msg

  write(unit,*) msg
  call exit(1)
end subroutine die
