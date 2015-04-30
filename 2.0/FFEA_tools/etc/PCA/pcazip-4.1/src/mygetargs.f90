! Mygetargs: general argument processing routines
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

module mygetargs
!
! can now use 'mygetarg' as a generic function in place
! of getiarg, getrarg, or getcarg
!
interface mygetarg
  module procedure getiarg, getrarg, getcarg
end interface

contains

  function getiarg(flag,ivalue)
    implicit none
    character(*)  :: flag
    character(160):: arg
    integer       :: ivalue
    logical       :: getiarg

    integer       :: l,i,j

    getiarg=.false.
    i=0
    do j=1,iargc()-1
      call getarg(j,arg)
      if(arg.eq.flag) i=j
    end do
    if(i.gt.0) then
      call getarg(i+1,arg)
      read(arg,*) ivalue
      getiarg=.true.
    endif
    end function getiarg

  function getrarg(flag,value)
    implicit none
    character(*)  :: flag
    character(160):: arg
    real          :: value
    logical       :: getrarg

    integer       :: l,i,j

    getrarg=.false.
    i=0
    do j=1,iargc()-1
      call getarg(j,arg)
      if(arg.eq.flag) i=j
    end do
    if(i.gt.0) then
      call getarg(i+1,arg)
      read(arg,*) value
      getrarg=.true.
    endif
    end function getrarg

    function getcarg(flag,cvalue)
    implicit none
    character(*)  :: flag,cvalue
    character(160) :: arg
    logical       :: getcarg

    integer       :: l,i,j

    getcarg=.false.
    i=0
    do j=1,iargc()-1
      call getarg(j,arg)
      if(arg.eq.flag) i=j
    end do
    if(i.gt.0) then
      call getarg(i+1,cvalue)
      getcarg=.true.
    endif
    end function getcarg

    function getflag(flag)
    implicit none
    character(*)  :: flag
    character(160) :: arg
    logical       :: getflag

    integer       :: l,i,j

    getflag=.false.
    i=0
    do j=1,iargc()
      call getarg(j,arg)
      if(arg.eq.flag) i=j
    end do
    if(i.gt.0) getflag=.true.
    end function getflag
end module mygetargs
