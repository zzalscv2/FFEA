! Trajutils: general utilities that operate on MD trajectory streams
! Copyright (C) Charlie Laughton 2006
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
module trajutils
!
!  a number of basic trajectory file manipulation utilities
!
!  subroutine trajfit(x,ref,natoms,nframes,ierr)
!    fits the frames in x to coordinates in ref
!    input: x(3,natoms,nframes),ref(3,natoms)
!    output: fitted coordinates in x
!
!  subroutine trajavg(x,xavg,natoms,nframes,ierr)
!    calculates average coordinates
!  
  contains
  
  subroutine trajfit(x,ref,natoms,nframes,ierr,w)
  implicit none
!
!  least-squares fits the frames in x to the structure ref
!  with optional argumnet w, does mass-weighted fitting
!
  real    ::  x(3,natoms,nframes),ref(3,natoms)
  real, optional :: w(natoms)
  real    ::  r(3,3),v(3)
  logical ::  findmove
  integer ::  natoms,nframes,ierr,i,j
  real    ::  xt,yt,zt,rmsd

  findmove=.true.
  do i=1,nframes
    if(present(w)) then
      call matfitw(natoms,ref,x(1,1,i),r,v,rmsd,findmove,w)
    else
      call matfit(natoms,ref,x(1,1,i),r,v,rmsd,findmove)
    endif
    do j=1,natoms
      xt=r(1,1)*x(1,j,i)+r(1,2)*x(2,j,i)+r(1,3)*x(3,j,i)+v(1)
      yt=r(2,1)*x(1,j,i)+r(2,2)*x(2,j,i)+r(2,3)*x(3,j,i)+v(2)
      zt=r(3,1)*x(1,j,i)+r(3,2)*x(2,j,i)+r(3,3)*x(3,j,i)+v(3)
      x(1,j,i)=xt
      x(2,j,i)=yt
      x(3,j,i)=zt
    end do
  end do

  end subroutine trajfit

  subroutine trajavg(x,xavg,natoms,nframes,ierr)
  implicit none
!
!  average coordinates in the trajectory
!
  real     ::  x(3,natoms,nframes), xavg(3,natoms)
  integer  ::  natoms,nframes,ierr,i,j
  real     ::  rframes

  rframes=1.0/nframes
  xavg=0.0
  do i=1,nframes
    do j=1,natoms
      xavg(1,j)=xavg(1,j)+x(1,j,i)
      xavg(2,j)=xavg(2,j)+x(2,j,i)
      xavg(3,j)=xavg(3,j)+x(3,j,i)
    end do
  end do
  xavg=xavg*rframes
  end subroutine trajavg
end module trajutils
