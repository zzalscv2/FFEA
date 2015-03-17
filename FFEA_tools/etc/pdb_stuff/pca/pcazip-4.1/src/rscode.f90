! rscode: base-64 style encoding/decoding for .pcz files
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
module rscode
!
! functions that convert between reals and a base-64 style
! 5-character word. The format supports reals in the range
! -0.9999999e-7 to 0.9999999e7 - quite good enough for
! pcazip. Underflows are set to zero, overflows result in an
! exit
! The real number is encoded into 30 bits. From highest to lowest
! these are:
! 30-7: mantissa
! 6:    sign of mantissa (1==negative)
! 5-2:  exponent
! 1:    sign of exponent (1==negative)
!
  implicit none

contains
  function rtos(r) result(w)
  real         :: r
  character*5  :: w
!
! convert real number into 5-byte base-64 type word
!
  integer      :: m,e,b,i,j
  character*14 :: s
  character*64 :: code
!
! standard MIME base-64 encoding:
  code="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"

  write(s,'(e14.7)') r
  read(s,'(3x,i7,1x,i3)') m,e
  if(e<-7) then
   m=0
   e=0
  end if
  if(e>7) then
    write(0,*) 'rtos: overflow error: ',r,' too big'
    call exit(1)
  endif
  m=m*32
  if(r<0) m=m+16
  e=e*2
  if(e<0) e=-e+1
  b=m+e
  do j=5,1,-1
    i=iand(b,63)+1
    w(j:j)=code(i:i)
    b=ishft(b,-6)
  end do
end function rtos

function stor(w) result(r)
  real         :: r
  character*5  :: w
  integer      :: b,i,j,m,e
  character*64 :: code

  code="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"

  b=0
  do j=1,5
    i=index(code,w(j:j))-1
    b=b*64+i
  end do
  e=iand(b,15)
  if(iand(e,1)==1) e=-e
  e=e/2
  m=b/16
  if(iand(m,1)==1) m=-m
  m=m/2
  r=(10.0**(e-7))*real(m)
end function stor

end module rscode
