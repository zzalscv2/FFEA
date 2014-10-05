! eigenutils: general eigenvector and eigenvalue routines
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
module eigenutils
  implicit none
contains
  function evals(a) result(w)
!
! returns the eigenvalues of a
!
  real                           :: a(:,:)
  real                           :: w(size(a,1))

  integer, allocatable          :: isuppz(:), iwork(:)
  real             , allocatable :: alocal(:,:)
  real             , allocatable :: z(:),work(:)
  integer                       :: lwork,liwork
  integer                       :: il,iu,info,m
  real                           :: abstol,vl,vu
  real                           :: wtmp(size(a,1))
!
!  workspace calculation
!
  vl=0.0
  vu=0.0
  il=1
  iu=size(a,1)
  abstol=0.0
  lwork=-1
  liwork=-1
  allocate(work(1),iwork(1),isuppz(2*iu),z(1))
  call ssyevr('V','I','U',iu,a,iu,vl,vu,il,iu,abstol, &
     m,w,z,iu,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    print*,'SSYEVR returns with info= ',info
    call exit(1)
  endif
  lwork=int(work(1))
  liwork=iwork(1)
!
! now calculate all eigenvalues
!
  deallocate(work,iwork)
  allocate(work(lwork),iwork(liwork),alocal(iu,iu))
!
! copy a into alocal as it will be destroyed
! in the call to dsyevr (naughty!)
!
  alocal = a
  call ssyevr('N','A','U',iu,alocal,iu,vl,vu,il,iu,abstol, &
       m,wtmp,z,iu,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    print*,'SSYEVR returns with info= ',info
    call exit(1)
  endif
  deallocate(work,iwork,alocal)
!
! re-order the eigenvectors into descending order
!
  iu=size(a,1)
  do il=1,iu
    w(il)=wtmp(iu-il+1)
  end do
  end function evals

  function evecs(a,n) result(z)
!
! calculate top n eigenvectors of a
!
  real                           :: a(:,:)
  integer                       :: n
  real                           :: z(size(a,1),n)

  integer, allocatable          :: isuppz(:), iwork(:)
  real             , allocatable :: alocal(:,:)
  real             , allocatable :: w(:),work(:)
  integer                       :: lwork,liwork
  integer                       :: il,iu,info,m
  real                           :: abstol,vl,vu
  real                           :: ztmp(size(a,1),n)
!
!  workspace calculation
!
  vl=0.0
  vu=0.0
  il=1
  iu=size(a,1)
  abstol=0.0
  lwork=-1
  liwork=-1
  allocate(work(1),iwork(1),isuppz(2*iu),w(iu),alocal(iu,iu))
  call ssyevr('V','I','U',iu,a,iu,vl,vu,il,iu,abstol, &
     m,w,z,iu,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    print*,'DSYEVR returns with info= ',info
    call exit(1)
  endif
  lwork=int(work(1))
  liwork=iwork(1)
  deallocate(work,iwork)
  allocate(work(lwork),iwork(liwork))
!
!  set il and iu to correct values
!
  il=iu-n+1
!
!  copy a into alocal
!
  alocal = a
  call ssyevr('V','I','L',iu,alocal,iu,vl,vu,il,iu,abstol, &
     m,w,ztmp,iu,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    print*,'DSYEVR returns with info= ',info
    call exit(1)
  endif
  deallocate(work,iwork,alocal)
!
! reorder eigenvectors into descending order
!
  iu=size(a,1)
  do il=1,n
    z(:,il)=ztmp(:,n-il+1)
  end do
  end function evecs
end module eigenutils
