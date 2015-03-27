module geom
contains
  function rmsd(a,b,n) result(rms)
!
! finds the rmsd between the coordinates in a and those in b
!
  implicit none
  integer :: n
  real    :: a(3*n),b(3*n)
  real    :: r(3,3),v(3),rms
  logical :: findmove

  findmove=.false.
  r=0.0
  v=0.0
  call matfit(n,a,b,r,v,rms,findmove)
  end function rmsd

  function lsqfit(a,b,n,w) result(c)
!
! fits the coordinate set b to set a, returning the result in c
! w is an optional weighting array
!
  implicit none
  real    :: a(3*n),b(3*n),c(3*n)
  real,optional :: w(n)
  real    :: r(3,3),v(3),rmsd
  logical :: findmove
  integer :: i,j,k,n

  findmove=.true.
  r=0.0
  v=0.0

  if(present(w)) then
    call matfitw(n,a,b,r,v,rmsd,findmove,w)
  else
    call matfit(n,a,b,r,v,rmsd,findmove)
  endif
  c=transform(b,n,r,v)
  end function lsqfit
        
  FUNCTION TAU(R1,R2,R3,R4)
!
! calculates the torsion angle defined by r1-r2-r3-r4
!
  implicit none
  real    ::  R1(3),R2(3),R3(3),R4(3),U1(3),U2(3),U3(3)
  real    ::  V1(3),V2(3),V3(3),ST(3),TT(3)
  integer :: j
  real    :: s,t,tau
  DO J=1,3
    U1(J)=R2(J)-R1(J)
    U2(J)=R3(J)-R2(J)
    U3(J)=R4(J)-R3(J)
  end do
  CALL XPROD(U1,U2,V1)
  CALL XPROD(U2,U3,V2)
  CALL XPROD(V1,V2,V3)
  DO J=1,3
    ST(J)=V1(J)*V2(J)
    TT(J)=U2(J)*V3(J)
  end do
  S=ST(1)+ST(2)+ST(3)
  T=TT(1)+TT(2)+TT(3)
  S=MIN(MAX(S,-1.),1.)
  TAU=SIGN(180.0/ACOS(-1.0),T)*ACOS(S)
  END function tau

  SUBROUTINE XPROD(U,V,W)
!
! returns cross producr of u and v in w
!
  implicit none
  real     ::  U(3),V(3),W(3)
  real     :: s
  integer  ::  i,j,k
  DO I=1,3
    J=MOD(I,3)+1
    K=MOD(J,3)+1
    W(I)=U(J)*V(K)-U(K)*V(J)
  end do
  S=W(1)*W(1)+W(2)*W(2)+W(3)*W(3)
  IF(S.EQ.0.)RETURN
  S=1./SQRT(S)
  DO I=1,3
    W(I)=W(I)*S
  end do
  END subroutine xprod

  SUBROUTINE VANGLE(P,Q,N,THETA)
!
! returns angle (in degrees) between vectors p and q
!
  implicit none
  real     :: P(N),Q(N)
  real     :: D1,D2,R,C,THETA
  integer  :: I,N

  D1=0.0
  D2=0.0
  R=0.0
  DO I=1,N
    D1=D1+P(I)*P(I)
    D2=D2+Q(I)*Q(I)
    R=R+P(I)*Q(I)
  end do
  D1=SQRT(D1)
  D2=SQRT(D2)
  C=R/(D1*D2)
  IF (C.GT.1.0) C=1.0
  IF (C.LT.-1.0) C=-1.0
  THETA=ACOS(C)*180.0/ACOS(-1.0)
  END subroutine vangle

  real function torang(x,i)
!
! returns the angle (in degrees) between the four atoms in x
! indexed by i
!
  implicit none

  real    :: x(:,:)
  integer :: i(4)
  real    :: r1(3),r2(3),r3(3),r4(3)
  integer j

  do j=1,3
    r1(j)=x(j,i(1))
    r2(j)=x(j,i(2))
    r3(j)=x(j,i(3))
    r4(j)=x(j,i(4))
  end do
  torang=tau(r1,r2,r3,r4)
  end function torang

  real function cosd(t)
!
! cos of an angle in degrees
!
  implicit none
  real    :: t

  cosd=cos(t*acos(-1.0)/180.0)
  end function cosd

  real function sind(t)
!
! sine of an angle in degrees
!
  implicit none
  real    :: t

  sind=sin(t*acos(-1.0)/180.0)
  end function sind

  real function atan2d(c,s)
  implicit none
  real   :: c,s

  atan2d=atan2(c,s)*180.0/acos(-1.0)
  end function atan2d
   

  function rmat(i,t) result(r)
!
!  creates rotation matrix for rotation angle t about axis i
!  x=1, y=2, z=3
!
  real     ::  r(3,3),st,ct,t
  integer  :: i

  st=sind(t)
  ct=cosd(t)

  if(i.eq.1) then
    r(1,1)=1.0
    r(1,2)=0.0
    r(1,3)=0.0
    r(2,1)=0.0
    r(2,2)=ct
    r(2,3)=st
    r(3,1)=0.0
    r(3,2)=-st
    r(3,3)=ct
  else if(i.eq.2) then
    r(1,1)=ct
    r(1,2)=0.0
    r(1,3)=-st
    r(2,1)=0.0
    r(2,2)=1.0
    r(2,3)=0.0
    r(3,1)=st
    r(3,2)=0.0
    r(3,3)=ct
  else
    r(1,1)=ct
    r(1,2)=st
    r(1,3)=0.0
    r(2,1)=-st
    r(2,2)=ct
    r(2,3)=0.0
    r(3,1)=0.0
    r(3,2)=0.0
    r(3,3)=1.0
  endif
  return
  end function rmat

  function transform(xyz1,ntot,r,v) result(xyz2)
!
!  applies the rotation matrix r and vector v to the coordinates in xyz1
!
  real    ::  xyz1(3*ntot),xyz2(3*ntot),r(3,3),v(3)
  integer :: i,j,k,ji,ki

  do i=1,ntot
    do j=1,3
      ji=(i-1)*3+j
      xyz2(ji)=v(j)
      do k=1,3
        ki=(i-1)*3+k
        xyz2(ji)=xyz2(ji)+r(j,k)*xyz1(ki)
      end do
    end do
  end do
  end function transform

end module geom 
