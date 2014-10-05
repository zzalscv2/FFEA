! Pczcomp: compare two '.pcz' format compressed MD trajectory files
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

program pczcomp

!
!  comparison of two .pcz files: rmsd between average structures
!  plus dotproduct matrix and subspace overlap
!
use pcz_io
use mygetargs
use trajutils
implicit none
  integer, parameter    :: STDERR=0,XUNIT=1,YUNIT=2
  real,allocatable      :: xavg(:),yavg(:),xydif(:),proj(:)
  real, allocatable     :: evalx(:),evaly(:),evx(:,:),evy(:,:)
  real, allocatable     :: dp(:,:),mxx(:),myy(:),mxy(:),myx(:)
  real                  :: r(3,3),v(3),rmsd,sm,smx,smy,smx2,smy2,xt,yt,zt,ss,pj
  real                  :: rmsxx,rmsxy,rmsyy,rmsyx,tx,ty
  integer               :: natoms,nvecs,nvcomp,stat,i,j,k,l
  character*160         :: xfilename,yfilename
  character*80          :: title
  logical               :: ok
  type(pczfile)         :: pczfx,pczfy

  ok=.true.
  if(.not.mygetarg('-x',xfilename)) then
    write(STDERR,*) 'Specify first pcz file with -x'
    ok=.false.
  endif
  if(.not.mygetarg('-y',yfilename)) then
    write(STDERR,*) 'Specify second pcz file with -y'
    ok=.false.
  endif
  if(.not.ok) call exit(1)
!
!  unless told to the contrary, look at the top 10 eigenvectors
!
  if(.not.mygetarg('-nv',nvcomp)) nvcomp=10
!
! open the two .pcz files, and extract the data
!
  pczfx=pczopen(XUNIT,xfilename)
  if(pczfx%status/=0) call die(STDERR,'Error opening '//trim(xfilename))
  pczfy=pczopen(YUNIT,yfilename)
  if(pczfy%status/=0) call die(STDERR,'Error opening '//trim(yfilename))
!
! bit of sanity checking
!
  if(pczfx%natoms/=pczfy%natoms) call die(STDERR,'Error: atom number mismatch')
  natoms=pczfx%natoms
  nvecs=max(pczfx%nvecs,pczfy%nvecs)
  if (nvcomp.gt.pczfx%nvecs) then
     write(STDERR,"('Error: ',a,'only contains',i3,' eigenvectors')") &
       trim(xfilename),pczfx%nvecs
    ok=.false.
  endif
  if (nvcomp.gt.pczfy%nvecs) then
     write(STDERR,"('Error: ',a,'only contains',i3,' eigenvectors')") &
       trim(yfilename),pczfy%nvecs
     ok=.false.
  endif
  if(.not.ok) call exit(1)
!
! ok, ready to go. Begin by reading average structures
! and calculating rmsd, plus r and v
!
  write(*,"('Comparison of X: ',a,/,'and Y: ',a)") trim(xfilename),trim(yfilename)
  allocate(xavg(3*natoms),yavg(3*natoms),xydif(3*natoms),STAT=stat)
  if(stat/=0) call die(STDERR,'Memory allocation error')
  xavg=pczavg(pczfx)
  yavg=pczavg(pczfy)
!
!  now read in sufficient sets of eigenvectors and eigenvalues
!
  allocate(evx(3*natoms,nvcomp),evy(3*natoms,nvcomp),STAT=stat)
  if(stat/=0) call die(STDERR,'Memory allocation error')
  allocate(evalx(pczfx%nvecs),evaly(pczfy%nvecs),STAT=stat)
  if(stat/=0) call die(STDERR,'Memory allocation error')
  do i=1,nvcomp
    evx(:,i)=pczevec(pczfx,i)
    evy(:,i)=pczevec(pczfy,i)
  end do
  evalx=pczevals(pczfx)
  evaly=pczevals(pczfy)
!
!  calculate mahalanobis distances
!
!  X from <X>, and Y from <Y>:
!
   allocate(mxx(pczfx%nframes),myy(pczfy%nframes))
   mxx=pczmaha(pczfx,nvcomp)
   myy=pczmaha(pczfy,nvcomp)
!
!  now y from <x>:
!
  allocate(myx(pczfy%nframes),mxy(pczfx%nframes))
  do k=1,pczfy%nframes
    xydif=pczsnap(pczfy,k)
    call trajfit(xydif,xavg,natoms,1,stat)
    xydif=xydif-xavg
    sm=0.0
    do j=1,nvcomp
      pj=0.0
      do i=1,natoms*3
        pj=pj+evx(i,j)*xydif(i)
      end do
      sm=sm+(pj*pj)/evalx(j)
    end do
    myx(k)=sqrt(sm)
  end do
!
!  now of x from <y>:
!
  do k=1,pczfx%nframes
    xydif=pczsnap(pczfx,k)
    call trajfit(xydif,yavg,natoms,1,stat)
    xydif=xydif-yavg
    sm=0.0
    do j=1,nvcomp
      pj=0.0
      do i=1,natoms*3
        pj=pj+evy(i,j)*xydif(i)
      end do
      sm=sm+(pj*pj)/evaly(j)
    end do
    mxy(k)=sqrt(sm)
  end do
!
! find averages and sds...
!
  smx=sum(mxx)/pczfx%nframes
  smy=sum(myy)/pczfy%nframes
  smx2=sum(mxx*mxx)/pczfx%nframes
  smy2=sum(myy*myy)/pczfy%nframes
  rmsxx=sqrt(smx2-smx*smx)
  rmsyy=sqrt(smy2-smy*smy)
!
! now do standard rmsd and dot product stuff
!
  ok=.true.
  call matfit(natoms,xavg,yavg,r,v,rmsd,ok)
  write(*,"('Rmsd between <x> and <Y>: ',f6.2)") rmsd
  write(*,"('Mahalanobis distances:')")
  write(*,"('X from <X>: ',f6.2,' +/-',f6.2)") smx,rmsxx
  write(*,"('Y from <Y>: ',f6.2,' +/-',f6.2)") smy,rmsyy
!
! convert to deviations from the mean...
!
  mxx=abs(mxx-smx)
  myy=abs(myy-smy)
! now cross-distances
  smx=sum(mxy)/pczfx%nframes
  smy=sum(myx)/pczfy%nframes
  smx2=sum(mxy*mxy)/pczfx%nframes
  smy2=sum(myx*myx)/pczfy%nframes
  rmsxx=sqrt(smx2-smx*smx)
  rmsyy=sqrt(smy2-smy*smy)
  write(*,"('X from <Y>: ',f6.2,' +/-',f6.2)") smx,rmsxx
  write(*,"('Y from <X>: ',f6.2,' +/-',f6.2)") smy,rmsyy
!
! convert to deviations from the mean...
!
  mxy=abs(mxy-smx)
  myx=abs(myx-smy)
!
! find 90% envelope
!
  call sort_ascending(mxx,pczfx%nframes)
  call sort_ascending(myy,pczfy%nframes)
  call sort_ascending(mxy,pczfx%nframes)
  call sort_ascending(myx,pczfy%nframes)
  tx=mxx(int(0.9*pczfx%nframes))
  ty=myy(int(0.9*pczfy%nframes))
  i=1
  do while((mxy(i)<tx).and.i<=pczfx%nframes)
    i=i+1
  end do
  write(*,"(i3,'% of X lies within the 90% envelope of Y')") &
         int(100*i/pczfx%nframes)
  
  i=1
  do while((myx(i)<ty).and.i<=pczfy%nframes)
    i=i+1
  end do
  write(*,"(i3,'% of Y lies within the 90% envelope of X')") &
         int(100*i/pczfy%nframes)
!
! now rotate y's eigenvectors
!
  do i=1,nvcomp
    do j=1,natoms*3,3
      k=j+1
      l=j+2
      xt=r(1,1)*evy(j,i)+r(1,2)*evy(k,i)+r(1,3)*evy(l,i)
      yt=r(2,1)*evy(j,i)+r(2,2)*evy(k,i)+r(2,3)*evy(l,i)
      zt=r(3,1)*evy(j,i)+r(3,2)*evy(k,i)+r(3,3)*evy(l,i)
      evy(j,i)=xt
      evy(k,i)=yt
      evy(l,i)=zt
    end do
  end do
!
! now calculate dot product matrix
!
  allocate(dp(nvcomp,nvcomp),STAT=stat)
  if(stat/=0) call die(STDERR,'Memory allocation error')
  ss=0.0
  do i=1,nvcomp
    do j=1,nvcomp
      dp(j,i)=abs(dot_product(evy(:,j),evx(:,i)))
      ss=ss+dp(j,i)*dp(j,i)
    end do
  end do
  ss=ss/nvcomp
!
!  now write everything out
!
  write(*,'(i5)') nvcomp
  do i=1,nvcomp
    write(*,'(20f5.2)') (dp(j,i),j=1,nvcomp)
  end do
  write(*,"('Subspace overlap: ',f6.3)") sqrt(ss)
  write(*,"('Average maximum dot product: ',f6.3)") &
    (sum(maxval(dp,DIM=1))+sum(maxval(dp,DIM=2)))/(2*nvcomp)
!
!  all done
!
  close(XUNIT)
  close(YUNIT)
  stat=pczclose(pczfx)
  stat=pczclose(pczfy)

end program pczcomp

subroutine die(unit,msg)
implicit none
  integer            :: unit
  character(len=*)   :: msg

  write(unit,*) msg
  call exit(1)
end subroutine die

subroutine sort_ascending(a,n)
implicit none
  real       :: a(n)
  integer    :: l(1),n,i,idx(n)
  logical    :: mask(n)

  mask=.true.

  do i=1,n
    l=minloc(a,mask)
    idx(i)=l(1)
    mask(l(1))=.false.
  end do
  a=a(idx)
  return
end subroutine sort_ascending
