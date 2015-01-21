!#######################################################################
!#                                                                     #
!#              pyDockCG branch: new protein-protein code              #
!#                                                                     #
!#######################################################################
      subroutine superimpose(imov,xm,xf,cm,cf,t)
c ----------------------------------------------------------------------PDBS0002
c                               PROGRAM PDBSUP                          PDBS0001
c             Determines rotation matrix and translation vector for     PDBS0002
c              best fit superimposition of two pdb files by solving     PDBS0003
c                      the quaternion eigenvalue problem.               PDBS0004
c                      Code by B.Rupp and S.Parkin (1996)               PDBS0005
c               Lawrence Livermore National Laboratory, BBRP            PDBS0006
c             Method by S.K.Kearsley, Acta Cryst. A45, 208 (1989)       PDBS0007
c              See http://www-structure.llnl.gov for details.           PDBS0008
c ----------------------------------------------------------------------PDBS0009
      implicit none
                                                                        PDBS0014
      include 'DIMENSIONS2'
      integer imov,n,ns,nat,nmrot
      parameter (nat=maxres)       
      real*8 dm(4),vm(4,4),cm(3),cf(3)     
      real*8 tr(3),t(3,3),q(4,4)  
      real*8 dxp(nat,3),dxm(nat,3),xf(nat,3),xm(nat,3)!,x(3)
      real*8 sumf(3),summ(3) 
      integer i,j,k
      real*8 zero, zero2
      
      zero=10E-30          
      zero2=0.d0
c --- initialize all --                                                 PDBS0124
      do i=1,3                                                          PDBS0125
         sumf(i)=0                                                      PDBS0126
         summ(i)=0                                                      PDBS0127
      end do                                                            PDBS0128
      call filmat(4,4,q,zero2)                                          PDBS0129
                                                                        PDBS0130
c --- sum up all coordinates (in dble precision) to find centre ---     PDBS0131
      do k=1,imov                                                       PDBS0132
         do i=1,3                                                       PDBS0133
            sumf(i)=sumf(i)+xf(k,i)                                     PDBS0134
            summ(i)=summ(i)+xm(k,i)                                     PDBS0135
         end do                                                         PDBS0136
      end do                                                            PDBS0137
      do i=1,3                                                          PDBS0138
         cm(i)=sngl(summ(i)/imov)                                       PDBS0139
         cf(i)=sngl(sumf(i)/imov)                                       PDBS0140
         tr(i)=cf(i)-cm(i)                                              PDBS0141
      end do                                                            PDBS0142
                                                                        PDBS0143
c --- create coordinate differences delta x plus (dxp) and minus (dxm)  PDBS0149
      do k=1,imov                                                       PDBS0150
         do j=1,3                                                       PDBS0151
            dxm(k,j)=xm(k,j)-cm(j)-(xf(k,j)-cf(j))                      PDBS0152
            dxp(k,j)=xm(k,j)-cm(j)+(xf(k,j)-cf(j))                      PDBS0153
         end do                                                         PDBS0154
      end do                                                            PDBS0155
                                                                        PDBS0156
c --- fill upper triangle of (symmetric) quaternion matrix --           PDBS0157
      do k=1,imov                                                       PDBS0159
c ---    diags are sums of squared cyclic coordinate differences        PDBS0160
         q(1,1)=q(1,1)+dxm(k,1)**2+dxm(k,2)**2+dxm(k,3)**2              PDBS0161
         q(2,2)=q(2,2)+dxp(k,2)**2+dxp(k,3)**2+dxm(k,1)**2              PDBS0162
         q(3,3)=q(3,3)+dxp(k,1)**2+dxp(k,3)**2+dxm(k,2)**2              PDBS0163
         q(4,4)=q(4,4)+dxp(k,1)**2+dxp(k,2)**2+dxm(k,3)**2              PDBS0164
c ---    cross differences                                              PDBS0165
         q(1,2)=q(1,2)+dxp(k,2)*dxm(k,3)-dxm(k,2)*dxp(k,3)              PDBS0166
         q(1,3)=q(1,3)+dxm(k,1)*dxp(k,3)-dxp(k,1)*dxm(k,3)              PDBS0167
         q(1,4)=q(1,4)+dxp(k,1)*dxm(k,2)-dxm(k,1)*dxp(k,2)              PDBS0168
         q(2,3)=q(2,3)+dxm(k,1)*dxm(k,2)-dxp(k,1)*dxp(k,2)              PDBS0169
         q(2,4)=q(2,4)+dxm(k,1)*dxm(k,3)-dxp(k,1)*dxp(k,3)              PDBS0170
         q(3,4)=q(3,4)+dxm(k,2)*dxm(k,3)-dxp(k,2)*dxp(k,3)              PDBS0171
      end do                                                            PDBS0172
c --- fill the rest by transposing it onto itself                       PDBS0173
      call trpmat(4,q,q)                                                PDBS0174
                                                                        PDBS0180
c --- orthogonalization by jacobi rotation = solution of EV -problem -- PDBS0181
      n=4                                                               PDBS0183
      ns=4                                                              PDBS0184
      call jacobi(q,n,ns,dm,vm,nmrot)                                   PDBS0185
c --- sort eigenvectors after eigenvalues, descending --                PDBS0186
      call eigsrt(dm,vm,n,ns)                                           PDBS0188
                                                                        PDBS0205
c --- fill the rotation matrix which is made of elements from 4th EV    PDBS0206
      t(1,1)=vm(1,4)**2+vm(2,4)**2-vm(3,4)**2-vm(4,4)**2                PDBS0207
      t(2,1)=2*(vm(2,4)*vm(3,4)+vm(1,4)*vm(4,4))                        PDBS0208
      t(3,1)=2*(vm(2,4)*vm(4,4)-vm(1,4)*vm(3,4))                        PDBS0209
      t(1,2)=2*(vm(2,4)*vm(3,4)-vm(1,4)*vm(4,4))                        PDBS0210
      t(2,2)=vm(1,4)**2+vm(3,4)**2-vm(2,4)**2-vm(4,4)**2                PDBS0211
      t(3,2)=2*(vm(3,4)*vm(4,4)+vm(1,4)*vm(2,4))                        PDBS0212
      t(1,3)=2*(vm(2,4)*vm(4,4)+vm(1,4)*vm(3,4))                        PDBS0213
      t(2,3)=2*(vm(3,4)*vm(4,4)-vm(1,4)*vm(2,4))                        PDBS0214
      t(3,3)=vm(1,4)**2+vm(4,4)**2-vm(2,4)**2-vm(3,4)**2                PDBS0215
                                                                        PDBS0216
      !Aquesta és la matriu que volem (t(i,j))
      return 
      end                                                               PDBS0354
                                                                        PDBS0355
      subroutine trpmat(n,t,tr)                                         TRPM0001
c --- transpose matrix -------------------------------------------------TRPM0002
      implicit none
      integer i,j,n
      real*8 t(n,n), tr(n,n)                                              TRPM0003
      do i=1,n                                                          TRPM0004
         do j=1,n                                                       TRPM0005
            tr(j,i)=t(i,j)                                              TRPM0006
         end do                                                         TRPM0007
      end do                                                            TRPM0008
      return                                                            TRPM0009
      end                                                               TRPM0010
                                                                        TRPM0011
      subroutine filmat(n,m,r,ifil)                                     FILM0001
      implicit none
      integer i,j,n,m
      real*8 ifil
      real*8 r(n,m)                                                       FILM0002
c --- initialize matrix ------------------------------------------------FILM0003
      do i=1,n                                                          FILM0004
         do j=1,m                                                       FILM0005
            r(i,j)=ifil                                                 FILM0006
         end do                                                         FILM0007
      end do                                                            FILM0008
      return                                                            FILM0009
      end                                                               FILM0010
                                                                        FILM0011
      SUBROUTINE eigsrt(d,v,n,np)                                       EIGS0001
c ----------------------------------------------------------------------EIGS0002
      implicit none
      INTEGER n,np                                                      EIGS0003
      REAL*8 d(np),v(np,np)                                               EIGS0004
      INTEGER i,j,k                                                     EIGS0005
      REAL*8 p                                                            EIGS0006
      do 13 i=1,n-1                                                     EIGS0007
        k=i                                                             EIGS0008
        p=d(i)                                                          EIGS0009
        do 11 j=i+1,n                                                   EIGS0010
          if(d(j).ge.p)then                                             EIGS0011
            k=j                                                         EIGS0012
            p=d(j)                                                      EIGS0013
          endif                                                         EIGS0014
11      continue                                                        EIGS0015
        if(k.ne.i)then                                                  EIGS0016
          d(k)=d(i)                                                     EIGS0017
          d(i)=p                                                        EIGS0018
          do 12 j=1,n                                                   EIGS0019
            p=v(j,i)                                                    EIGS0020
            v(j,i)=v(j,k)                                               EIGS0021
            v(j,k)=p                                                    EIGS0022
12        continue                                                      EIGS0023
        endif                                                           EIGS0024
13    continue                                                          EIGS0025
      return                                                            EIGS0026
      END                                                               EIGS0027
C  (C) Copr. 1986-92 Numerical Recipes Software A2.Q2$2500.             EIGS0028
                                                                        EIGS0029
      SUBROUTINE jacobi(a,n,np,d,v,nrot)                                JACO0001
c ----------------------------------------------------------------------JACO0002
c     modified from numerical recipes book                              JACO0003
c     one needs to set the threshold for sm from sm.eq.0 to sm.lt.10E-30JACO0004
c     (anything in this range would be ok) due to underflow errors on   JACO0005
c     some computers/compilers.                                         JACO0006
c ----------------------------------------------------------------------JACO0007
      implicit none
      integer nmax
      PARAMETER (nmax=500)                                              JACO0008
                                                                        JACO0009
      INTEGER n,np,nrot                                                 JACO0010
      REAL*8 a(np,np),d(np),v(np,np)                                      JACO0011
      INTEGER i,ip,iq,j,maxrot                                          JACO0012
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax),zero          JACO0013
                                                                        JACO0014
c --- zero set and iteration maximum                                    JACO0015
      zero=10E-30                                                       JACO0016
      maxrot=50                                                         JACO0017
                                                                        JACO0018
      do 12 ip=1,n                                                      JACO0019
        do 11 iq=1,n                                                    JACO0020
          v(ip,iq)=0.                                                   JACO0021
11      continue                                                        JACO0022
        v(ip,ip)=1.                                                     JACO0023
12    continue                                                          JACO0024
      do 13 ip=1,n                                                      JACO0025
        b(ip)=a(ip,ip)                                                  JACO0026
        d(ip)=b(ip)                                                     JACO0027
        z(ip)=0.                                                        JACO0028
13    continue                                                          JACO0029
      nrot=0                                                            JACO0030
      do 24 i=1,maxrot                                                  JACO0031
        sm=0.                                                           JACO0032
        do 15 ip=1,n-1                                                  JACO0033
          do 14 iq=ip+1,n                                               JACO0034
            sm=sm+abs(a(ip,iq))                                         JACO0035
14        continue                                                      JACO0036
15      continue                                                        JACO0037
c ---   modified convergence threshold ---                              JACO0038
        if(sm.lt.zero)return                                            JACO0039
        if(i.lt.4)then                                                  JACO0040
          tresh=0.2*sm/n**2                                             JACO0041
        else                                                            JACO0042
          tresh=0.                                                      JACO0043
        endif                                                           JACO0044
        do 22 ip=1,n-1                                                  JACO0045
          do 21 iq=ip+1,n                                               JACO0046
            g=100.*abs(a(ip,iq))                                        JACO0047
            if((i.gt.4).and.(abs(d(ip))+                                JACO0048
     &g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then            JACO0049
              a(ip,iq)=0.                                               JACO0050
            else if(abs(a(ip,iq)).gt.tresh)then                         JACO0051
              h=d(iq)-d(ip)                                             JACO0052
              if(abs(h)+g.eq.abs(h))then                                JACO0053
                t=a(ip,iq)/h                                            JACO0054
              else                                                      JACO0055
                theta=0.5*h/a(ip,iq)                                    JACO0056
                t=1./(abs(theta)+sqrt(1.+theta**2))                     JACO0057
                if(theta.lt.0.)t=-t                                     JACO0058
              endif                                                     JACO0059
              c=1./sqrt(1+t**2)                                         JACO0060
              s=t*c                                                     JACO0061
              tau=s/(1.+c)                                              JACO0062
              h=t*a(ip,iq)                                              JACO0063
              z(ip)=z(ip)-h                                             JACO0064
              z(iq)=z(iq)+h                                             JACO0065
              d(ip)=d(ip)-h                                             JACO0066
              d(iq)=d(iq)+h                                             JACO0067
              a(ip,iq)=0.                                               JACO0068
              do 16 j=1,ip-1                                            JACO0069
                g=a(j,ip)                                               JACO0070
                h=a(j,iq)                                               JACO0071
                a(j,ip)=g-s*(h+g*tau)                                   JACO0072
                a(j,iq)=h+s*(g-h*tau)                                   JACO0073
16            continue                                                  JACO0074
              do 17 j=ip+1,iq-1                                         JACO0075
                g=a(ip,j)                                               JACO0076
                h=a(j,iq)                                               JACO0077
                a(ip,j)=g-s*(h+g*tau)                                   JACO0078
                a(j,iq)=h+s*(g-h*tau)                                   JACO0079
17            continue                                                  JACO0080
              do 18 j=iq+1,n                                            JACO0081
                g=a(ip,j)                                               JACO0082
                h=a(iq,j)                                               JACO0083
                a(ip,j)=g-s*(h+g*tau)                                   JACO0084
                a(iq,j)=h+s*(g-h*tau)                                   JACO0085
18            continue                                                  JACO0086
              do 19 j=1,n                                               JACO0087
                g=v(j,ip)                                               JACO0088
                h=v(j,iq)                                               JACO0089
                v(j,ip)=g-s*(h+g*tau)                                   JACO0090
                v(j,iq)=h+s*(g-h*tau)                                   JACO0091
19            continue                                                  JACO0092
              nrot=nrot+1                                               JACO0093
            endif                                                       JACO0094
21        continue                                                      JACO0095
22      continue                                                        JACO0096
        do 23 ip=1,n                                                    JACO0097
          b(ip)=b(ip)+z(ip)                                             JACO0098
          d(ip)=b(ip)                                                   JACO0099
          z(ip)=0.                                                      JACO0100
23      continue                                                        JACO0101
24    continue                                                          JACO0102
      print *,  'too many iterations in jacobi'                         JACO0103
      return                                                            JACO0104
      END                                                               JACO0105
C  (C) Copr. 1986-92 Numerical Recipes Software A2.Q2$2500.             JACO0106
