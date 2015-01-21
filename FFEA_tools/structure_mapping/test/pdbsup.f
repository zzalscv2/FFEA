      program pdbsup                                                    PDBS0001
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
c !MS$DEBUG                                                             PDBS0010
                                                                        PDBS0011
c --- number of atoms read -                                            PDBS0012
      parameter (nat=10000)                                             PDBS0013
                                                                        PDBS0014
      character answ,matyp(nat)*2,fatyp(nat)*2,atyp*3                   PDBS0015
      character probe*50,a30*30,a12*12,s*8,fres(nat)*3,mres(nat)*3      PDBS0016
      character line*80,target*50,result*50,name*50,keywrd*6,res*3      PDBS0017
      real dm(4),vm(4,4),cm(3),cf(3)                                    PDBS0018
      real tr(3),t(3,3),q(4,4)                                          PDBS0019
      real dxp(nat,3),dxm(nat,3),xf(nat,3),xm(nat,3),x(3)               PDBS0020
      real*8 sumf(3),summ(3)                                            PDBS0021
                                                                        PDBS0022
      deg = 57.29577951                                                 PDBS0023
      zero=10E-30                                                       PDBS0024
                                                                        PDBS0025
c --- write header message                                              PDBS0026
 0003      format(a)                                                    PDBS0027
      write(*,3)' -----------------------------------------------------'PDBS0028
      write(*,3)'                    PROGRAM PDBSUP '                   PDBS0001
      write(*,3)' Determines rotation matrix and translation vector for'PDBS0002
      write(*,3)'  best fit superimposition of two pdb files by solving'PDBS0003
      write(*,3)'          the quaternion eigenvalue problem.'          PDBS0004
      write(*,3)'          Code by B.Rupp and S.Parkin (1996)'          PDBS0005
      write(*,3)'     Lawrence Livermore National Laboratory, BBRP'     PDBS0006
      write(*,3)'  Method by S.K.Kearsley, Acta Cryst. A45, 208 (1989)' PDBS0007
      write(*,3)'    See http://www-structure.llnl.gov for details.'    PDBS0008
      write(*,'(a,i6)')'   Max number of atoms per structure set to',natPDBS0009
      write(*,3)' -----------------------------------------------------'PDBS0010
                                                                        PDBS0011
 8001 write(*,'(/a$)')' Target input file (w/o *.pdb) : '               PDBS0012
      read(*,'(a)') name                                                PDBS0013
      target = name(1:index(name,' ')-1)//'.pdb'                        PDBS0014
      open (2,file=target, status='old',err=8001)                       PDBS0015
                                                                        PDBS0016
 8002 write(*,'(a$)')' Moving input file (w/o *.pdb) : '                PDBS0017
      read(*,'(a)') name                                                PDBS0018
      probe = name(1:index(name,' ')-1)//'.pdb'                         PDBS0019
      open (3,file=probe, status='old',err=8002)                        PDBS0020
                                                                        PDBS0021
      result = probe(1:index(probe,'.')-1)//'_mvd.pdb'                  PDBS0022
      open (6,file=result, status='unknown',err=9991)                   PDBS0023
                                                                        PDBS0024
      write(*,'(/a)') ' Target file name is   : '//target               PDBS0025
      write(*,'(a)')  ' Moving file name is   : '//probe                PDBS0026
      write(*,'(a/)') ' Moved  file name is   : '//result               PDBS0027
                                                                        PDBS0028
c --- read fixed target file coordinates ---                            PDBS0029
      ifix=0                                                            PDBS0030
      do i=1,nat                                                        PDBS0031
         read(2,'(a)',end=0001) line                                    PDBS0032
c     pdb atom format (6x,i5,1x,2a1,a2,a1,a3,1x,a1,i4,4x,3f8.3,2f6.2)   PDBS0033
         read(line,'(a6)') keywrd                                       PDBS0034
         if ((keywrd.eq.'ATOM  ').or.(keywrd.eq.'HETATM')) then         PDBS0035
            read(line,'(14x,a3,a3)')  atyp,res                          PDBS0036
c ---       disregard water molecules in any case                       PDBS0037
            if (res.ne.'HOH') then                                      PDBS0038
               ifix=ifix+1                                              PDBS0039
               read(line,'(30x,3f8.3)') (xf(ifix,j), j=1,3)             PDBS0040
               fres(ifix)=res                                           PDBS0041
c ---          fatyp is the atom descriptor INCLUDING branch indicator  PDBS0042
               fatyp(ifix)=atyp                                         PDBS0043
            end if                                                      PDBS0044
         end if                                                         PDBS0045
      end do                                                            PDBS0046
 0001 write(*,'(a,i5)')' Atoms read from target file : ',ifix           PDBS0047
                                                                        PDBS0048
c --- read moving probe coordinates --                                  PDBS0049
      imov=0                                                            PDBS0050
      do i=1,nat                                                        PDBS0051
         read(3,'(a)',end=0002) line                                    PDBS0052
         read(line,'(a6)') keywrd                                       PDBS0053
         if ((keywrd.eq.'ATOM  ').or.(keywrd.eq.'HETATM')) then         PDBS0054
            read(line,'(14x,a3,a3)')  atyp,res                          PDBS0055
            if (res.ne.'HOH') then                                      PDBS0056
               imov=imov+1                                              PDBS0057
               read(line,'(30x,3f8.3)') (xm(imov,j), j=1,3)             PDBS0058
               mres(imov)=res                                           PDBS0059
               matyp(imov)=atyp                                         PDBS0060
            end if                                                      PDBS0061
         end if                                                         PDBS0062
      end do                                                            PDBS0063
 0002 write(*,'(a,i5/)')' Atoms read from moving file : ',imov          PDBS0064
                                                                        PDBS0065
c --- whatever the case, imov will hold the number of atoms used        PDBS0066
      isok=1                                                            PDBS0067
      if (ifix.gt.imov) then                                            PDBS0068
         write(*,'(a)') ' *Warning* : number of atoms different '       PDBS0069
         write(*,'(a,i5,a/)') ' Will use only',imov,' atoms of probe'   PDBS0070
         ifix=imov                                                      PDBS0071
         isok=0                                                         PDBS0072
      else if (imov.gt.ifix) then                                       PDBS0073
         write(*,'(a)') ' *Warning* : number of atoms different '       PDBS0074
         write(*,'(a,i5,a/)') ' Will use only',ifix,' atoms of target'  PDBS0075
         imov=ifix                                                      PDBS0076
         isok=0                                                         PDBS0077
      end if                                                            PDBS0078
      if (isok.ne.0) then                                               PDBS0079
         write(*,'(a)')                                                 PDBS0080
     &   ' Check for identical number of atoms was ok'                  PDBS0081
      else                                                              PDBS0082
         write(*,'(a$)')' Do you want to continue (Y/N) : '             PDBS0083
         call inkey(answ)                                               PDBS0084
         if (answ.ne.'Y') goto 9992                                     PDBS0085
      end if                                                            PDBS0086
                                                                        PDBS0087
c --- check for pairwise correspondence of residues ---                 PDBS0088
      isok=1                                                            PDBS0089
      do i=1,imov                                                       PDBS0090
         if (fres(i).ne.mres(i)) then                                   PDBS0091
            write (*,'(a,i5,a,a)')' Mismatch at atom number',i,         PDBS0092
     &      ', offending residue ',mres(i)                              PDBS0093
            isok=0                                                      PDBS0094
         end if                                                         PDBS0095
      end do                                                            PDBS0096
      if (isok.ne.0) then                                               PDBS0097
         write(*,'(a)')                                                 PDBS0098
     &   ' Check for pairwise correspondence of residues was ok'        PDBS0099
      else                                                              PDBS0100
         write(*,'(a$)')' Do you want to continue (Y/N) : '             PDBS0101
         call inkey(answ)                                               PDBS0102
         if (answ.ne.'Y') goto 9992                                     PDBS0103
      end if                                                            PDBS0104
                                                                        PDBS0105
c --- check for pairwise correspondence of atoms ---                    PDBS0106
      isok=1                                                            PDBS0107
      do i=1,imov                                                       PDBS0108
         if (fatyp(i).ne.matyp(i)) then                                 PDBS0109
            write (*,'(a,i5,a,a)')' Mismatch at atom number',i,         PDBS0110
     &      ', offending atom type ',matyp(i)                           PDBS0111
            isok=0                                                      PDBS0112
         end if                                                         PDBS0113
      end do                                                            PDBS0114
      if (isok.ne.0) then                                               PDBS0115
         write(*,'(a)')                                                 PDBS0116
     &   ' Check for pairwise correspondence of atom types was ok'      PDBS0117
      else                                                              PDBS0118
         write(*,'(a$)')' Do you want to continue (Y/N) : '             PDBS0119
         call inkey(answ)                                               PDBS0120
         if (answ.ne.'Y') goto 9992                                     PDBS0121
      end if                                                            PDBS0122
                                                                        PDBS0123
c --- initialize all --                                                 PDBS0124
      do i=1,3                                                          PDBS0125
         sumf(i)=0                                                      PDBS0126
         summ(i)=0                                                      PDBS0127
      end do                                                            PDBS0128
      call filmat(4,4,q,0)                                              PDBS0129
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
      write(*,'(/a,3f8.3)')' Centre of target molecule  =',(cf(i),i=1,3)PDBS0144
      write(*,'(a,3f8.3)') ' Centre of moving molecule  =',(cm(i),i=1,3)PDBS0145
      write(*,'(a,3f8.3/)')' T - vector probe -> target =',(tr(i),i=1,3)PDBS0146
                                                                        PDBS0147
      write (*,'(a)')' Creating coordinate differences.......'          PDBS0148
c --- create coordinate differences delta x plus (dxp) and minus (dxm)  PDBS0149
      do k=1,imov                                                       PDBS0150
         do j=1,3                                                       PDBS0151
            dxm(k,j)=xm(k,j)-cm(j)-(xf(k,j)-cf(j))                      PDBS0152
            dxp(k,j)=xm(k,j)-cm(j)+(xf(k,j)-cf(j))                      PDBS0153
         end do                                                         PDBS0154
      end do                                                            PDBS0155
                                                                        PDBS0156
c --- fill upper triangle of (symmetric) quaternion matrix --           PDBS0157
      write (*,'(a)')' Filling quaternion matrix ............'          PDBS0158
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
      write (*,'(/a)')                                                  PDBS0175
     &     '       q(1)         q(2)         q(3)        q(4)'          PDBS0176
      do i=1,4                                                          PDBS0177
         write(*,'(4e13.5)') (q(i,j),j=1,4)                             PDBS0178
      end do                                                            PDBS0179
                                                                        PDBS0180
c --- orthogonalization by jacobi rotation = solution of EV -problem -- PDBS0181
      write (*,'(/a)')' Jacobi orthogonalization ..........'            PDBS0182
      n=4                                                               PDBS0183
      ns=4                                                              PDBS0184
      call jacobi(q,n,ns,dm,vm,nmrot)                                   PDBS0185
c --- sort eigenvectors after eigenvalues, descending --                PDBS0186
      write (*,'(a/)')' Sorting eigenvalues/vectors .......'            PDBS0187
      call eigsrt(dm,vm,n,ns)                                           PDBS0188
      write (*,'(a,i2,a)')' Eigenvalues and Eigenvectors (',            PDBS0189
     & nmrot,' Jacobi rotations)'                                       PDBS0190
      write (*,'(a)') '      e(1)        e(2)        e(4)        e(4)'  PDBS0191
      write (*,'(4e12.5,i5)') (dm(j),j=1,4)                             PDBS0192
      write (*,'(a)') '      ev(1)       ev(2)       ev(3)       ev(4)' PDBS0193
      do i=1,4                                                          PDBS0194
         write(*,'(4f12.6)') (vm(i,j),j=1,4)                            PDBS0195
      end do                                                            PDBS0196
                                                                        PDBS0197
c --- the smallest eigenvector contains best fit srs                    PDBS0198
      rmsd=sqrt(abs(dm(4)/imov))                                        PDBS0199
      write(*,'(/a/)')                                                  PDBS0200
     & ' The smallest eigenvalue represents s.r.s. of best fit'         PDBS0201
      write(*,'(a)')                                                    PDBS0202
     & ' Constructing the best fit rotation matrix from associated'     PDBS0203
      write(*,'(a/)') ' eigenvector elements (last column).....'        PDBS0204
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
      do i=1,3                                                          PDBS0217
         write(*,'(3f11.5)') (t(i,j),j=1,3)                             PDBS0218
      end do                                                            PDBS0219
                                                                        PDBS0220
c --- reset dxm to store the individual rmsd's in it now -              PDBS0221
      call filmat(nat,3,dxm,0)                                          PDBS0222
                                                                        PDBS0223
c --- xm and xf are not translated                                      PDBS0224
      do k=1,imov                                                       PDBS0225
c ---    subtract cm                                                    PDBS0226
         xm(k,1:3)=xm(k,1:3)-cm                                         PDBS0227
c ---    rotate it                                                      PDBS0228
         call rotvec(3,xm(k,1:3),t)                                     PDBS0229
c ---    now add cf                                                     PDBS0230
         xm(k,1:3)=xm(k,1:3)+cf                                         PDBS0231
         do i=1,3                                                       PDBS0232
            dxm(k,i)=sqrt((xf(k,i)-xm(k,i))**2)                         PDBS0233
         end do                                                         PDBS0234
      end do                                                            PDBS0235
                                                                        PDBS0236
c --- write header to moved file                                        PDBS0237
      write(*,'(/a,a,a/)')' Writing new pdb file ',                     PDBS0238
     & result(1:index(result,' ')),'..........'                         PDBS0239
      s='REMARK  '                                                      PDBS0240
      write(6,3)'TITLE   '//result(1:index(result,' '))                 PDBS0241
      write(6,3)s//'This pdb file has been created by superposition of 'PDBS0242
      write(6,3)s//'Target : '//target(1:index(target,' '))//' and'     PDBS0243
      write(6,3)s//'Source : '//probe(1:index(probe,' '))//' applying'  PDBS0244
      write(6,3)s//'the quaternion method S.K.Kearsley Acta Cryst '     PDBS0245
     & //'A45 208 (1989)'                                               PDBS0246
      write(6,3)s//'implemented in code PDBSUP by B.Rupp and S.Parkin'  PDBS0247
      write(6,3)s//'Lawrence Livermore National Laboratory, BBRP, 1996' PDBS0248
      write(6,3)s//'See http://www-structure.llnl.gov for details'      PDBS0249
      write(6,3)s                                                       PDBS0250
      write(6,3)s//'Translation vector: source -> target (Angstroem):'  PDBS0251
      write(6,'(a,3f12.4)')s, (tr(i),i=1,3)                             PDBS0252
      write(6,3)s//'Rotation matrix   : source -> target :'             PDBS0253
      do i=1,3                                                          PDBS0254
         write(6,'(a,3f12.6)')s, (t(i,j),j=1,3)                         PDBS0255
      end do                                                            PDBS0256
      write(6,'(a,f12.4,a)')                                            PDBS0257
     &  s//'The overall best fit r.m.s. deviation is : ',rmsd,' A'      PDBS0258
      write(6,3)s//'Last 3 columns contain rmsx,rmsy,rmsz,rmsd'         PDBS0259
      write(6,3)s                                                       PDBS0260
      rewind(3)                                                         PDBS0261
      idone=0                                                           PDBS0262
      do k=1,nat                                                        PDBS0263
         read(3,'(a)',end=9002) line                                    PDBS0264
         read(line,'(a6)') keywrd                                       PDBS0265
         if ((keywrd.eq.'ATOM  ').or.(keywrd.eq.'HETATM')) then         PDBS0266
            read(line,'(17x,a3)') res                                   PDBS0267
            if (res.ne.'HOH') then                                      PDBS0268
               idone=idone+1                                            PDBS0269
               if (idone.gt.imov) goto 9002                             PDBS0270
c     pdb atom format (6x,i5,1x,2a1,a2,a1,a3,1x,a1,i4,4x,3f8.3,2f6.2)   PDBS0271
               read(line,'(a30,3f8.3,a12)') a30,(x(i),i=1,3),a12        PDBS0272
c ---          center atom in origin subtract cm                        PDBS0273
               do i=1,3                                                 PDBS0274
                  x(i)=x(i)-cm(i)                                       PDBS0275
               end do                                                   PDBS0276
c ---          the thing is centred in origin now, rotate it            PDBS0277
               call rotvec (3,x,t)                                      PDBS0278
c ---          now add cf                                               PDBS0279
                  do i=1,3                                              PDBS0280
                  x(i)=x(i)+cf(i)                                       PDBS0281
               end do                                                   PDBS0282
c ---          r.m.s.d. = distance between target and probe             PDBS0283
               rmsd=sqrt(dxm(idone,1)+dxm(idone,2)+dxm(idone,3))        PDBS0284
               write(6,'(a30,3f8.3,a12,4f6.3)')                         PDBS0285
     &         a30,(x(i),i=1,3),a12,(dxm(idone,j),j=1,3),rmsd           PDBS0286
            end if                                                      PDBS0287
         else                                                           PDBS0288
            write(6,'(a)') line                                         PDBS0289
         end if                                                         PDBS0290
      end do                                                            PDBS0291
 9002 continue                                                          PDBS0292
                                                                        PDBS0293
      write(*,'(a,i5)')' Atoms read from probe file '//                 PDBS0294
     &probe(1:index(probe,' '))//':',imov                               PDBS0295
      write(*,'(a,i5/)') ' Atoms written to results file '//            PDBS0296
     &result(1:index(result,' '))//':',idone                            PDBS0297
      close(3)                                                          PDBS0298
      close(2)                                                          PDBS0299
      close(6)                                                          PDBS0300
                                                                        PDBS0301
c --- do we want to apply the same transformation to another molecule?  PDBS0302
      write(*,'(a$)')                                                   PDBS0303
     & ' Apply the same transformation to another molecule (Y/N) : '    PDBS0304
      call inkey(answ)                                                  PDBS0305
      if (answ.eq.'Y') then                                             PDBS0306
 8003    write(*,'(a$)')' File to be moved (w/o *.pdb) : '              PDBS0307
         read(*,'(a)') name                                             PDBS0308
         probe = name(1:index(name,' ')-1)//'.pdb'                      PDBS0309
         open (3,file=probe, status='old',err=8003)                     PDBS0310
         result = probe(1:index(probe,'.')-1)//'_mvd.pdb'               PDBS0311
         open (6,file=result, status='unknown',err=9991)                PDBS0312
         write(*,'(a)')  ' Moving file name is   : '//probe             PDBS0313
         write(*,'(a/)') ' Moved  file name is   : '//result            PDBS0314
         do k=1,nat                                                     PDBS0315
            read(3,'(a)',end=9003) line                                 PDBS0316
            read(line,'(a6)') keywrd                                    PDBS0317
            if ((keywrd.eq.'ATOM  ').or.(keywrd.eq.'HETATM')) then      PDBS0318
               read(line,'(17x,a3)') res                                PDBS0319
               if (res.ne.'HOH') then                                   PDBS0320
                  idone=idone+1                                         PDBS0321
c     pdb atom format (6x,i5,1x,2a1,a2,a1,a3,1x,a1,i4,4x,3f8.3,2f6.2)   PDBS0322
                  read(line,'(a30,3f8.3,a12)') a30,(x(i),i=1,3),a12     PDBS0323
c ---             center atom in origin subtract cm                     PDBS0324
                  do i=1,3                                              PDBS0325
                     x(i)=x(i)-cm(i)                                    PDBS0326
                  end do                                                PDBS0327
c ---             the thing is centred in origin now, rotate it         PDBS0328
                  call rotvec (3,x,t)                                   PDBS0329
c ---             now add cf                                            PDBS0330
                  do i=1,3                                              PDBS0331
                    x(i)=x(i)+cf(i)                                     PDBS0332
                  end do                                                PDBS0333
                  write(6,'(a30,3f8.3,a12,4f6.3)')                      PDBS0334
     &            a30,(x(i),i=1,3),a12,(dxm(idone,j),j=1,3),rmsd        PDBS0335
               end if                                                   PDBS0336
            else                                                        PDBS0337
               write(6,'(a)') line                                      PDBS0338
            end if                                                      PDBS0339
         end do                                                         PDBS0340
 9003    continue                                                       PDBS0341
         write(*,'(a,i5/)') ' Atoms written to results file '//         PDBS0342
     &   result(1:index(result,' '))//':',idone                         PDBS0343
         close(3)                                                       PDBS0344
         close(6)                                                       PDBS0345
      end if                                                            PDBS0346
      stop 'PDBSUP terminated normally'                                 PDBS0347
                                                                        PDBS0348
 9992 close(3)                                                          PDBS0349
      close(2)                                                          PDBS0350
      stop 'PDBSUP terminated  - atom mismatch'                         PDBS0351
                                                                        PDBS0352
 9991 stop ' Cannot create output file - tilt'                          PDBS0353
      end                                                               PDBS0354
                                                                        PDBS0355
      subroutine inkey (answ)                                           INKE0001
c ----------------------------------------------------------------------INKE0002
c     reads a key as answer                                             INKE0003
c ----------------------------------------------------------------------INKE0004
      character answ                                                    INKE0005
      read (*,'(a1)') answ                                              INKE0006
      call upstrg (answ,1)                                              INKE0007
      if ((answ.eq.' ').or.(answ.eq.char(13))) answ='Y'                 INKE0008
      return                                                            INKE0009
      end                                                               INKE0010
                                                                        INKE0011
      subroutine upstrg(strg,istrln)                                    UPST0001
c ----------------------------------------------------------------------UPST0002
c converts string str$ of lenght istrlen to upcase                      UPST0003
c ----------------------------------------------------------------------UPST0004
      character strg(istrln)                                            UPST0005
      integer      iascii                                               UPST0006
C-----change to upper case                                              UPST0007
         do 3031 i=1,istrln                                             UPST0008
            if (ichar(strg(i)).ge.95.and.ichar(strg(i)).le.122) then    UPST0009
               iascii=ichar(strg(i))-32                                 UPST0010
               strg(i)=char(iascii)                                     UPST0011
            endif                                                       UPST0012
 3031    continue                                                       UPST0013
      return                                                            UPST0014
      end                                                               UPST0015
                                                                        UPST0016
      subroutine trpmat(n,t,tr)                                         TRPM0001
c --- transpose matrix -------------------------------------------------TRPM0002
      real t(n,n), tr(n,n)                                              TRPM0003
      do i=1,n                                                          TRPM0004
         do j=1,n                                                       TRPM0005
            tr(j,i)=t(i,j)                                              TRPM0006
         end do                                                         TRPM0007
      end do                                                            TRPM0008
      return                                                            TRPM0009
      end                                                               TRPM0010
                                                                        TRPM0011
      subroutine filmat(n,m,r,ifil)                                     FILM0001
      real r(n,m)                                                       FILM0002
c --- initialize matrix ------------------------------------------------FILM0003
      do i=1,n                                                          FILM0004
         do j=1,m                                                       FILM0005
            r(i,j)=ifil                                                 FILM0006
         end do                                                         FILM0007
      end do                                                            FILM0008
      return                                                            FILM0009
      end                                                               FILM0010
                                                                        FILM0011
      subroutine rotvec (n,v,t)                                         ROTV0001
c --- multiply vector with matrix --------------------------------------ROTV0002
      real t(n,n), v(n),s(n)                                            ROTV0003
                                                                        ROTV0004
      do i=1,n                                                          ROTV0005
         s(i)=v(i)                                                      ROTV0006
         v(i)=0.0                                                       ROTV0007
      end do                                                            ROTV0008
      do i=1,n                                                          ROTV0009
         do j=1,n                                                       ROTV0010
            v(i)=v(i)+s(j)*t(i,j)                                       ROTV0011
         end do                                                         ROTV0012
      end do                                                            ROTV0013
      return                                                            ROTV0014
      end                                                               ROTV0015
                                                                        ROTV0016
      SUBROUTINE eigsrt(d,v,n,np)                                       EIGS0001
c ----------------------------------------------------------------------EIGS0002
      INTEGER n,np                                                      EIGS0003
      REAL d(np),v(np,np)                                               EIGS0004
      INTEGER i,j,k                                                     EIGS0005
      REAL p                                                            EIGS0006
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
      PARAMETER (nmax=500)                                              JACO0008
                                                                        JACO0009
      INTEGER n,np,nrot                                                 JACO0010
      REAL a(np,np),d(np),v(np,np)                                      JACO0011
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
      write(*,*) 'too many iterations in jacobi'                             JACO0103
      return                                                            JACO0104
      END                                                               JACO0105
C  (C) Copr. 1986-92 Numerical Recipes Software A2.Q2$2500.             JACO0106
