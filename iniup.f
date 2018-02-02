C  ****************************** INIUP **********************
C     THIS ROUTINE IS TO GIVE INITIAL FLOW FIELDS 
C     ON A PARABOLIC PROFILES
C     THIS SUBROUTINE IS MODIFIED BY JICHOI 99/07/27
C     TO MAKE CONSTANT FLOW RATE AT EACH PLANE
C     U IN YZ-PLANE, W IN XY-PLANE
C     THE FILE IS PART OF COMPACT-3D

      SUBROUTINE INIUP(U,P,HH,SS,PRESG,PRESG3)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/PARA/RE
      COMMON/VPERIN/VPER
      COMMON/SIZE/ALX,ALY,ALZ,VOL

      REAL U(0:M1,0:M2,0:M3,3)
      REAL HH(0:M1,0:M2,0:M3,3)
      REAL SS(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL RDNUM(3)
!      INTEGER ISEED(4)
      REAL NXZ

!      ISEED(1)=1001
!      ISEED(2)=2001
!      ISEED(3)=3001
!      ISEED(4)=4001

C     Random Number Generation in MKL
C     Math Kernel Library Manual pp.5-182 ~5-183
C     call dlarnv(idist, iseed, n, x)
C     idist INTEGER.iseed On exit, the seed is updated.
C          =1:uniform(0,1)
C          =2:uniform(-1,1)
C          =3:normal(0,1)
C     iseed INTEGER, Array, DIMENSION(4)
C           On entry, the seed of the random number generator;     
C           the array elements
C           must be between 0 and 4095, and iseed(4) must be odd.
C     n integer. The number of random numbers to be generated.
C     x DOUBLE PRECISION for dlarnv, Array, DIMENSION(n)
C           The generated random numbers.
C     initialization U
!      U=0.
!      PI=ACOS(-1.) 
!      DO 10 K=1,N3M
!      DO 10 J=1,N2M
!      DO 10 I=1,N1M
!      CALL DLARNV(2,ISEED,3,RDNUM)
!      U(I,J,K,1)=VPER*RDNUM(1)
!      U(I,J,K,2)=VPER*RDNUM(2)
!      U(I,J,K,3)=VPER*RDNUM(3)
!10    CONTINUE

      call random_seed ()           
      U=0.
      PI=ACOS(-1.) 
!$omp parallel do private(ru11,ru21,ru31,ru1,ru2,ru3)
      DO 10 I=1,N1M
      DO 10 J=1,N2M
      DO 10 K=1,N3M
      call random_number (ru11)
      ru1=-1.0+2.0*ru11
      U(I,J,K,1)=ru1 !*VPER
      call random_number (ru21)
      ru2=-1.0+2.0*ru21
      U(I,J,K,2)=ru2 !*VPER
      call random_number (ru31)
      ru3=-1.0+2.0*ru31
      U(I,J,K,3)=ru3 !*VPER 
10    CONTINUE



C     ELIMINATE MEAN QUANTITIES OF RANDOM FLUCTUATIONS
!$omp parallel do private(V1M)
      DO 21 I=1,N1M
      V1M=0.
      DO 22 K=1,N3M
      DO 22 J=1,N2M
      V1M=V1M+U(I,J,K,1)*DY(J)/DX3
22    CONTINUE
      V1M=V1M/ALY/ALZ
      DO 23 K=1,N3M
      DO 23 J=1,N2M
      U(I,J,K,1)=U(I,J,K,1)-V1M
23    CONTINUE
21    CONTINUE

!$omp parallel do private(V2M)
      DO 31 J=2,N2M
      V2M=0.
      DO 32 K=1,N3M
      DO 32 I=1,N1M
      V2M=V2M+U(I,J,K,2)/DX1/DX3
32    CONTINUE
      V2M=V2M/ALX/ALZ
      DO 33 K=1,N3M
      DO 33 I=1,N1M
      U(I,J,K,2)=U(I,J,K,2)-V2M
33    CONTINUE
31    CONTINUE

!$omp parallel do private(V3M)
      DO 41 K=1,N3M
      V3M=0.
      DO 42 J=1,N2M
      DO 42 I=1,N1M
      V3M=V3M+U(I,J,K,3)*DY(J)/DX1
42    CONTINUE
      V3M=V3M/ALX/ALY
      DO 43 J=1,N2M
      DO 43 I=1,N1M
      U(I,J,K,3)=U(I,J,K,3)-V3M
43    CONTINUE
41    CONTINUE

C     IMPOSE LAMINAR VELOCITY PROFILES IN U VELOCITIESA

      RET=0.09*RE**0.88  ! FROM REC TO RE_TAU
      UTAU_INI=RET/RE

!$omp parallel do private (YH,JP,REALU)
      DO 30 K=1,N3M
      DO 30 J=1,N2M
      JP=J+1
      DO 30 I=1,N1M
      YH=0.5*(Y(J)+Y(JP))
      REALU=YH*(ALY-YH)*(ALY-1.0)+YH*(2*ALY-YH)*(2.0-ALY)  ! laminar flow, full channel 0828 (NON-DIMENSIONAL)
! !     REALU=   !laminar flow, half channel 0828
!      IF(J.LE.N2M/2)THEN
!      REALU=1.0-UTAU_INI*(-1.0/0.41*ALOG(YH)+0.2)  !log law
!      ELSE
!      REALU=1.0-UTAU_INI*(-1.0/0.41*ALOG((ALY-YH))+0.2) !log law
!      ENDIF      

      !U(I,J,K,1)=U(I,J,K,1)+REALU

      U(I,J,K,1)=U(I,J,K,1)*VPER*REALU+REALU
      U(I,J,K,2)=U(I,J,K,2)*VPER   
      U(I,J,K,3)=U(I,J,K,3)*VPER
30    CONTINUE


C     CHECK FLOW RATE
      FLOW1=0.0
c      FLOW3=0.0
!$omp parallel do reduction(+:FLOW1)
      DO 40 K=1,N3M
      DO 40 J=1,N2M
      DO 40 I=1,N1M
      FLOW1=FLOW1+U(I,J,K,1)*DY(J)/DX1/DX3
c      FLOW3=FLOW3+U(I,J,K,3)*DY(J)/DX1/DX3
40    CONTINUE
      FLOW1=FLOW1/ALX/ALZ
c      FLOW3=FLOW3/ALX/ALZ/ALY

      FLOWR=(2*ALY)**3.0/6.0/2.0*(2.0-ALY)+(ALY)**3.0/6.0*(ALY-1.0) ! LAMINAR MASS FLOW RATE IN X1 DIRECTION
!$omp parallel do
      DO 45 K=1,N3M
      DO 45 J=1,N2M
      DO 45 I=1,N1M
      U(I,J,K,1)=FLOWR/FLOW1*U(I,J,K,1)
      U(I,J,K,3)=FLOWR/FLOW1*U(I,J,K,3)
   45 CONTINUE

C     IMPOSE ZERO VELOCITY AT BOUNDARY
!$omp parallel do
      DO 15 K=1,N3M  
      DO 15 I=1,N1M
      U(I,0,K,1)=0.0
      U(I,N2,K,1)=U(I,N2M,K,1)*(2.0-ALY)
      U(I,1,K,2)=0.0
      U(I,N2,K,2)=0.0
      U(I,0,K,3)=0.0
      U(I,N2,K,3)=U(I,N2M,K,3)*(2.0-ALY)
15    CONTINUE

C     IMPOSE ZERO-PRESSURE FLUCTUATIONS
!$omp parallel do
      DO 60 K=1,N3M
      DO 60 J=1,N2M
      DO 60 I=1,N1M
      P(I,J,K)=0.0
   60 CONTINUE
      
C     IMPOSE ZERO HH

      DO 65 NV=1,3
!$omp parallel do     
      DO 65 K=1,N3M
      DO 65 J=1,N2M
      DO 65 I=1,N1M
      HH(I,J,K,NV)=0.0
      SS(I,J,K,NV)=0.0
   65 CONTINUE
      
C     INITIAL MEAN PRESSURE GRADIENT AT LAMINAR FLOW FIELD
      PRESG=0.0  
      PRESG3=0.0      ! Z-direction pressure gradient 

      OPEN(201,FILE='INIU.PLT')
      WRITE(201,*)'VARIABLES="X","Y","Z","U","V","W"'
      WRITE(201,*)'ZONE T="',1,'" I=',N1M ,' J=',N2+1,' K=',N3M
      DO K=1,N3M
      DO J=0,N2
      DO I=1,N1M
      IF(J.NE.N2)THEN
      WRITE(201,333) (I-1)/DX1,0.5*(Y(J)+Y(J+1)),(K-1)/DX3,
     >            U(I,J,K,1),U(I,J,K,2),U(I,J,K,3)
      ENDIF
      IF(J.EQ.N2)THEN
      WRITE(201,333) (I-1)/DX1,Y(J),(K-1)/DX3,
     >            U(I,J,K,1),U(I,J,K,2),U(I,J,K,3)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
333   FORMAT(6(E12.5,2X))
      CLOSE(201)
 
      
      
      RETURN
      END
