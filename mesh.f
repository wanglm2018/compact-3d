c***************** MESH ***********************  
c the part of compact-3d
c This file is part of lesgo.
      SUBROUTINE MESH
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/MESH4/DYM(M2),DYC(M2),DYP(M2)
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/MESH02/ XG(0:M1),YG(0:M2),ZG(0:M3),
     >	           XC(0:M1),YC(0:M2),ZC(0:M3)
      COMMON/FILENAME/fileini,filegrd,fileout,fileavg
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/WALLCO/HPW(M2),HMW(M2),HCW(M2)
      COMMON/FINOUT/INCODE,IDTOPT,NWRITE,NREAD,IAVG,NPRN,INSF,NINS,FILES
      COMMON/ROUGHNESS/Y0
      
      CHARACTER*20 fileini,filegrd,fileout,fileavg

c      CREATE THE UNIFORM GRID IN X2 DIRECTION      
!$omp parallel do
      DO I=1,N1M
      XG(I)=DBLE(I-1)*ALX/DBLE(M1-1)
      ENDDO
      XG(0)=0.0
      XG(M1)=ALX
      
!$omp parallel do      
      DO K=1,N3M
      ZG(K)=DBLE(K-1)*ALZ/DBLE(M3-1)
      ENDDO
      ZG(0)=0.0
      ZG(M3)=ALZ 
      
CCCCC    
      
      CKS=2.5      
      IF (WMODEL.EQ.0.0) THEN
      Y(0)=0.0
      YC(0)=0.0
      YG(1)=0.0
      ELSE
      Y(0)=Y0
      YC(0)=Y0
      YG(1)=Y0
      ENDIF      

      DO J=1,N2     
      IF (WMODEL.EQ.0.0) THEN 
      IF (ALY.EQ.2.) THEN
      Y(J)=1-TANH(CKS*(1-ALY*(J-1)/(N2-1)))/TANH(CKS)
	ELSE
	Y(J)=1-TANH(CKS*(1-2.*ALY*(J-1)/(2*N2-1)))/TANH(CKS)
	ENDIF
      ELSE
      Y(J)=(ALY-Y0)*(J-1)/(N2-1)+Y0
      ENDIF
      YG(J)=Y(J)
      ENDDO      
      ALY = Y(N2)

!$omp parallel do       
      DO I=1,N1M
      XC(I)=0.5*(XG(I+1)+XG(I))
      ENDDO
      XC(0)=0.0
      XC(M1)=XG(M1)  
   
!$omp parallel do
      DO K=1,N3M
      ZC(K)=0.5*(ZG(K+1)+ZG(K))
      ENDDO
      ZC(0)=0.0
      ZC(M3)=ZG(M3)

!$omp parallel do      
      DO J=1,N2M
      YC(J)=0.5*(Y(J+1)+Y(J))
      ENDDO
      YC(N2)=Y(N2)
      
      VOL=ALX*ALY*ALZ 
      DX1=DBLE(N1M)/ALX
      DX3=DBLE(N3M)/ALZ
      DX1Q=DX1**2.0
      DX3Q=DX3**2.0
      
      DY(1)=Y(2)
!$omp parallel do
      DO 20 J=2,N2M
      DY(J)=Y(J+1)-Y(J)
      H(J)=0.5*(DY(J)+DY(J-1))
 20   CONTINUE
      H(1)=0.5*DY(1)
      H(N2)=0.5*DY(N2M)

!$omp parallel do      
      DO 30 J=2,N2M-1
      HP(J)=1.0/H(J+1)/DY(J)
      HC(J)=(H(J+1)+H(J))/H(J+1)/H(J)/DY(J)
      HM(J)=1.0/H(J)/DY(J)
 30   CONTINUE 
 
!     AKSELVOLL&MOIN(1995)'S SPECIAL TREATMENT AT THE WALL.
!     SEE THE PAGE 33-34 OF REPORT NO. TF-63 (1995)
      HP(1)=4.0*H(1)/H(2)/(H(1)+H(2))/DY(1)
      HC(1)=4.0/H(2)/DY(1)
      HM(1)=4.0/(H(1)+H(2))/DY(1)
      HP(N2M)=4.0/(H(N2M)+H(N2))/DY(N2M)
      HC(N2M)=4.0/H(N2M)/DY(N2M)
      HM(N2M)=4.0*H(N2)/H(N2M)/(H(N2M)+H(N2))/DY(N2M)    
      

      IF(WMODEL.NE.0.)THEN
      HP(1)=1.0/H(2)/DY(1)
      HC(1)=1.0/H(2)/DY(1)
      HM(1)=0.0
      HP(N2M)=0.0
      HC(N2M)=1.0/H(N2M)/DY(N2M)
      HM(N2M)=1.0/H(N2M)/DY(N2M)
      ENDIF    
     
      DO 40 J=2,N2M
      DYP(J)=1.0/H(J)/DY(J)
      DYC(J)=1.0/H(J)*(1.0/DY(J)+1.0/DY(J-1))
      DYM(J)=1.0/H(J)/DY(J-1)
 40   CONTINUE
 
  
      OPEN(203,FILE='GRID.PLT',STATUS='UNKNOWN')
      DO J=0,N2
      IF(J.GT.0.AND.J.LT.N2)THEN
      DETAY=Y(J+1)-Y(J)
      ELSE 
      DETAY=0.0
      ENDIF
      WRITE(203,*)j,Y(J),DETAY
      END DO
      CLOSE(203)     
      
      RETURN
      END
c***************** INDICES ***********************     
      SUBROUTINE INDICES
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/FINDX3/FUP(M2),FDOWN(M2)

c!$omp parallel do
      DO 10 IC=1,N1M
      IPA(IC)=IC+1   ! I+1
 10   IMA(IC)=IC-1   ! I-1
      IPA(N1M)=1
      IMA(1)=N1M

c!$omp parallel do
      DO 20 KC=1,N3M
      KPA(KC)=KC+1
 20   KMA(KC)=KC-1
      KPA(N3M)=1
      KMA(1)=N3M

c!$omp parallel do
      DO 30 JC=1,N2M
      JPA(JC)=JC+1
      JMU(JC)=JC-1
 30   JMV(JC)=JC-1
      JPA(N2M)=N2M
      JMU(1)=1
      JMV(2)=2

c!$omp parallel do      
      DO 40 JC=1,N2M
      FJPA(JC)=JPA(JC)-JC   ! At boundary, these values are "0".
      FJMU(JC)=JC-JMU(JC)   ! Otherwise, these are "1".
 40   FJMV(JC)=JC-JMV(JC)

c!$omp parallel do                            
      DO 50 JC=1,N2M/2
      FDOWN(JC)=1
 50   FUP  (JC)=0

c!$omp parallel do 
      DO 60 JC=N2M/2+1,N2M
      FDOWN(JC)=0
 60   FUP  (JC)=1
                            
      RETURN
      END
