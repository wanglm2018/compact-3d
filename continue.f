C  ****************************** continue **********************
C     THIS ROUTINE IS TO WRITE AND READ THE FLOW  
C     THIS SUBROUTINE IS MODIFIED BY JICHOI 99/07/27
C     TO MAKE CONSTANT FLOW RATE AT EACH PLANE
C     U IN YZ-PLANE, W IN XY-PLANE
C     THE FILE IS PART OF COMPACT-3D
C  ************************CWRITE*************************
C     WRITE FLOW FIELD AND BOUNDARY CONDITIONS
C     AND MEAN PRESSURE GRADIENT

      SUBROUTINE CWRITE(U,P,HH,SS,PRESG,PRESG3,NTIME,TIME)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/FILENAME/fileini,filegrd,fileout,fileavg
      CHARACTER*20 fileini,filegrd,fileout,fileavg

      REAL U(0:M1,0:M2,0:M3,3)
      REAL HH(0:M1,0:M2,0:M3,3)
      REAL SS(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)

      OPEN(3,FILE='CONTINUE')
      WRITE(3,*) PRESG,PRESG3
      WRITE(3,*) NTIME
      WRITE(3,*) TIME
      DO K=1,N3M
      DO J=0,N2
      DO I=1,N1M
      WRITE(3,*) U(I,J,K,1),U(I,J,K,2),U(I,J,K,3)
      ENDDO
      ENDDO
      ENDDO
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
      WRITE(3,*) P(I,J,K)
      ENDDO
      ENDDO
      ENDDO
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
      WRITE(3,*) HH(I,J,K,1),HH(I,J,K,2),HH(I,J,K,3)
      WRITE(3,*) SS(I,J,K,1),SS(I,J,K,2),SS(I,J,K,3)
      ENDDO
      ENDDO
      ENDDO
      CLOSE(3)

      RETURN
      END   
C  ************************  READUP **********************
C     READ FLOW FIELD AND BOUNDARY CONDITIONS
C     AND MEAN PRESSURE GRADIENT

      SUBROUTINE READUP(U,P,HH,SS,PRESG,PRESG3)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/FILENAME/fileini,filegrd,fileout,fileavg
      CHARACTER*20 fileini,filegrd,fileout,fileavg

      REAL U(0:M1,0:M2,0:M3,3)
      REAL HH(0:M1,0:M2,0:M3,3)
      REAL SS(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)

      OPEN(3,FILE=fileini,FORM='UNFORMATTED',STATUS='OLD')
      READ(3) PRESG1,PRESG3
      READ(3) (((U(I,J,K,1),U(I,J,K,2),U(I,J,K,3)
     >            ,K=1,N3M),J=0,N2),I=1,N1M)
      READ(3) (((P(I,J,K),K=1,N3M),J=1,N2M),I=1,N1M)
      READ(3) (((HH(I,J,K,1),HH(I,J,K,2),HH(I,J,K,3)
     >            ,K=1,N3M),J=1,N2M),I=1,N1M)
      READ(3) (((SS(I,J,K,1),SS(I,J,K,2),SS(I,J,K,3)
     >            ,K=1,N3M),J=1,N2M),I=1,N1M)  
      CLOSE(3)

      RETURN
      END