c***************** SETUP ***********************  
c the part of compact-3d
c read the file of parame.ter
      SUBROUTINE SETUP
      INCLUDE 'PARAM.H'
      
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/VPERIN/VPER
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/FINOUT/INCODE,IDTOPT,NWRITE,NREAD,IAVG,NPRN,INSF,NINS,FILES
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/NEARWAL/NEARWALL
      !!!!!!!!!!!!!!!!!!
      
      Common/PAR1/ IFPAR,Nstep_P,Nprint_P,NTIMEP 
      COMMON/PARTICLE/DENRATIO,DIMP,PITIME,PDT,NUMPI
      COMMON/FILENAME/fileini,filegrd,fileout,fileavg
      CHARACTER*20 fileini,filegrd,fileout,fileavg
      CHARACTER*70 DUMMY
      integer NEARWALL
      
      OPEN(1,FILE='parame.ter',STATUS='OLD')
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,302) DUMMY,FILES
      WRITE(*,302) DUMMY,FILES
      READ (1,302) DUMMY,WMODEL
      WRITE(*,302) DUMMY,WMODEL
      READ (1,301) DUMMY,NEARWALL
      WRITE(*,301) DUMMY,NEARWALL
      READ (1,301) DUMMY,N1
      WRITE(*,301) DUMMY,N1
      READ (1,301) DUMMY,N2
      WRITE(*,301) DUMMY,N2
      READ (1,301) DUMMY,N3
      WRITE(*,301) DUMMY,N3
      READ (1,302) DUMMY,RE
      WRITE(*,302) DUMMY,RE
      READ (1,302) DUMMY,ALX
      WRITE(*,302) DUMMY,ALX
      READ (1,302) DUMMY,ALZ
      WRITE(*,302) DUMMY,ALZ
      READ (1,301) DUMMY,INCODE
      WRITE(*,301) DUMMY,INCODE
      READ (1,301) DUMMY,NTST
      WRITE(*,301) DUMMY,NTST
      READ (1,302) DUMMY,VPER
      WRITE(*,302) DUMMY,VPER
      READ (1,302) DUMMY,DT
      WRITE(*,302) DUMMY,DT
      READ (1,301) DUMMY,IDTOPT
      WRITE(*,301) DUMMY,IDTOPT
      READ (1,302) DUMMY,CFLMAX
      WRITE(*,302) DUMMY,CFLMAX
      READ (1,301) DUMMY,NWRITE
      WRITE(*,301) DUMMY,NWRITE
      READ (1,301) DUMMY,NREAD
      WRITE(*,301) DUMMY,NREAD
      READ (1,301) DUMMY,IAVG
      WRITE(*,301) DUMMY,IAVG
      READ (1,301) DUMMY,NPRN
      WRITE(*,301) DUMMY,NPRN
      READ (1,301) DUMMY,INSF
      WRITE(*,301) DUMMY,INSF
      READ (1,301) DUMMY,NINS
      WRITE(*,301) DUMMY,NINS
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,303) DUMMY,fileini
      WRITE(*,303) DUMMY,fileini
      READ (1,303) DUMMY,filegrd
      WRITE(*,303) DUMMY,filegrd
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      

  300 FORMAT(A65)    
  301 FORMAT(A45,I15)
  302 FORMAT(A45,E15.7)
  303 FORMAT(A45,A20)

C------------------------------------
C     PHYSICAL LENGTH
      IF(INCODE.EQ.1) THEN
      PI=ACOS(-1.0)
      ALX=5.0*PI     
      ALZ=2.0*PI
      ALY=2.0     
      ENDIF
C------------------------------------
      
      N1M=N1-1
      N2M=N2-1
      N3M=N3-1

      CALL MESH
      CALL INDICES 
      CALL INIWAVE
     
      RETURN
      END