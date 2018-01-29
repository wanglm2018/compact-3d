c     FOURTH COMPACT DIFFERENCE DISCRETIZATION OF  
c     THE INCOMPRESSIBLE NAVIER-STOKES EQUATIONS
c     FOR THREE DIMENSIONAL PERIODIC TURBULENT CHANNEL FLOW
c
c     ORIGINALLY WRITTEN BY FLOW CONTROL LAB AT KAIST, KOREA
c     MODIFIED BY XIDIAN UNIVERSITY & LANZHOU UNIVERSITY
c
c     09/19/17 iwmles module for rough surface by rfhu
!!    Copyright (C)  2018  LanZhou University
!!
!!    This file is part of compact-3d.                                                 
      PROGRAM MAIN
      INCLUDE 'COMMON1.FI'
      include 'omp_lib.h'  
! Main file for compact-3d solver
! Contains main time loop

      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/WMEANO/WMO(M2,3),WRMSO(M2,3)
      COMMON/PMEANO/PMO(M2),PRMSO(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)     
!      COMMON/PARA/RE
      COMMON/FLOWRATE/FLOW1,FLOW3
c      COMMON/SIZE/ALX,ALY,ALZ,VOL
      
      COMMON/SPEC/XSPCTR(M1,M2,3),ZSPCTR(M3,M2,3)
      COMMON/SPECO/XSPCTRO(M1,M2,3),ZSPCTRO(M3,M2,3)
      COMMON/SPECTAU/SPECTAUW(M1,M3,2),SPECTAUB(M1,M3,2)
      COMMON/SPECTAUO/SPECTAUWO(M1,M3,2),SPECTAUBO(M1,M3,2)
      COMMON/CORR/RXX(M1,M2,3),RZZ(M3,M2,3)
      COMMON/CORRO/RXXO(M1,M2,3),RZZO(M3,M2,3)
      COMMON/CORR2/RXY(M1,M2,3)
      COMMON/CORR2O/RXYO(M1,M2,3)
      COMMON/QUAD/Q(M2,4),QP(M2,4)     
      COMMON/QUADO/QO(M2,4),QPO(M2,4)   

c-----for turbulent kinetic budgets 
      COMMON/VMEAN2O/DUDYO(M2),DVDYO(M2),DWDYO(M2)
      COMMON/PVMEANO/PVO(M2)
      COMMON/VHIGHO/VSKEWO(M2,3),VFLATO(M2,3),U2VO(M2),VW2O(M2)
      COMMON/PSTRO/PDUDXO(M2),PDVDYO(M2),PDWDZO(M2)
      COMMON/DISSUO/DUDX2O(M2),DUDY2O(M2),DUDZ2O(M2)
      COMMON/DISSVO/DVDX2O(M2),DVDY2O(M2),DVDZ2O(M2)
      COMMON/DISSWO/DWDX2O(M2),DWDY2O(M2),DWDZ2O(M2)      
      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)
      
c-----for SGS stress and dissipation 
      COMMON/SGSMEANO/STRMO(M2,6),SGSNUTO(M2),SGSDSO(M2),SMAGCO(M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)
      COMMON/WALLTAO/TAOW
      COMMON/WALLTAOAVE/TAOWO 
      COMMON/PARTICLE/DENRATIO,DIMP,PITIME,PDT     
      
!      COMMON/ROUGHNESS/Y0
c      COMMON/SIZESCALE/HHH
      
      REAL U(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL SGSVIS(0:M1,0:M2,0:M3)
      REAL IFWALL
      REAL HH(0:M1,0:M2,0:M3,3)  ! RIGHT-HAND SIDE CONVECTIVE TERM      
      REAL SS(0:M1,0:M2,0:M3,3)   ! RIGHT-HAND SIDE SUBGRID VISCOS TERM

      REAL UBULK,TBULK,TPLUS
      
      REAL    SEED
	INTEGER ISEED
      CALL CPU_TIME(SEED)
	ISEED=-NINT(RAN2(SEED)*1.E08)
      
!      nthds=num_parthds()
      nthds = omp_get_max_threads()
      
      tt0=rtc()      
      
      HHH=0.15      
      Y0=2.7E-5/HHH
      
      OPEN(31,FILE='WALLSS.plt',STATUS='UNKNOWN')
      OPEN(32,FILE='FLOW.plt',STATUS='UNKNOWN')
      OPEN(33,FILE='CFL.plt',STATUS='UNKNOWN')
      OPEN(34,FILE='RSH.plt',STATUS='UNKNOWN')
      OPEN(315,FILE='Retau.plt',STATUS='UNKNOWN')
      
      CALL SETUP  ! read inputs, mesh, index, etc.
      
      IF(NREAD.EQ.0) Call INIUP(U,P,HH,SS,PRESG,PRESG3)
      IF(NREAD.EQ.1) Call READUP(U,P,PRESG,PRESG3,NTIME,TIME) 
      
      TIME=0.0
      TBULK=0.0
      
c      WRITE(*,*) DIVMAX      
      
      IMORE=0                ! Index of written file
      IPLUS=0                ! Index of written file
      IFWALL=0.0
      NAVG=1
      NAVP=1
      DTR=DT
      IF(NREAD.EQ.1) THEN
      NTIME=NTIME
      ELSE
      NTIME=1
      ENDIF
      
      Init_P=1
      
      CALL DIVCHECK(U,TIME,DIVMAX,NTIME)  ! check divergence
      
      CALL CHKMF(U,TIME)  ! calculate flow rate  
      
!      DO 10 NTIME=1,NTST
20    t0=rtc()
      t1=mclock()
      t2 = omp_get_wtime()
      
      CALL CFL(U,CFLM) ! calc maximum cfl number
      DT=DTR
      IF (CFLM*DT.GE.CFLMAX.AND.IDTOPT.EQ.1) DT=CFLMAX/CFLM
c      IF (CFLM*DTR.LE.CFLMAX.OR.IDTOPT.NE.1)  DT=DTR
c      DT=DTR      
      IF(NTIME.LE.100.AND.NREAD.EQ.0) DT=0.001    
      
      TIME=TIME+DT      
      
      IF(IAVG.EQ.1.AND.NTIME.GT.75000.AND.MOD(NTIME,1).EQ.0) THEN ! NAVG pure flow field time average
          
      NAVG=NAVG+1     
      CALL ENERGY(U,P,NAVG)   ! x-z plane average of velocity, vorticity, pressure, etc       
      CALL SGSTRESS(U,SGSVIS)     !  x-z plane average of sgs stress       
      CALL TAVER(U,NAVG)    ! time-average of velocity, etc
      
      ELSE IF(NP_ALL.GT.0.AND.NTIME.GT.140000.AND.MOD(NTIME,2).EQ.0)
     > THEN         ! NAVP including particle flow field time average
      
      NAVP=NAVP+1
      CALL ENERGY(U,P,NAVP)     
      CALL SGSTRESS(U,SGSVIS)
      CALL TAVER(U,NAVP)
      
      ENDIF    
      
!      IF(NTIME.GT.90000.AND.MOD(NTIME,200).EQ.0) THEN
!      CALL WRITEUP(U,P,PRESG,PRESG3,NTIME) 
!      ENDIF

      IF(WMODEL.EQ.0.) CALL WALLSS(U,P,TIME)    !  x-z plane average of wall shear stress
      IF(WMODEL.NE.0.) TAOWO=TAOW              
      
      CALL GETUP(U,P,HH,SS,NTIME,PRESG,PRESG3,SGSVIS)     
      
      CALL WALLSS_WRITE(U,P,PRESG,TIME)    
      
      CALL DIVCHECK(U,TIME,DIVMAX,NTIME) 
      
      CALL CHKMF(U,TIME)      
      
      CALL CFL(U,CFLM)        

      IF(WMODEL.EQ.0)THEN
      IF(IAVG.EQ.1.AND.NTIME.LE.75000) THEN
      UTAU=SQRT(WSM(1))
      ELSE
      UTAU=SQRT(WSMO(1))
      ENDIF
      ELSE
      UTAU=SQRT(TAOWO)
      ENDIF 

      UBULK=FLOW1/ALZ/ALY
      TBULK=TIME*UBULK/ALX
      TPLUS=TIME*UTAU
      
      
      IF(MOD(NTIME,100).EQ.0) THEN
      WRITE(*,100) NTIME,PRESG,DT,TBULK,TPLUS,DIVMAX,CFLM*DT,NAVG,
     >              UTAU,UTAU*RE
      WRITE(*,'(A,I8,F7.4)') 'PARTICLE NUMBER AND FRICTION VELOCITY :'
     > ,NP_ALL,UTAU*RE*1.5/100000.0/HHH
c      t0=rtc()-t0
c      t1=mclock()-t1
      t3 = omp_get_wtime()
      write(6,*) '-----------------------------------'
      write(6,*) '-----------------------------------'
      write(6,*) 'Elapsed time',(t3-t2)*100
c      write(6,*) 'CPU     time',t1*1.e-2
      END IF     
      
      
      WRITE(33,*) TIME,DT,CFLM*DT
      
      IF(NREAD.EQ.0.AND.NTIME.EQ.100000) THEN
      CALL WRITEUP(U,P,PRESG,PRESG3,NTIME) 
      CALL WRITEAVG(U,P,TIME,NAVG,IPLUS,NTIME)      
      CALL CWRITE(U,P,PRESG,PRESG3,NTIME,TIME)
      ELSEIF(NTIME.GT.140000.AND.MOD(NTIME,10000).EQ.0) THEN
      CALL WRITEUP(U,P,PRESG,PRESG3,NTIME)
      CALL WRITEAVG(U,P,TIME,NAVP,IPLUS,NTIME)
      ENDIF
      
      NTIME=NTIME+1
 
      IF (NTIME.GT.NTST) GOTO 10
      GOTO 20


 10   CONTINUE

      CALL WRITEUP(U,P,PRESG,PRESG3,NTIME) 
      CALL WRITEAVG(U,P,TIME,NAVG,IPLUS,NTIME)      
!      CALL CWRITE(U,P,PRESG,PRESG3,NTIME,TIME)
      
C      IF (IAVG.EQ.1) CALL WRITEAVG(U,P,TIME,NAVG,IPLUS,NTIME)
C      CALL PROFILE(U,P,NTIME)

 100  FORMAT('STEP=',I6,2X,'DP=',E10.3,2X,'DT=',E12.5,2X,'TBULK=',E12.5,
     >        2X,'TPLUS=',E12.5,2X,'DIVMAX=',E10.3,2X,'CFL=',E10.3,2X,
     >        'NAVG=',I6.2,2X,'UTAU=',E10.3,2X,'RETAU=',E10.3)

      STOP
      END
