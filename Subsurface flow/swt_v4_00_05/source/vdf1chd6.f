C
      SUBROUTINE VDF1CHD6AD(NCHDS,MXCHD,CHDS,HNEW,HOLD,PERLEN,PERTIM,
     1            NCOL,NROW,NLAY,NCHDVL,IOUT,IBOUND,BOTM,NBOTM)
C
C-----VERSION 11JAN2000 GWF1CHD6FM
C     ******************************************************************
C     COMPUTE HEAD FOR TIME STEP AT EACH TIME-VARIANT SPECIFIED HEAD
C     CELL
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE VDFMODULE,   ONLY: PS,ELEV
      DOUBLE PRECISION HNEW,DZERO,HB
      DIMENSION CHDS(NCHDVL,MXCHD),HNEW(NCOL,NROW,NLAY),
     1          HOLD(NCOL,NROW,NLAY)
C--SEAWAT: NEED IBOUND FOR CHDOPT 2
      DIMENSION IBOUND(NCOL,NROW,NLAY),BOTM(NCOL,NROW,0:NBOTM)
C--SEAWAT: INCLUDE THE CHDAUX TO CHECK FOR CHDDENSOPT
      COMMON /CHDCOM/CHDAUX(5)
      COMMON /DISCOM/LBOTM(999)
      CHARACTER*16 CHDAUX
C     ------------------------------------------------------------------
      DZERO=0.
C
C1------IF NCHDS<=0 THEN THERE ARE NO TIME VARIANT SPECIFIED-HEAD CELLS.
C1------RETURN.
      IF(NCHDS.LE.0) RETURN
C
C--SEAWAT: CHECK FOR PRESENCE OF CHDDENSOPT
      LOCCHDDENSOPT=0
      LOCCHDDEN=0
      DO I=1,5
        IF(CHDAUX(I).EQ.'CHDDENSOPT') LOCCHDDENSOPT=I+5
        IF(CHDAUX(I).EQ.'CHDDEN') LOCCHDDEN=I+5
      ENDDO
C6------INITIALIZE HNEW TO 0 AT SPECIFIED-HEAD CELLS.
      DO 50 L=1,NCHDS
      IL=CHDS(1,L)
      IR=CHDS(2,L)
      IC=CHDS(3,L)
      HNEW(IC,IR,IL)=DZERO
   50 CONTINUE
C
C2------COMPUTE PROPORTION OF STRESS PERIOD TO CENTER OF THIS TIME STEP
      IF (PERLEN.EQ.0.0) THEN
        FRAC=1.0
      ELSE
        FRAC=PERTIM/PERLEN
      ENDIF
C
C2------PROCESS EACH ENTRY IN THE SPECIFIED-HEAD CELL LIST (CHDS)
      DO 100 L=1,NCHDS
C
C3------GET COLUMN, ROW AND LAYER OF CELL CONTAINING BOUNDARY
      IL=CHDS(1,L)
      IR=CHDS(2,L)
      IC=CHDS(3,L)
C
      IF (PERLEN.EQ.0.0 .AND. CHDS(4,L).NE.CHDS(5,L)) THEN
        WRITE(IOUT,200)IL,IR,IC
 200    FORMAT(/,' ***WARNING***  FOR CHD CELL (',I3,',',I5,',',I5,
     &'), START HEAD AND END HEAD DIFFER',/,
     &' FOR A STRESS PERIOD OF ZERO LENGTH --',/,
     &' USING ENDING HEAD AS CONSTANT HEAD',
     &' (GWF1CHD6AD)',/)
      ENDIF
C5------COMPUTE HEAD AT CELL BY LINEAR INTERPOLATION.
      HB=CHDS(4,L)+(CHDS(5,L)-CHDS(4,L))*FRAC
C--SEAWAT: CHECK FOR CHD OPTIONS AND DETERMINE DENSITY VALUE
C--SEAWAT: SET DENSE = PS FOR DEFAULT OPTION
      DENSE=PS(IC,IR,IL)
      IF(LOCCHDDENSOPT.EQ.0) GOTO 80
      ICHDDENSOPT=INT(CHDS(LOCCHDDENSOPT,L))
      IF(ICHDDENSOPT.EQ.3) GOTO 90 !HB IS READ AS REFERENCE HEAD--DO NOTHING
      IF(ICHDDENSOPT.EQ.1) THEN
        IF(LOCCHDDEN.EQ.0) THEN
            WRITE(IOUT,*) 'STOPPING: CHD DENSITY OPTION SPECIFIED AS 1, 
     +BUT CHDDEN AUXILIARY VARIABLE NOT FOUND'
            WRITE(*,*) 'STOPPING: CHD DENSITY OPTION SPECIFIED AS 1, BUT
     + CHDDEN AUXILIARY VARIABLE NOT FOUND'
            STOP
        ENDIF
        DENSE=CHDS(LOCCHDDEN,L)
      ENDIF
C--SEAWAT: CALCULATE DENSITY VALUE FOR OPTION 2
      IF(ICHDDENSOPT.EQ.2) THEN  !THIS OPTION ASSUMES THAT HB IS THE ENVIRONMENTAL HEAD
        TOTDEPTH=HB-ELEV(IC,IR,IL)
        DENSE=PS(IC,IR,IL)
        IF(TOTDEPTH.LE.0.0.OR.IL.EQ.1) GOTO 80   !TOP CELL OR HB IS BELOW CELL CENTER
        DO K=1,IL   !FIND HIGHEST ACTIVE CELL
            IF(IBOUND(IC,IR,K).NE.0) THEN
                KTOP=K  
                EXIT
            ENDIF
        ENDDO
        IF(HB.GT.BOTM(IC,IR,LBOTM(1)-1)) THEN    !HB ABOVE TOP OF MODEL, ASSUME HB-TOP IS OF DENSE OF HIGHEST ACTIVE CELL
            WT=HB-BOTM(IC,IR,LBOTM(1)-1)  
            DENSE=WT*PS(IC,IR,KTOP)
        ELSE
            DENSE=0.
        ENDIF
        DO K=IL,1,-1   !CALCULATE AVERAGE DENSITY IN COLUMN OF WATER
            TOP=BOTM(IC,IR,LBOTM(K)-1)
            BOT=BOTM(IC,IR,LBOTM(K))
            IF(K.EQ.IL) BOT=ELEV(IC,IR,IL)
            WT=TOP-BOT
            IF(HB.LT.TOP) THEN
                IF (HB.GT.BOT) THEN
                    WT=HB-BOT                 !HB IN THIS CELL, WT IS SAT THICKNESS
                ELSE
                    WT=0.                     !HB BELOW THIS CELL, NO CONTRIBUTION TO AVERAGE DENSITY
                ENDIF
            ENDIF
            IF(IBOUND(IC,IR,K).EQ.0) THEN            !INACTIVE CELL DOES NOT CONTRIBUTE TO AVERAGE DENSITY OF COLUMN
                TOTDEPTH=TOTDEPTH-(TOP-BOT)
                CYCLE
            ENDIF
            DENSE=DENSE+WT*PS(IC,IR,K)
        ENDDO
   70   CONTINUE
        DENSE=DENSE/TOTDEPTH
      ENDIF
   80 CONTINUE
C--SEAWAT: CONVERT HB TO EQUIVALENT FRESHWATER HEAD
      HB=FEHEAD(HB,DENSE,ELEV(IC,IR,IL))
   90 CONTINUE
C
C6------UPDATE THE APPROPRIATE HNEW VALUE
      HNEW(IC,IR,IL)=HNEW(IC,IR,IL)+HB
      HOLD(IC,IR,IL)=HNEW(IC,IR,IL)
  100 CONTINUE
C
C7------RETURN
      RETURN
      END
