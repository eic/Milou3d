*CMZ :          28/05/2004  10.21.59  by  H1 Saclay
*CMZ :  1.06/00 16/06/95  11.14.33  by  Benno List
*CMZ :  1.01/01 05/09/94  14.48.47  by  Benno List
*CMZ :  1.01/00 01/09/94  13.40.14  by  Benno List
*-- Author :    Benno List   01/09/94
************************************************************************
      SUBROUTINE DECDEL (IDELTA, IDEC)

      INTEGER IDELTA, IDEC
************************************************************************
*
*     Purpose:    Let a Delta(1232) decay
*     --------    Delta++ -> p pi+ (100%)
*                 Delta+  -> p pi0 (67%), n pi+ (33%)
*                 Delta0  -> n pi0 (67%), p pi- (33%)
*                 Delta-  -> n pi- (100%)
*
*     Input:      IDELTA:  Number of Delta in /CPROD/
*     ------      from common block /CPROD/:
*                 NPART, PPCMS8, ITYPE
*
*     Output:     IDEC:    Number of first decay particle in /CPROD/
*     -------     to common block /CPROD/:
*                 NPART, PPCMS8, ITYPE, ISTAT, MOHEP, IDAHEP
*
*     Called by:  DECNST
*     ----------
*
*     Author:     Benno List, 1.9.94
*     -------
*
*     Changed by:
*     -----------
*
*     Calling:  From CERN Library:
*     --------  U102 LORENB8 LORENTZ boost back to lab system (REAL*8)
*
*               From H1UTIL:
*               PMASS
*               H1RN         Random number generator
*
*               RAMBO        RAndom Momenta BOoster (phase space gen.)
*               DECPI0       pi0 decay
*
************************************************************************

      include 'cprod.common'

      REAL*8 XM (100), P (4, 100), WT, SUM
      REAL*8 PDPCMS (1:4), PDPRS (1:4)

C     PRINT *, '--- DECDEL ----'
C     CALL PR5V8 ('p_Del ', PPCMS8 (1, IDELTA))
C     PRINT *, ITYPE (IDELTA), ' decays to'

      IF (ISTAT (IDELTA) .NE. 1) THEN
        PRINT *, 'DECDEL error: ISTAT (', IDELTA, ') = ', ISTAT (IDELTA)
        CALL ERRLOG (210, 'W: DECDEL: Particle not undecayed!')
        IDEC = 0
        RETURN
      END IF

*
* Fill ITYPE
*

      IPI0 = 0
      ISIG = SIGN (1, ITYPE (IDELTA))
C     PRINT *, 'ISIG = ', ISIG
      IF (ABS (ITYPE (IDELTA)) .EQ. 2224) THEN
        ITYPE (NPART+1) = 2212*ISIG
        ITYPE (NPART+2) =  211
      ELSE IF (ABS (ITYPE (IDELTA)) .EQ. 2214) THEN
        IF (H1RN () .LE. 0.333) THEN
          ITYPE (NPART+1) = 2112*ISIG
          ITYPE (NPART+2) =  211*ISIG
        ELSE
          ITYPE (NPART+1) = 2212*ISIG
          ITYPE (NPART+2) =  111
        END IF
      ELSE IF (ABS (ITYPE (IDELTA)) .EQ. 2114) THEN
        IF (H1RN () .LE. 0.333) THEN
          ITYPE (NPART+1) = 2212*ISIG
          ITYPE (NPART+2) = -211*ISIG
        ELSE
          ITYPE (NPART+1) = 2112*ISIG
          ITYPE (NPART+2) =  111
        END IF
      ELSE IF (ABS (ITYPE (IDELTA)) .EQ. 1114) THEN
        ITYPE (NPART+1) = 2112*ISIG
        ITYPE (NPART+2) = -211*ISIG
      ELSE
        IDEC = 0
        PRINT *, '*** DECDEL Error: ', ITYPE (IDELTA), ' is not a Delta'
        CALL ERRLOG (211, 'S: DECDEL: Wrong particle!')
        RETURN
      END IF

C     PRINT *, ITYPE (NPART+1), ITYPE (NPART+2)

*
* Fill the mass array XM
*

      SUM = 0D0
      DO 10, I = 1, 2
        XM (I) = DBLE (PMASS (ITYPE (NPART+I)))
        SUM = SUM + XM (I)
 10   CONTINUE

      IF (SUM .GT. PPCMS8 (5, IDELTA)) THEN
        PRINT *, 'DECDEL error: mass ', PPCMS8 (5, IDELTA), 'too small!'
        CALL ERRLOG (212, 'S: DECDEL: Mass of Delta too small!')
        IDEC = 0
        RETURN
      END IF

*
* Call RAMBO: distribute particles uniformely in phase space
*

      CALL RAMBO (2, PPCMS8 (5, IDELTA), XM, P, WT, 1, IOK)

      IF (IOK .NE. 0) THEN
        CALL ERRLOG (213, 'S: DECDEL: RAMBO failed!')
      END IF

*
* Update /CPROD/
*

      IDEC = NPART + 1

      ISTAT (IDELTA) = 2
      IDAHEP (1, IDELTA) = NPART + 1
      IDAHEP (2, IDELTA) = NPART + 2

      DO 50, I = 1, 2
        PDPRS (4) = XM (I)**2
        DO 30, J = 1, 3
          PDPRS (J) = P (J, I)
          PDPRS (4) = PDPRS (4) + PDPRS (J)**2
 30     CONTINUE
        PDPRS (4) = SQRT (PDPRS (4))
        CALL LORENB8 (PPCMS8(5,IDELTA), PPCMS8 (1,IDELTA),PDPRS, PDPCMS)
        DO 40, J = 1, 4
          PPCMS8 (J, NPART+I) = PDPCMS (J)
 40     CONTINUE
        PPCMS8 (5, NPART+I) = XM (I)
        ISTAT (NPART + I) = 1
        IDAHEP (1, NPART + I) = 0
        IDAHEP (2, NPART + I) = 0
        MOHEP (1, NPART + I) = IDELTA
        MOHEP (2, NPART + I) = 0
C       CALL PR4V8 ('PDPRS ', PDPRS)
C       CALL PR4V8 ('PDPCMS', PDPCMS)
C       CALL PR5V8 ('PPCMS8', PPCMS8 (1, NPART+I))
 50   CONTINUE

      NPART = NPART + 2

*
* Perform pi0 decay
*

      IF (ITYPE (IDEC+1) .EQ. 111) THEN
        CALL DECPI0 (NPART)
      END IF

      RETURN
      END
