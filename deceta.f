*CMZ :          28/05/2004  10.24.31  by  H1 Saclay
*CMZ :  1.06/00 16/06/95  11.18.45  by  Benno List
*CMZ :  0.23/01 07/02/94  17.27.15  by  Benno List
*CMZ :  0.23/00 03/02/94  11.31.56  by  Benno List
*CMZ :  0.20/00 11/10/93  11.23.32  by  Benno List
*CMZ :  0.13/00 15/06/93  10.12.04  by  Benno List
*CMZ :  0.08/00 20/05/93  14.01.54  by  Benno List
*CMZ :  0.07/00 20/05/93  12.54.35  by  Benno List
*CMZ :  0.06/01 20/05/93  12.35.12  by  Benno List
*CMZ :  0.06/00 19/05/93  15.48.55  by  Benno List
*CMZ :  0.05/06 18/05/93  15.23.31  by  Benno List
*CMZ :  0.00/04 18/02/93  13.57.37  by  Benno List
*CMZ :  0.00/00 15/02/93  17.17.14  by  Benno List
*-- Author :    Benno List   25/01/93
************************************************************************
      SUBROUTINE DECETA (IETA)

      INTEGER IETA
************************************************************************
*
*     Purpose:    Let a eta0 decay
*     --------    eta0 -> 2 gamma (39.1%), 3 pi0 (32.2%),
*                           pi+ pi- pi0 (23.8%), pi+ pi- gamma (4.9%)
*
*     Input:      IETA:    Number of eta0 in /CPROD/
*     ------      from common block /CPROD/:
*                 NPART, PPCMS8, ITYPE
*
*     Output:     to common block /CPROD/:
*     -------     NPART, PPCMS8, ITYPE, ISTAT, MOHEP, IDAHEP
*
*     Called by:  DECPHI
*     ----------
*
*     Author:     Benno List, 7.1.94, 25.1.93 (DECAYVM)
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
      LOGICAL ACCEPT

C     PRINT *, '--- DECETA ----'
C     CALL PR5V8 ('p_eta0', PPCMS8 (1, IETA))

      IF (ITYPE (IETA) .NE. 221) THEN
        PRINT *, 'DECETA error: ITYPE (', IETA, ') = ', ITYPE (IETA)
        CALL ERRLOG (190, 'S: DECETA: Wrong particle!')
        IDEC = 0
        RETURN
      END IF
      IF (ISTAT (IETA) .NE. 1) THEN
        PRINT *, 'DECETA error: ISTAT (', IETA, ') = ', ISTAT (IETA)
        CALL ERRLOG (191, 'W: DECETA: eta not undecayed!')
        IDEC = 0
        RETURN
      END IF

*
* Choose decay channel
*

      R = H1RN ()
      IF (R .LT. 0.049) THEN
        N = 3
        ITYPE (NPART+1) =  211
        ITYPE (NPART+2) = -211
        ITYPE (NPART+3) =   22
      ELSE IF (R .LT. 0.287) THEN
        N = 3
        ITYPE (NPART+1) =  211
        ITYPE (NPART+2) = -211
        ITYPE (NPART+3) =  111
      ELSE IF (R .LT. 0.287) THEN
        N = 3
        ITYPE (NPART+1) =  111
        ITYPE (NPART+2) =  111
        ITYPE (NPART+3) =  111
      ELSE
        N = 2
        ITYPE (NPART+1) =   22
        ITYPE (NPART+2) =   22
      END IF

*
* Fill the mass array XM
*

      SUM = 0D0
      DO 10, I = 1, N
        XM (I) = DBLE (PMASS (ITYPE (NPART+I)))
        SUM = SUM + XM (I)
 10   CONTINUE

      IF (SUM .GT. PPCMS8 (5, IETA)) THEN
        PRINT *, 'DECETA error: mass ', PPCMS8 (5, IETA), 'too small!'
        CALL ERRLOG (192, 'S: DECETA: mass of eta too small!')
        IDEC = 0
        RETURN
      END IF

*
* Iterate RAMBO call until cos (theta*) fits distribution;
* maximally 100 iterations are allowed.
*

      ITER = 0

 20   CONTINUE
        ITER = ITER + 1

*
* Call RAMBO: distribute particles uniformely in phase space
*

        CALL RAMBO (N, PPCMS8 (5, IETA), XM, P, WT, 1, IOK)
        IF (IOK .NE. 0) THEN
          CALL ERRLOG (193, 'S: DECETA: RAMBO failed!')
        END IF
*
* Calculate cos (theta*)
*

        IF (N .EQ. 2) THEN
          ACCEPT = .TRUE.
        ELSE
          ACCEPT = .TRUE.
        END IF

      IF (.NOT. ACCEPT .AND. ITER .LE. 100) GO TO 20

      IF (ITER .GT. 100) THEN
        PRINT *, '*** DECETA warning: more than 100 RAMBO iterations'
        CALL ERRLOG (194, 'S: DECETA: mass of eta too small!')
      END IF

C     PRINT *, 'Iterations: ', ITER

*
* Update /CPROD/
*

      IDEC = NPART + 1

      ISTAT (IETA) = 2
      IDAHEP (1, IETA) = NPART + 1
      IDAHEP (2, IETA) = NPART + N

      DO 50, I = 1, N
        PDPRS (4) = XM (I)**2
        DO 30, J = 1, 3
          PDPRS (J) = P (J, I)
          PDPRS (4) = PDPRS (4) + PDPRS (J)**2
 30     CONTINUE
        PDPRS (4) = SQRT (PDPRS (4))
        CALL LORENB8 (PPCMS8 (5, IETA), PPCMS8 (1, IETA), PDPRS, PDPCMS)
        DO 40, J = 1, 4
          PPCMS8 (J, NPART+I) = PDPCMS (J)
 40     CONTINUE
        PPCMS8 (5, NPART+I) = XM (I)
        ISTAT (NPART + I) = 1
        IDAHEP (1, NPART + I) = 0
        IDAHEP (2, NPART + I) = 0
        MOHEP (1, NPART + I) = IETA
        MOHEP (2, NPART + I) = 0
C       CALL PR4V8 ('PDPRS ', PDPRS)
C       CALL PR4V8 ('PDPCMS', PDPCMS)
C       CALL PR5V8 ('PPCMS8', PPCMS8 (1, NPART+I))
 50   CONTINUE

      NPART = NPART + N

*
* Perform pi0 decay
*
      N1 = NPART - N + 1
      N2 = NPART
      DO 60, I = N1, N2
        IF (ITYPE (I) .EQ. 111) THEN
          CALL DECPI0 (I)
        END IF
 60   CONTINUE

      RETURN
      END
