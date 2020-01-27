*CMZ :          28/05/2004  10.22.45  by  H1 Saclay
*CMZ :  1.06/01 29/06/95  13.43.01  by  Benno List
*CMZ :  1.06/00 16/06/95  11.16.10  by  Benno List
*CMZ :  1.05/01 10/05/95  15.24.34  by  Benno List
*CMZ :  1.00/01 18/07/94  14.58.20  by  Benno List
*CMZ :  0.23/01 07/02/94  20.55.25  by  Benno List
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
      SUBROUTINE DECRHO (IRHO, HELI, IDEC)

      INTEGER IRHO, HELI, IDEC
************************************************************************
*
*     Purpose:    Let a rho decay
*     --------    rho+ -> pi+ pi0
*                 rho- -> pi- pi0
*                 rho0 -> pi+ pi- (99.0%), pi+ pi- gamma (1.0%)
*
*     Input:      IRHO:    Number of rho in /CPROD/
*     ------      HELI:    Helicity of the rho meson
*                 from common block /CPROD/:
*                 NPART, PPCMS8, ITYPE
*
*     Output:     IDEC:    Number of first decay particle in /CPROD/
*     -------     to common block /CPROD/:
*                 NPART, PPCMS8, ITYPE, ISTAT, MOHEP, IDAHEP
*
*     Called by:  DECVM, DECPHI
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

C     PRINT *, '--- DECRHO ----'
C     CALL PR5V8 ('p_rho ', PPCMS8 (1, IRHO))

      IF (ABS (ITYPE (IRHO)) .NE. 213 .AND. ITYPE (IRHO) .NE. 113) THEN
        PRINT *, 'DECRHO error: ITYPE (', IRHO, ') = ', ITYPE (IRHO)
        CALL ERRLOG (160, 'S: DECRHO: Wrong particle!')
        IDEC = 0
        RETURN
      END IF
      IF (ISTAT (IRHO) .NE. 1) THEN
        PRINT *, 'DECRHO error: ISTAT (', IRHO, ') = ', ISTAT (IRHO)
        CALL ERRLOG (161, 'W: DECRHO: rho not undecayed!')
        IDEC = 0
        RETURN
      END IF

*
* Determine parameters for decay angle distribution
*

      IF (HELI .EQ. 0) THEN
        ALPHA = 0.0
        BETA  = 1.0
      ELSE
        ALPHA = 1.0
        BETA  = 0.0
      END IF

*
* Choose decay channel
*

      R = H1RN ()
      IPI0 = 0
      IF (ITYPE (IRHO) .EQ. 113) THEN
        IF (R .LT. 0.010) THEN
          N = 3
          ITYPE (NPART+1) =  211
          ITYPE (NPART+2) = -211
          ITYPE (NPART+3) =   22
        ELSE
          N = 2
          ITYPE (NPART+1) =  211
          ITYPE (NPART+2) = -211
        END IF
      ELSE
        N = 2
        ITYPE (NPART+1) =  SIGN (211, ITYPE (IRHO))
        ITYPE (NPART+2) =  111
        IPI0 = NPART + 2
      END IF

*
* Fill the mass array XM
*

      SUM = 0D0
      DO 10, I = 1, N
        XM (I) = DBLE (PMASS (ITYPE (NPART+I)))
        SUM = SUM + XM (I)
 10   CONTINUE

      IF (SUM .GT. PPCMS8 (5, IRHO)) THEN
        PRINT *, 'DECRHO error: mass ', PPCMS8 (5, IRHO), 'too small!'
        CALL ERRLOG (162, 'S: DECRHO: mass of rho too small!')
        IDEC = 0
        RETURN
      END IF

*
* Iterate RAMBO call until cos (theta*) fits distribution;
* maximally 100 iterations are allowed.
*

      ITER = 0
      PVM2 = PPCMS8 (1, IRHO)**2 + PPCMS8 (2, IRHO)**2 +
     +       PPCMS8 (3, IRHO)**2

 20   CONTINUE
        ITER = ITER + 1

*
* Call RAMBO: distribute particles uniformely in phase space
*

        CALL RAMBO (N, PPCMS8 (5, IRHO), XM, P, WT, 1, IOK)
        IF (IOK .NE. 0) THEN
          CALL ERRLOG (163, 'S: DECRHO: RAMBO failed!')
        END IF
*
* Calculate cos (theta*)
*

        IF (N .EQ. 2) THEN
          CTHST = SNGL (PPCMS8 (1, IRHO)*P (1, 1) +
     +                  PPCMS8 (2, IRHO)*P (2, 1) +
     +                  PPCMS8 (3, IRHO)*P (3, 1)) /
     +          SQRT (SNGL ((P (1, 1)**2+P (2, 1)**2+P (3, 1)**2)*PVM2))

C         PRINT *, 'cos (theta*): ', CTHST

          IF (CTHST .LT. -1.0) THEN
            PRINT *, '*** DECRHO warning: CTHST < -1: ', CTHST
            CALL ERRLOG (164, 'W: DECRHO: CTHST < -1!')
          ELSE IF (CTHST .GT. 1.0) THEN
            PRINT *, '*** DECRHO warning: CTHST > 1: ', CTHST
            CALL ERRLOG (165, 'W: DECRHO: CTHST > -1!')
          END IF
          ACCEPT = (H1RN () .LT. ALPHA + (BETA-ALPHA)*CTHST**2)
        ELSE
          ACCEPT = .TRUE.
        END IF

      IF (.NOT. ACCEPT .AND. ITER .LE. 100) GO TO 20

      IF (ITER .GT. 100) THEN
        PRINT *, '*** DECRHO warning: more than 100 RAMBO iterations'
        CALL ERRLOG (166, 'W: DECRHO: More than 100 RAMBO iterations!')
      END IF

C     PRINT *, 'Iterations: ', ITER

*
* Update /CPROD/
*

      IDEC = NPART + 1

      ISTAT (IRHO) = 2
      IDAHEP (1, IRHO) = NPART + 1
      IDAHEP (2, IRHO) = NPART + N

      DO 50, I = 1, N
        PDPRS (4) = XM (I)**2
        DO 30, J = 1, 3
          PDPRS (J) = P (J, I)
          PDPRS (4) = PDPRS (4) + PDPRS (J)**2
 30     CONTINUE
        PDPRS (4) = SQRT (PDPRS (4))
        CALL LORENB8 (PPCMS8 (5, IRHO), PPCMS8 (1, IRHO), PDPRS, PDPCMS)
        DO 40, J = 1, 4
          PPCMS8 (J, NPART+I) = PDPCMS (J)
 40     CONTINUE
        PPCMS8 (5, NPART+I) = XM (I)
        ISTAT (NPART + I) = 1
        IDAHEP (1, NPART + I) = 0
        IDAHEP (2, NPART + I) = 0
        MOHEP (1, NPART + I) = IRHO
        MOHEP (2, NPART + I) = 0
C       CALL PR4V8 ('PDPRS ', PDPRS)
C       CALL PR4V8 ('PDPCMS', PDPCMS)
C       CALL PR5V8 ('PPCMS8', PPCMS8 (1, NPART+I))
 50   CONTINUE

      NPART = NPART + N

*
* Perform pi0 decay
*

      IF (IPI0 .GT. 0) THEN
        CALL DECPI0 (IPI0)
      END IF

      RETURN
      END
