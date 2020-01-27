*CMZ :          28/05/2004  10.22.25  by  H1 Saclay
*CMZ :  1.06/01 29/06/95  13.48.52  by  Benno List
*CMZ :  1.06/00 16/06/95  10.55.56  by  Benno List
*CMZ :  0.23/01 07/02/94  15.10.15  by  Benno List
*CMZ :  0.23/00 03/02/94  14.34.14  by  Benno List
*-- Author :    Benno List   03/02/94
************************************************************************
      SUBROUTINE DECPI0 (IPI0)

      INTEGER IPI0
************************************************************************
*
*     Purpose:   Let a pi0 decay isotropically into photons
*     --------
*
*     Input:     IPI0:   Number of pi0 in /CPROD/
*                from common block /CPROD/:
*     ------     NPART, PPCMS8, ITYPE, ISTAT
*
*     Output:    to common block /CPROD/
*     -------    NPART, PPCMS8, MOHEP, IDAHEP, ITYPE, ISTAT
*
*     Author:    Benno List  3.2.94
*     -------
*
*     Called by: FRAGPX, FRAGVX, CASCAD, DECVM
*     ----------
*
*     Calling:   LORENB8
*     --------   From H1UTIL:
*                H1RN
*
************************************************************************

      include 'cprod.common'

      REAL*8 PDPRS (4), CTHETA, PHI, E
      LOGICAL I

      IF (ITYPE (IPI0) .NE. 111) THEN
        PRINT *, 'DECPI0 error: ITYPE (', IPI0, ') = ', ITYPE (IPI0)
        CALL ERRLOG (140, 'S: DECPI0: Wrong particle')
        RETURN
      END IF
      IF (ISTAT (IPI0) .NE. 1) THEN
        PRINT *, 'DECPI0 error: ISTAT (', IPI0, ') = ', ISTAT (IPI0)
        CALL ERRLOG (141, 'W: DECPI0: pi0 not undecayed')
        RETURN
      END IF

      ISTAT (IPI0) = 2
      IDAHEP (1, IPI0) = NPART+1
      IDAHEP (2, IPI0) = NPART+2
      CTHETA = 1.99998D0*H1RN () - 0.99999D0
      PHI = 6.28318530718D0 * H1RN ()
      E = 0.5D0 * PPCMS8 (5, IPI0)
      PDPRS (1) = E*COS (PHI)*CTHETA
      PDPRS (2) = E*SIN (PHI)*CTHETA
      PDPRS (3) = E*SQRT (1D0 - CTHETA**2)
      PDPRS (4) = E
      N = NPART+1
      CALL LORENB8 (PPCMS8(5,IPI0), PPCMS8(1,IPI0), PDPRS, PPCMS8(1,N))
      PPCMS8 (5, N) = 0D0
      ISTAT (N) = 1
      ITYPE (N) = 22
      IDAHEP (1, N) = 0
      IDAHEP (2, N) = 0
      MOHEP (1, N) = IPI0
      MOHEP (2, N) = 0
      PDPRS (1) = -PDPRS (1)
      PDPRS (2) = -PDPRS (2)
      PDPRS (3) = -PDPRS (3)
      N = N+1
      CALL LORENB8 (PPCMS8(5,IPI0), PPCMS8(1,IPI0), PDPRS, PPCMS8(1,N))
      PPCMS8 (5, N) = 0D0
      ISTAT (N) = 1
      ITYPE (N) = 22
      IDAHEP (1, N) = 0
      IDAHEP (2, N) = 0
      MOHEP (1, N) = IPI0
      MOHEP (2, N) = 0

      NPART = N

      RETURN
      END
