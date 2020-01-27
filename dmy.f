*CMZ :          27/05/2004  18.18.31  by  H1 Saclay
*-- Author :    H1 Saclay   27/05/2004



************************************************************************
      REAL*8 FUNCTION DMY(ym)
************************************************************************
*
*     Purpose:    Generate hadronic mass between Mmin and Mmax
*     --------    (mass distribution may depend on t)for p-vertex
*
*     Superseds:  MASSREJ1 by R. LAUSEN
*     ----------
*
*     Input:      MMIN:   minimal allowed mass
*     ------      MMAX:   maximal allowed mass
*                 T:      MANDELSTAM variable t
*                 from common block GDIF: EPSILM
*
*     Output:     MPX:    hadronic mass in GeV
*     -------     Y:      Y > 1. corresponds to an enhancement of
*                     the mass spectrum over the 1/mx**2 spectrum
*
*     Called by:  GENVM
*     ----------
*
*     Author:     Benno List, 14.1.92
*     -------
*
*     Changed by:
*     -----------
*
*     Calling:    From H1UTIL:
*     --------
*                 H1RN:  Random generator
*
*     Remarks:    T is up to now unused
*     --------
*
************************************************************************

      include 'dvcs.common'

      REAL*8 MPX, MMIN, MMAX, T, Y, ENHA

      REAL*8 M2, LMIN, DELTA, M2MIN, FACT
      REAL*8 ym

      EPSILM = stEPSM
      MMIN   = dble((stMYMIN))
      MMAX   = dble((stMYMAX))

      IF (ABS (EPSILM) .LT. 0.001) THEN
        LMIN = 2.0D0*DLOG (MMIN)
        DELTA = 2.0D0*DLOG (MMAX/MMIN)
      ELSE
        M2MIN = MMIN**(-2*EPSILM)
        FACT = (MMAX**(-2*EPSILM) - M2MIN)
      END IF

      ITER = 0

ccc 10   CONTINUE
ccc        ITER = ITER + 1

ccc        IF (ABS (EPSILM) .LT. 0.001) THEN
ccc* Basic specrum: 1/M^2
ccc          M2 = DEXP (DBLE (H1RN ())*DELTA + LMIN)
ccc        ELSE
ccc* Basic spectrum: 1/M^2(1+epsilon)
ccc          M2 = (FACT*H1RN () + M2MIN)**(-1.0/EPSILM)
ccc        END IF


        M2 = ym**2

        IF (M2 .LT. MMIN**2) THEN
          DMY = 0.d0
          goto 999
        ELSE IF (M2 .GT. MMAX**2) THEN
          DMY = 0.d0
          goto 999
        END IF
*
*       Old version with enhancements in lower mass region
*
        IF (M2 .GE. 4.00D0) THEN
          Y = 1.0D0
        ELSE IF (M2 .GE. 3.10D0) THEN
          Y = 1.64D0 - 0.16*M2
        ELSE IF (M2 .GE. 2.65D0) THEN
          Y = M2*(0.47D0 - 0.42D0*(M2-2.65D0)**2)
        ELSE IF (M2 .GE. 2.25D0) THEN
          Y = M2*(0.47D0 + 0.46D0*(M2-2.65D0)**2)
        ELSE IF (M2 .GE. 2.02D0) THEN
          Y = M2*(0.76D0 - 2.69D0*(M2-2.02D0)**2)
        ELSE IF (M2 .GE. 1.72D0) THEN
          Y = M2*(0.76D0 - 1.98D0*(M2-2.02D0)**2)
        ELSE
          Y = 1.05*(M2 - 1.165D0)
        END IF

       DMY = Y/(ym**(2*(1.+EPSILM)))

       xJACOB = 2.d0*ym

       DMY = DMY/xJACOB

 999   continue

       write(67,*) ym**2,DMY


      RETURN
      END

