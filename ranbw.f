*CMZ :  0.23/00 28/01/94  13.00.51  by  Benno List
*CMZ :  0.19/01 20/07/93  10.19.53  by  Benno List
*CMZ :  0.18/00 18/07/93  16.50.05  by  Benno List
*-- Author :    Benno List   18/07/93
************************************************************************
      REAL FUNCTION RANBW (ER, GAMMA, EMIN, EMAX)

      REAL*8 ER, GAMMA, EMIN, EMAX
************************************************************************
*
*     Purpose:   Generate random number with BREIT-WIGNER distribution
*     --------
*
*     Input:     ER:     Maximum of distribution
*     ------     GAMMA:  Width of distribution
*                EMIN:   Minimal value of RANBW
*                EMAX:   Maximal value of RANBW
*
*     Output:    Function value: Random number between EMIN and EMAX
*     -------    with BREIT-WIGNER distribution:
*                1 / ((E - ER)**2 + GAMMA**2/4)
*
*     Author:    B. List,  18.7.93
*     -------
*
************************************************************************

      IF (GAMMA .LT. 1.0E-3*ER) THEN
        RANBW = ER
      ELSE
        A = ATAN (2.0*(EMAX - ER)/GAMMA)
        B = ATAN (2.0*(EMIN - ER)/GAMMA)
        E = ER + 0.5*GAMMA*TAN (H1RN ()*(A - B) + B)
        IF (E .LE. EMAX) THEN
          RANBW = E
        ELSE
          RANBW = EMAX
        END IF
      END IF

      RETURN
      END
