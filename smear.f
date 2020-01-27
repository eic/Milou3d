CC-----------------------------
      SUBROUTINE SMEAR(SM1,SIG)
CC-----------------------------

      IMPLICIT NONE

      REAL SM1, SIG, RNDM
      INTEGER I

      DO I = 1, 12
         SM1 = SM1 + RNDM(1.)
      ENDDO

      SM1 = SIG*(SM1-6.0)

      RETURN
      END

