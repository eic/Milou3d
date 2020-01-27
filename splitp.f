*CMZ :          30/05/2004  07.50.33  by  H1 Saclay
*CMZ :  1.06/02 03/04/96  21.30.13  by  Benno List
*CMZ :  1.06/00 16/06/95  10.51.03  by  Benno List
*CMZ :  0.23/00 31/01/94  14.49.55  by  Benno List
*CMZ :  0.22/00 25/01/94  15.13.43  by  Benno List
*-- Author :    Benno List   20/01/94
************************************************************************
      SUBROUTINE SPLITP (P, Q, PQ, PREMN, IERR)

      REAL*8 P (5), Q (5), PQ (5), PREMN (5)
      INTEGER IERR
************************************************************************
*
*     Purpose:  A quark in a hadron with 5-momentum P is struck
*     --------  by a boson with 5-momentum Q;
*               Q is to be negative for spacelike particles.
*               The hadron splits into
*               a remnant with 5-momentum (1-x)*P(5) and
*               a quark with mass P (5) and 3-momentum x*P(1:3)+Q(1:3).
*               x is determined so that 4-vector conservation holds.
*
*     Input:    P:         5-momentum of the hadron
*     ------    Q:         5-momentum of the boson
*
*     Output:   PQ:        5-vector of the quark
*     -------   PREMN:     5-vector of the remnant
*               IERR:      0 if everything is O.K.,
*                         -1 if x<0 or x>1.
*                          1 if no solution was found.
*
*     Author:   Benno List, 20.1.94
*     -------
*
************************************************************************

      REAL*8 X, D, P2, Q2

      real*8 x_main,q_main,phi_main,t_main,ym_main
      common/dvcs_VAR/x_main,q_main,phi_main,t_main,ym_main


      P2= P (5)*ABS (P (5))
      Q2= Q (5)*ABS (Q (5))
      D = (P(4)+Q(4))**2-(P(1)+Q(1))**2-(P(2)+Q(2))**2-(P(3)+Q(3))**2
     +     - P2 - Q2
      IF (D .EQ. 0D0) THEN
        IERR = 1
        PRINT *, '### SPLITP error: D = 0!'
        CALL ERRLOG (130, 'S: SPLITP: D = 0!')
        PRINT *, 'P = (', P(1),',',P(2),',',P(3),',',P(4),'), ',P(5)
        PRINT *, 'Q = (', Q(1),',',Q(2),',',Q(3),',',Q(4),'), ',Q(5)
      ELSE
        X = -Q2 / D
        IF (X .LT. 0D0 .OR. X .GT. 1D0) THEN
          IERR = -1
          CALL ERRLOG (131, 'W: SPLITP: x<0 or x>1')
        ELSE
          IERR = 0
        END IF

cep << --
cep      write(6,*) 'x x_main ',x,x_main
cep      x = x_main
cep -- >>



        DO 10, I = 1, 4
          PREMN (I) = (1D0 - X)*P (I)
 10     CONTINUE
        PREMN (5) = ABS (1D0 - X)*P (5)
        PQ (5) = ABS (X)*P (5)
        PQ (4) = X**2 * P2
        DO 20, I = 1, 3
          PQ (I) = X*P (I) + Q (I)
          PQ (4) = PQ (4) + PQ (I)**2
 20     CONTINUE
        IF (PQ (4) .LT. 0D0) THEN
          IERR = 1
        END IF
        PQ (4) = SQRT (ABS (PQ (4)))

C       PRINT *, '   SPLITP: X =             ', X
C       PRINT *, '   SPLITP: (P+Q)**2 =      ', (P (4) + Q (4))**2
C    +      -(P (1)+Q (1))**2-(P (2)+Q (2))**2-(P (3)+Q (3))**2
C       PRINT *, '   SPLITP: (PQ+PREMN)**2 = ',(PQ (4) + PREMN (4))**2
C    +      -(PQ(1)+PREMN(1))**2-(PQ(2)+PREMN(2))**2-(PQ(3)+PREMN(3))**2
C       PRINT *, '   SPLITP: PQ (5) =        ', PQ (5)
C       PRINT *, '   SPLITP: PREMN (5) =     ', PREMN (5)
      END IF

      RETURN
      END
