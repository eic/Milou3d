*CMZ :          29/05/2004  19.17.17  by  H1 Saclay
*CMZ :  1.06/00 16/06/95  11.14.57  by  Benno List
*CMZ :  1.03/00 13/02/95  15.13.20  by  Benno List
*CMZ :  1.02/00 03/11/94  10.00.50  by  Benno List
*CMZ :  1.01/01 21/10/94  16.21.23  by  Benno List
*CMZ :  1.01/00 01/09/94  17.56.16  by  Benno List
*-- Author :    Benno List   01/09/94
************************************************************************
      SUBROUTINE DECNST (INSTAR, IDEC, IDECN)

      include 'dvcs.common'
      include 'cprod.common'

      INTEGER INSTAR, IDEC
************************************************************************
*
*     Purpose:    Let a N* decay
*     --------    Possible decay modes (also charge conjugates)
*                 N*+ -> p pi0, n pi+
*                 N*+ -> Delta++ pi-, Delta+ pi0, Delta0 pi+
*                 N*+ -> p (pi+ pi-)_I=0, p (pi0 pi0)_I=0
*                 N*+ -> p rho0, n rho+
*                 N*+ -> p eta
*                 N*+ -> Lambda K+
*                 N*0 -> p pi-, n pi0
*                 N*0 -> Delta+ pi-, Delta0 pi0, Delta- pi+
*                 N*0 -> n (pi+ pi-)_I=0, n (pi0 pi0)_I=0
*                 N*0 -> p rho-, n rho0
*                 N*0 -> n eta
*                 N*0 -> Lambda K0S, Lambda K0L
*                 Branching ratios between lines taken from
*                   Rev. Part. Prop. 1994: Phys. Rev. D50, Pt. I, 1173.
*                 Branching ratios between particles in one line
*                 calculated from CLEBSCH-GORDAN coefficients
*                 N* can be:
*                 name     PDG code      name     PDG code
*                 N(1440)+    12212  N(1440)0    12112
*                 N(1520)+     2124  N(1520)0     1214
*                 N(1535)+    22212  N(1535)0    22112
*                 N(1650)+    32212  N(1650)0    32112
*                 N(1675)+     2216  N(1675)0     2116
*                 N(1680)+    12216  N(1680)0    12116
*                 N(1700)+    22124  N(1700)0    21214
*                 N(1710)+    42212  N(1710)0    42112
*                 N(1720)+    32124  N(1720)0    31214
*
*     Input:      INSTAR:  Number of N* in /CPROD/
*     ------      from common block /CPROD/:
*                 NPART, PPCMS8, ITYPE
*
*     Output:     IDEC:    Number of first decay particle in /CPROD/
*     -------     IDECN:   Number of stable decay nucleon in /CPROD/
*                 to common block /CPROD/:
*                 NPART, PPCMS8, ITYPE, ISTAT, MOHEP, IDAHEP
*
*     Called by:  GENEVT
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
*               DECRHO       rho decay
*               DECDEL       Delta decay
*               DECETA       eta decay
*
************************************************************************

      PARAMETER(NMAX = 4000)

      COMMON /LUJETS/ NLUJETS,KLU(4000,5),PLU(4000,5),VLU(4000,5)
      integer benno_code


      PARAMETER (NCODES = 9)
      PARAMETER (NDECAY = 6)

* Legal N*+ codes
      INTEGER NPLCOD (1:NCODES)
     +        /12212, 2124, 22212, 32212, 2216, 12216, 22124,
     +         42212, 32124/
* Legal N*0 codes
      INTEGER N0COD (1:NCODES)
     +        /12112, 1214, 22112, 32112, 2116, 12116, 21214,
     +         42112, 31214/
* Cummulative branching ratios for the channels
* N* -> N pi, Delta pi, N rho, N (pi pi)_I=0, N eta, Lambda K
* for the different mass states
      REAL DECCHN (1:NDECAY-1, 1:NCODES)
     +        /0.67, 0.92, 0.92, 1.00, 1.00,
     +         0.60, 0.80, 1.00, 1.00, 1.00,
     +         0.50, 0.50, 0.50, 0.50, 1.00,
     +         0.75, 0.80, 0.90, 0.90, 0.90,
     +         0.45, 1.00, 1.00, 1.00, 1.00,
     +         0.67, 0.78, 0.89, 1.00, 1.00,
     +         0.10, 0.45, 0.65, 1.00, 1.00,
     +         0.15, 0.40, 0.55, 0.85, 0.85,
     +         0.15, 0.15, 0.90, 0.90, 0.90/
* Cummulative branching ratios for different particle charges
      REAL CGCOEF (1:2, 1:NDECAY)
     +        /0.33, 1.00,
     +         0.50, 0.83,
     +         0.33, 1.00,
     +         0.67, 1.00,
     +         1.00, 1.00,
     +         1.00, 1.00/
* Decay particles for N*+ decays
      INTEGER  NPLDC1 (1:3, 1:NDECAY), NPLDC2 (1:3, 1:NDECAY),
     +         NPLDC3 (1:3, 1:NDECAY)
      DATA ((NPLDC1(I,J), NPLDC2(I,J), NPLDC3(I,J), I=1, 3), J=1,NDECAY)
* N* -> N pi
     +        /2212,  111,    0,
     +         2112,  211,    0,
     +            0,    0,    0,
* N* -> Delta pi
     +         2224, -211,    0,
     +         2214,  111,    0,
     +         2114,  211,    0,
* N* -> N rho
     +         2212,  113,    0,
     +         2112,  213,    0,
     +            0,    0,    0,
* N** -> N (pi pi)_I=0
     +         2212,  211, -211,
     +         2212,  111,  111,
     +            0,    0,    0,
* N* -> N eta
     +         2212,  221,    0,
     +            0,    0,    0,
     +            0,    0,    0,
* N* -> Lambda K
     +         3122,  321,    0,
     +            0,    0,    0,
     +            0,    0,    0/

* Decay particles for N*0 decays
      INTEGER N0DC1 (1:3, 1:NDECAY), N0DC2 (1:3, 1:NDECAY),
     +        N0DC3 (1:3, 1:NDECAY)
      DATA ((N0DC1(I,J), N0DC2(I,J), N0DC3(I,J), I=1, 3), J=1, NDECAY)
* N* -> N pi
     +        /2112,  111,    0,
     +         2212, -211,    0,
     +            0,    0,    0,
* N* -> Delta pi
     +         1114,  211,    0,
     +         2114,  111,    0,
     +         2214, -211,    0,
* N* -> N rho
     +         2112,  113,    0,
     +         2212, -213,    0,
     +            0,    0,    0,
* N** -> N (pi pi)_I=0
     +         2112,  211, -211,
     +         2112,  111,  111,
     +            0,    0,    0,
* N* -> N eta
     +         2112,  221,    0,
     +            0,    0,    0,
     +            0,    0,    0,
* N* -> Lambda K
     +         3122,  311,    0,
     +            0,    0,    0,
     +            0,    0,    0/



      REAL*8 XM (100), P (4, 100), WT, SUM, WIDTH, XMMIN
      REAL*8 PDPCMS (1:4), PDPRS (1:4)

C     PRINT *, '--- DECNST ----'
C     PRINT *, 'ITYPE (INSTAR) = ', ITYPE (INSTAR)
C     CALL PR5V8 ('p_N*  ', PPCMS8 (1, INSTAR))


      if (stIDEBUG.ge.2) then
       write(6,*) ' --   begin of DECNST'
      endif


c --- Fill the CPROD common :   -----------------------------------------

      NPART = 1
      if (stIRAD.eq.0) then
       iline = 7
      else
       iline = 9
      endif
      PPART(1,1) = PLU(iline,1)
      PPART(2,1) = PLU(iline,2)
      PPART(3,1) = PLU(iline,3)
      PPART(4,1) = PLU(iline,4)
      PPART(5,1) = PLU(iline,5)
      do i=1,5
       PPART8(i,1) = dble(PPART(i,1))
      enddo
c ---    for the while, in the resonance rest frame
      PPCMS8(1,1) = PPART8(1,1)
      PPCMS8(2,1) = PPART8(2,1)
      PPCMS8(3,1) = PPART8(3,1)
      PPCMS8(4,1) = PPART8(4,1)
      PPCMS8(5,1) = PPART8(5,1)
cep      ITYPE(1) = KLU(iline,2)
      ITYPE(1) = benno_code (KLU(iline,2))
      MOHEP(1,1) = 0
      MOHEP(2,1) = 0
      IDAHEP(1,1) = 0
      IDAHEP(2,1) = 0
      ISTAT(1) = 1


       do i=1,NPARTMAX
        hasdecayed(i) = 0
       enddo

c -----------------------------------------------------------------------
      IF (ISTAT (INSTAR) .NE. 1) THEN
        PRINT *, 'DECNST error: ISTAT (', INSTAR, ') = ', ISTAT (INSTAR)
        CALL ERRLOG (200, 'W: DECNST: Particle not undecayed!')
        IDEC  = 0
        IDECN = 0
        RETURN
      END IF

*
* Look up particle code
*
      NNSTAR = 0
      NCHARG = 0
      DO 1, I = 1, NCODES
        IF (ABS (ITYPE (INSTAR)) .EQ. NPLCOD (I)) THEN
          NNSTAR = I
          NCHARG = 1
        ELSE IF (ABS (ITYPE (INSTAR)) .EQ. N0COD (I)) THEN
          NNSTAR = I
          NCHARG = 0
        END IF
 1    CONTINUE

      if (stIDEBUG.ge.2) then
       write(6,*) 'NNSTAR = ', NNSTAR, ', NCHARG = ', NCHARG
      endif

      IF (NNSTAR .EQ. 0) THEN
        PRINT *, '*** DECNST error: Don''t know particle ',
     +           ITYPE (INSTAR)
        CALL ERRLOG (201, 'S: DECNST: Wrong particle!')
        IDEC  = 0
        IDECN = 0
        RETURN
      END IF

*
* Choose decay channel
*
 7    CONTINUE

      R = H1RN ()
      NDEC1 = 1
 2    CONTINUE
        IF (NDEC1 .EQ. NDECAY .OR. DECCHN (NDEC1, NNSTAR) .GE. R) GOTO 3
        NDEC1 = NDEC1 + 1
      GO TO 2
 3    CONTINUE

C     IF (NDEC1 .LT. NDECAY) THEN
C       PRINT *, 'R = ', R, ', DECCHN (', NDEC1, ', ', NNSTAR, ') = ',
C    +           DECCHN (NDEC1, NNSTAR)
C     ELSE
C       PRINT *, 'R = ', R, ', DECCHN (', NDEC1-1, ', ', NNSTAR, ') = ',
C    +           DECCHN (NDEC1-1, NNSTAR)
C     END IF

      R = H1RN ()
      NDEC2 = 1
 4    CONTINUE
        IF (NDEC2 .EQ. 3 .OR. CGCOEF (NDEC2, NDEC1) .GE. R) GO TO 5
        NDEC2 = NDEC2 + 1
      GO TO 4
 5    CONTINUE

      IF (NDEC2 .LT. 3) THEN
C       PRINT *, 'R = ', R, ', CGCOEF (', NDEC2, ', ', NDEC1, ') = ',
C    +           CGCOEF (NDEC2, NDEC1)
      ELSE
C       PRINT *, 'R = ', R, ', CGCOEF (', NDEC2-1, ', ', NDEC1, ') = ',
C    +           CGCOEF (NDEC2-1, NDEC1)
      END IF


*
* Fill the mass array XM
*

      IDEC  = NPART + 1
      IDECN = NPART + 1
      IF (NCHARG .EQ. 1) THEN
        ITYPE (IDEC)   = NPLDC1 (NDEC2, NDEC1)
        ITYPE (IDEC+1) = NPLDC2 (NDEC2, NDEC1)
        ITYPE (IDEC+2) = NPLDC3 (NDEC2, NDEC1)
      ELSE
        ITYPE (IDEC)   = N0DC1 (NDEC2, NDEC1)
        ITYPE (IDEC+1) = N0DC2 (NDEC2, NDEC1)
        ITYPE (IDEC+2) = N0DC3 (NDEC2, NDEC1)
      END IF
      IF (ITYPE (INSTAR) .LT. 0) THEN
        ITYPE (IDEC) = -ITYPE (IDEC)
        IF (ITYPE (IDEC+1) .NE. 111 .AND. ITYPE (IDEC+1) .NE. 221) THEN
          ITYPE (IDEC+2) = -ITYPE (IDEC+2)
          ITYPE (IDEC+3) = -ITYPE (IDEC+3)
        END IF
      END IF
c Decay to K0: 50% K0S, 50% K0L
      IF (ABS (ITYPE (IDEC+1)) .EQ. 311) THEN
        IF (H1RN () .LT. 0.50) THEN
          ITYPE (IDEC+1) = 130
        ELSE
          ITYPE (IDEC+1) = 310
        END IF
      END IF

      IF (ITYPE (IDEC+2) .EQ. 0) THEN
        N = 2
        if (stIDEBUG.ge.2) write(6,*) ' Decay to ', ITYPE (IDEC), ITYPE (IDEC+1)
      ELSE
        N = 3
        if (stIDEBUG.ge.2) write(6,*) ' Decay to ',
     +     ITYPE(IDEC), ITYPE(IDEC+1), ITYPE(IDEC+2)
      END IF

      ITER = 0
 9    CONTINUE
      ITER = ITER + 1
      SUM = 0D0
      DO 10, I = 1, N
        IPART = ITYPE (NPART+I)
        WIDTH = DBLE (PWIDTH (IPART))
        XM (I) = DBLE (PMASS (IPART))
        XMMIN = XM (I) - 2D0*WIDTH
        IF (ABS (IPART) .EQ. 1114 .OR. ABS (IPART) .EQ. 2114 .OR.
     +      ABS (IPART) .EQ. 2214 .OR. ABS (IPART) .EQ. 2224) THEN
          XMMIN = 1.1D0
        END IF
        IF (WIDTH .GT. 1D-6) THEN
          XM (I) = RANBW (XM (I), WIDTH, XMMIN, XM (I) + 2D0*WIDTH)
        END IF
        SUM = SUM + XM (I)
 10   CONTINUE

c      IF (SUM .GT. PPCMS8 (5, INSTAR) .AND. ITER .LT. 100) then
c       write(6,*) 'pb SUM  ITER = ',ITER
c        write(6,*) 'SUM, PPCMS8 (5, INSTAR) ',SUM, PPCMS8 (5, INSTAR)
c      endif

      IF (SUM .GT. PPCMS8 (5, INSTAR) .AND. ITER .LT. 100) GO TO 9

C     PRINT *, ' ITER: ', ITER

      IF (SUM .GT. PPCMS8 (5, INSTAR)) GO TO 7

*
* Call RAMBO: distribute particles uniformely in phase space
*

      CALL RAMBO (N, PPCMS8 (5, INSTAR), XM, P, WT, 1, IOK)
      IF (IOK .NE. 0) THEN
        CALL ERRLOG (202, 'S: DECNST: RAMBO failed!')
      END IF

*
* Update /CPROD/
*

      ISTAT (INSTAR) = 2
      IDAHEP (1, INSTAR) = NPART + 1
      IDAHEP (2, INSTAR) = NPART + N

      DO 50, I = 1, N
        PDPRS (4) = XM (I)**2
        DO 30, J = 1, 3
          PDPRS (J) = P (J, I)
          PDPRS (4) = PDPRS (4) + PDPRS (J)**2
 30     CONTINUE
        PDPRS (4) = SQRT (PDPRS (4))
        CALL LORENB8 (PPCMS8(5,INSTAR), PPCMS8(1,INSTAR), PDPRS, PDPCMS)
        DO 40, J = 1, 4
          PPCMS8 (J, NPART+I) = PDPCMS (J)
 40     CONTINUE
        PPCMS8 (5, NPART+I) = XM (I)
        ISTAT (NPART + I) = 1
        IDAHEP (1, NPART + I) = 0
        IDAHEP (2, NPART + I) = 0
        MOHEP (1, NPART + I) = INSTAR
        MOHEP (2, NPART + I) = 0
C       CALL PR4V8 ('PDPRS ', PDPRS)
C       CALL PR4V8 ('PDPCMS', PDPCMS)
C       CALL PR5V8 ('PPCMS8', PPCMS8 (1, NPART+I))
 50   CONTINUE

      NPART = NPART + N



       if (stIRFRA.eq.0) then
*
* Perform Delta, pi0, rho, eta decay
*

      IF (ABS (ITYPE (IDEC)) .EQ. 1114 .OR.
     +    ABS (ITYPE (IDEC)) .EQ. 2114 .OR.
     +    ABS (ITYPE (IDEC)) .EQ. 2214 .OR.
     +    ABS (ITYPE (IDEC)) .EQ. 2224) THEN
        CALL DECDEL (IDEC, IDECN)
        hasdecayed(idec) = 1
      END IF

      IF (ITYPE (IDEC+1) .EQ. 111) THEN
        CALL DECPI0 (IDEC+1)
        hasdecayed(idec+1)= 1
      ELSE IF (ITYPE (IDEC+1) .EQ. 113 .OR.
     +         ABS (ITYPE(IDEC+1)) .EQ. 213) THEN
        IHELI = INT (3.0*H1RN()) - 1
        CALL DECRHO (IDEC+1, IHELI, IDUMMY)
        hasdecayed(idec+1)= 1
      ELSE IF (ITYPE (IDEC+1) .EQ. 221) THEN
        CALL DECETA (IDEC+1)
        hasdecayed(idec+1)= 1
      END IF

      IF (ITYPE (IDEC+2) .EQ. 111) THEN
        CALL DECPI0 (IDEC+2)
        hasdecayed(idec+2)= 1
      END IF

      endif



      n0 = NLUJETS

      do i=2,NPART
       NLUJETS = NLUJETS +1

       PLU(NLUJETS,1) = sngl(PPCMS8(1,i))
       PLU(NLUJETS,2) = sngl(PPCMS8(2,i))
       PLU(NLUJETS,3) = sngl(PPCMS8(3,i))
       PLU(NLUJETS,4) = sngl(PPCMS8(4,i))
       PLU(NLUJETS,5) = sngl(PPCMS8(5,i))

       KLU(NLUJETS,2) = ITYPE(i)
       KLU(NLUJETS,1) = 1
       if (hasdecayed(i).eq.1) KLU(NLUJETS,1) = 11

c       KLU(NLUJETS,3) = iline
        if (MOHEP (1,i).gt.1) then
         KLU(NLUJETS,3) = n0 + MOHEP (1,i) -1
        else
         KLU(NLUJETS,3) = iline
        endif
      enddo

      if (stIDEBUG.ge.1) call LULIST(1)

      RETURN
      END
