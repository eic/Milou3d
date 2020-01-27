*CMZ :          30/05/2004  08.45.31  by  H1 Saclay
*CMZ :  1.06/03 04/04/96  12.06.39  by  Benno List
*CMZ :  1.06/02 04/04/96  12.00.41  by  Benno List
*CMZ :  1.06/00 16/06/95  10.53.07  by  Benno List
*CMZ :  1.01/00 01/09/94  15.34.08  by  Benno List
*CMZ :  1.00/02 25/07/94  18.50.37  by  Benno List
*CMZ :  0.23/01 04/02/94  17.02.53  by  Benno List
*CMZ :  0.23/00 02/02/94  13.14.59  by  Benno List
*CMZ :  0.22/00 25/01/94  13.30.59  by  Benno List
*CMZ :  0.21/00 19/01/94  11.31.41  by  Benno List
*-- Author :    Benno List   17/01/94
************************************************************************
      SUBROUTINE FRAGPJ


      include 'dvcs.common'

      INTEGER IBP, IDIFP, ISCN
************************************************************************
*
*     Purpose:    Take hadronic mass at p vertex and let it decay
*     --------    via JETSET
*
*     Input:      IBP:    Number of beam proton in /CPROD/
*     ------      IDIFP:  Number of diffractive state in /CPROD/
*                 from common block /CGDIF/:
*                 IBEAMP
*                 from common block /CPROD/:
*                 NPART, PPCMS8
*                 from common block /CMASS/:
*                 DMP, DMPI0
*
*     Output:     ISCN:   Number of scattered nucleon /CPROD/
*     -------     to common block /CPROD/:
*                 NPART, PPCMS8, ITYPE, MOHEP, IDAHEP, ISTAT
*
*     Called by:  GENEVT
*     ----------
*
*     Author:     Benno List, 17.1.94
*     -------
*
*     Changed by:
*     -----------
*
*     Calling:    From JETSET 7.3:
*     --------    LU1ENT, LUGIVE, LUJOIN, LULIST, LUEXEC, PLU, KLU
*
************************************************************************


      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON /LUJETS/ NLUJETS,KLU(4000,5),PLU(4000,5),VLU(4000,5)
      COMMON /LUDAT1/ MSTU(200), PARU(200), MSTJ(200), PARJ(200)



      REAL*8 DPI
      PARAMETER (DPI = 3.14159265359D0)

      PARAMETER (PI = 3.14159265359E0)
      PARAMETER (RADDEG = 180E0/PI)

      logical SPLIT_POM


      include 'cprod.common'




      REAL*8 PQQ (5), PQ (5), PPOM (5)
      CHARACTER*100 C100
      INTEGER IARR12 (2) /1, 2/
      DATA NERR /0/

C     PRINT *, '--- FRAGPJ ---'
C     PRINT *, 'IBP   = ', IBP
C     PRINT *, 'IDIFP = ', IDIFP
C     PRINT *, 'NPART = ', NPART


      DMP = 0.93827d0
      DMPI0 = DBLE (PMAS (111,1))



c --- Fill the CPROD common :   -----------------------------------------

      NPART = 2

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
c ---    for the while, in the lab frame
      PPCMS8(1,1) = PPART8(1,1)
      PPCMS8(2,1) = PPART8(2,1)
      PPCMS8(3,1) = PPART8(3,1)
      PPCMS8(4,1) = PPART8(4,1)
      PPCMS8(5,1) = PPART8(5,1)
      ITYPE(1) = KLU(iline,2)
      MOHEP(1,1) = 0
      MOHEP(2,1) = 0
      IDAHEP(1,1) = 0
      IDAHEP(2,1) = 0
      ISTAT(1) = 1

      jline = 1
      PPART(1,2) = PLU(jline,1)
      PPART(2,2) = PLU(jline,2)
      PPART(3,2) = PLU(jline,3)
      PPART(4,2) = PLU(jline,4)
      PPART(5,2) = PLU(jline,5)
      do i=1,5
       PPART8(i,2) = dble(PPART(i,2))
      enddo
c ---    for the while, in the lab frame
      PPCMS8(1,2) = PPART8(1,2)
      PPCMS8(2,2) = PPART8(2,2)
      PPCMS8(3,2) = PPART8(3,2)
      PPCMS8(4,2) = PPART8(4,2)
      PPCMS8(5,2) = PPART8(5,2)

      ITYPE(2) = KLU(jline,2)
      MOHEP(1,2) = 0
      MOHEP(2,2) = 0
      IDAHEP(1,2) = 0
      IDAHEP(2,2) = 0
      ISTAT(2) = 1


       do i=1,NPARTMAX
        hasdecayed(i) = 0
       enddo


       IDIFP = 1
       IBP   = 2
       IBEAMP = 2212


c -----------------------------------------------------------------------

      IF (PPCMS8 (5, IDIFP) .LT. DMP+DMPI0) THEN
        PRINT *, '### FRAGPJ error: not enough energy! '
        PRINT *, 'IBP:   ', IBP
        PRINT *, 'IDIFP: ', IDIFP
        PRINT *, 'PPCMS8 (5, IDIFP): ', PPCMS8 (5, IDIFP)
        PRINT *, 'DMP + DMPI0:       ', DMP+DMPI0
        CALL ERRLOG (40, 'F: FRAGPJ: Not enough energy!')
        CALL H1STOP
      END IF


      if (stPROSPLIT.eq.1) then
       SPLIT_POM = .true.
      else if (stPROSPLIT.eq.0) then
       SPLIT_POM = .false.

      else
       proba = H1RN()
       if (proba.le.stPROSPLIT) then
        SPLIT_POM = .true.
       else
        SPLIT_POM = .false.
       endif

      endif

      if (SPLIT_POM) then
*
* Split proton into q - qq system: choose q and qq codes
* The pomeron couples with the same strength to u and d quarks.
* ud_1 and uu_1 are spin triplets, ud_0 is a spin singulett.
* Therefore q-qq is 1/7 u-ud_0, 3/7 u-ud_1, and 3/7 d-uu_1.
*
      R = 7.0*H1RN()
      IF (R .LE. 1.0) THEN
        IQ  = SIGN (2, IBEAMP)
        IQQ = SIGN (2101, IBEAMP)
      ELSE IF (R .LE. 4.0) THEN
        IQ  = SIGN (2, IBEAMP)
        IQQ = SIGN (2103, IBEAMP)
      ELSE
        IQ  = SIGN (1, IBEAMP)
        IQQ = SIGN (2203, IBEAMP)
      END IF

      else
*
* q-qq probabilities in the ratio of the quark charge**4
*
      R = 67.0*H1RN()
      IF (R .LE. 3.0) THEN         ! d (uu)_1
        IQ  = SIGN (1, IBEAMP)
        IQQ = SIGN (2203, IBEAMP)
      ELSE IF (R .LE. 19.0) THEN   ! u (ud)_0
        IQ  = SIGN (2, IBEAMP)
        IQQ = SIGN (2101, IBEAMP)
      ELSE                         ! u (ud)_1
        IQ  = SIGN (2, IBEAMP)
        IQQ = SIGN (2103, IBEAMP)
      END IF


      endif


C     PRINT *, 'IQ:  ', IQ,  ', mass: ', ULMASS (IQ)
C     PRINT *, 'IQQ: ', IQQ, ', mass: ', ULMASS (IQQ)

* Calculate the 4-vector of the proton remnant and of the struck quark
* To avoid numerical problems, calculate pomeron momentum from
* proton side
      DO 10, I = 1, 4
        PPOM (I) = PPCMS8 (I, IDIFP) - PPCMS8 (I, IBP)
 10   CONTINUE
      PPOM (5) = -SQRT (PPOM (1)**2+PPOM (2)**2+PPOM (3)**2-PPOM(4)**2)
      CALL SPLITP (PPCMS8 (1, IBP), PPOM, PQ, PQQ, IERR)

      NTRIAL = 0
 12   CONTINUE
        NTRIAL = NTRIAL+1


c --- diquark :

      NLUJETS = NLUJETS +1
      KLU(NLUJETS,2) = IQQ

      PLU(NLUJETS,1) = sngl(PQQ(1))
      PLU(NLUJETS,2) = sngl(PQQ(2))
      PLU(NLUJETS,3) = sngl(PQQ(3))
      PLU(NLUJETS,4) = sngl(PQQ(4))
      PLU(NLUJETS,5) = sngl(PQQ(5))

      KLU(NLUJETS,3) = iline
      KLU(NLUJETS,1) = 3

       KLU(NLUJETS,5) = NLUJETS+1
       KLU(NLUJETS,4) = 10000 * (NLUJETS+1)

c --- quark :

      NLUJETS = NLUJETS +1
      KLU(NLUJETS,2) = IQ

      PLU(NLUJETS,1) = sngl(PQ(1))
      PLU(NLUJETS,2) = sngl(PQ(2))
      PLU(NLUJETS,3) = sngl(PQ(3))
      PLU(NLUJETS,4) = sngl(PQ(4))
      PLU(NLUJETS,5) = sngl(PQ(5))

      KLU(NLUJETS,3) = iline
      KLU(NLUJETS,1) = 3

       KLU(NLUJETS,4) =  NLUJETS-1
       KLU(NLUJETS,5) =  10000 * (NLUJETS-1)

c      call lulist(1)


      goto 999

* Define string between qq and q
        IARR12(1) = NLUJETS-1
        IARR12(2) = NLUJETS
        CALL LUJOIN (2, IARR12)

C       CALL LULIST (1)

* Perform fragmentation with JETSET
        CALL LUEXEC

* Check for possible Error
        IF (MSTU (24) .NE. 0) THEN
          CALL ERRLOG (41, 'I: FRAGPJ: JETSET error')
          NERR = NERR + 1
          IF (NERR .LT. 100) THEN
            PRINT *, ' *** FRAGPJ: JETSET error', MSTU (24)
          END IF
        ELSE IF (MSTU (28) .NE. 0) THEN
          CALL ERRLOG (42, 'I: FRAGPJ: JETSET warning')
          NERR = NERR + 1
          IF (NERR .LT. 100) THEN
            PRINT *, ' *** FRAGPJ: JETSET warning', MSTU (28)
          END IF
        END IF

      IF (NTRIAL .LE. 10 .AND.
     +    (MSTU (24) .NE. 0 .OR. MSTU (28) .NE. 0)) GO TO 12

      IF (MSTU (24) .NE. 0) THEN
        PRINT *, ' *** FRAGPJ: Cannot prevent JETSET error'
        PRINT *, '     Mass: ', PPCMS8 (5, IDIFP)
        CALL ERRLOG (43, 'S: FRAGPJ: Cannot prevent JETSET error')
      ELSE IF (MSTU (28) .NE. 0) THEN
        PRINT *, ' *** FRAGPJ: Cannot prevent JETSET warning'
        PRINT *, '     Mass: ', PPCMS8 (5, IDIFP)
        CALL ERRLOG (44, 'W: FRAGPJ: Cannot prevent JETSET warning')
      END IF

C     CALL LULIST (1)

* Print control output:
C     PRINT *, 'FRAGPJ: Result'
C     PRINT *, 'ISCN: ', ISCN
C     CALL PR5V8 ('p_pdif', PPCMS8 (1, IDIFP))
C     CALL PR5S  ('JETSET',PLU(0,1),PLU(0,2),PLU(0,3),PLU(0,4),PLU(0,5))


999   continue
      RETURN
      END
