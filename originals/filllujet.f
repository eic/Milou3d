*CMZ :          30/05/2004  07.46.27  by  H1 Saclay
*-- Author :    Unknown   17/12/2003


*=========================================
      SUBROUTINE FILLLUJET(ifail)
*=========================================

      include 'dvcs.common'
      include 'forpaw.common'

      PARAMETER(NMAX = 4000)

      COMMON /LUJETS/ N,K(4000,5),P(4000,5),V(4000,5)
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON /LUDAT3/ MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)

      real elepi,EGAMR,sweight,srad,elepin,ehadi
      common /RADGEN/ elepi,EGAMR,sweight,srad,elepin,ehadi


* --- Common for steering, read via FFKEY

      REAL*4 SPACE
      COMMON /CFREAD/SPACE(1000)



* --- Common /dvcsvar1/ Filled in calc_kinem --------------------
*     -> in SEQ DVCSPaw.
* ------------------------------------------------------------


      if (stIDEBUG.ge.2) then
       write(6,*) 'debut de FILLLUJET'
      endif

C -> Initialize :
      N = 0
      MSTU(3) = 0
      CALL VZERO(K,NMAX)
      CALL VZERO(P,NMAX)
      CALL VZERO(V,NMAX)

C*
C*    save beam particles in the first two lines
C*    ==========================================================

      if (stIRAD.eq.0) then

C -> Initialize :
      N = 0
      MSTU(3) = 0
      CALL VZERO(K,NMAX)
      CALL VZERO(P,NMAX)
      CALL VZERO(V,NMAX)

C -> Incoming Proton :

      N = N+1
      K(N,1)= 21
      K(N,2)= 2212    ! a adapter pour cibles qcq
      K(N,3)=  0
      K(N,4)=  0
      K(N,5)=  0
      P(N,1)=  0.0
      P(N,2)=  0.0
      if (.not.stFIXED) then
       P(N,3)=  stETARG
       P(N,5)= ULMASS(ABS(K(N,2)))
       P(N,4)=  sqrt(P(N,3)**2 + P(N,5)**2)
      else
       P(N,3)= 0.
       P(N,5)= ULMASS(ABS(K(N,2)))
       P(N,4)= P(N,5)
      endif

      V(N,1)=  0.0
      V(N,2)=  0.0
      V(N,3)=  0.0
      V(N,4)=  0.0
      V(N,5)=  0.0

C -> Incoming Electron :

      N = N+1
      K(N,1)= 21
      K(N,2) = 11    !  a adapter pour muon
      if (stLCHAR > 0) K(N,2) = -K(N,2)
      K(N,3)=  0
      K(N,4)=  0
      K(N,5)=  0
      P(N,1)=  0.0
      P(N,2)=  0.0
      P(N,3)=  stELEP
      P(N,5)=  ULMASS(IABS(K(N,2)))
      P(N,4)=  sqrt(P(N,3)**2 + P(N,5)**2)
      V(N,1)=  0.0
      V(N,2)=  0.0
      V(N,3)=  0.0
      V(N,4)=  0.0
      V(N,5)=  0.0


C --- Final particles

C --------  Outcoming lepton

      N = N+1
      K(N,1)= 21
      K(N,2) = K(2,2)
      V(N,1)=  0.0
      V(N,2)=  0.0
      V(N,3)=  0.0
      V(N,4)=  0.0
      V(N,5)=  0.0
      K(N,3)=  0
      K(N,4)=  0
      K(N,5)=  0
      P(N,1) = plolab(1)
      P(N,2) = plolab(2)
      P(N,3) = plolab(3)
      P(N,4) = plolab(4)
c      write(6,*) P(N,4)**2-P(N,1)**2-P(N,2)**2-P(N,3)**2
      P(N,5)=  ULMASS(IABS(K(N,2)))
c      P(N,5) = sqrt(P(N,4)**2-P(N,1)**2-P(N,2)**2-P(N,3)**2)


C --------  Outcoming hadron

      N = N+1
      K(N,1)= 21

      if (stIELAS.eq.1) then
       K(N,2) = K(1,2)
      else
       call PDGCODE_RESO(ipout)
       if (iabs(ipout).ne.2210) then
        K(N,2) = KC_FREE(ipout)
       else
        K(N,2) = ipout
       endif
      endif

      V(N,1)=  0.0
      V(N,2)=  0.0
      V(N,3)=  0.0
      V(N,4)=  0.0
      V(N,5)=  0.0
      K(N,3)=  0
      K(N,4)=  0
      K(N,5)=  0
      P(N,1) = ppolab(1)
      P(N,2) = ppolab(2)
      P(N,3) = ppolab(3)
      P(N,4) = ppolab(4)
      P(N,5) = sqrt(P(N,4)**2-P(N,1)**2-P(N,2)**2-P(N,3)**2)
cep       P(N,5) = ym


C --------  Outcoming real gamma

      N = N+1
      K(N,1)= 21
      K(N,2) = 22
      V(N,1)=  0.0
      V(N,2)=  0.0
      V(N,3)=  0.0
      V(N,4)=  0.0
      V(N,5)=  0.0
      K(N,3)=  0
      K(N,4)=  0
      K(N,5)=  0
      P(N,1) = prglab(1)
      P(N,2) = prglab(2)
      P(N,3) = prglab(3)
      P(N,4) = prglab(4)
c      P(N,5) = sqrt(P(N,4)**2-P(N,1)**2-P(N,2)**2-P(N,3)**2)
      P(N,5) = 0.


C --- Copy outcoming particles

      do i=6,8
       n = n+1
       j = i-3
       do kk=1,5
        k(i,kk) = k(j,kk)
       enddo
       do kk=1,5
        v(i,kk) = v(j,kk)
        p(i,kk) = p(j,kk)
       enddo
      enddo

      k(6,1) =  1    ! lepton can radiate
      if (stIELAS.eq.1) then
       k(7,1) =  1
      else
       k(7,1) =  11
      endif

      k(8,1) =  1

      endif ! IRAD=0

C-----------------------------------------------------

      if (stIRAD.eq.1) then

C -> Initialize :
      N = 0
      MSTU(3) = 0
      CALL VZERO(K,NMAX)
      CALL VZERO(P,NMAX)
      CALL VZERO(V,NMAX)

C -> Incoming Proton :

      N = N+1     ! n=1
      K(N,1)= 21
      K(N,2)= 2212
      K(N,3)=  0
      K(N,4)=  0
      K(N,5)=  0
      P(N,1)=  0.0
      P(N,2)=  0.0
      if (.not.stFIXED) then
       P(N,3)=  stETARG
       P(N,5)= ULMASS(ABS(K(N,2)))
       P(N,4)=  sqrt(P(N,3)**2 + P(N,5)**2)
      else
       P(N,3)= 0.
       P(N,5)= ULMASS(ABS(K(N,2)))
       P(N,4)= P(N,5)
      endif

      V(N,1)=  0.0
      V(N,2)=  0.0
      V(N,3)=  0.0
      V(N,4)=  0.0
      V(N,5)=  0.0

C -> Incoming Electron :

      N = N+1     ! n=2
      K(N,1)= 21
      K(N,2) = 11
      if (stLCHAR > 0) K(N,2) = -K(N,2)
      K(N,3)=  0
      K(N,4)=  0
      K(N,5)=  0
      P(N,1)=  0.0
      P(N,2)=  0.0
      P(N,3)=  stELEP
      P(N,5)=  ULMASS(IABS(K(N,2)))
      P(N,4)=  sqrt(P(N,3)**2 + P(N,5)**2)
      V(N,1)=  0.0
      V(N,2)=  0.0
      V(N,3)=  0.0
      V(N,4)=  0.0
      V(N,5)=  0.0

C -> Interacting lepton

      N = N+1    ! n=3 e- incident apres radiation
      K(N,1)= 21
      K(N,2) = 11
      if (stLCHAR > 0) K(N,2) = -K(N,2)
      K(N,3)=  0
      K(N,4)=  0
      K(N,5)=  0
      P(N,1)=  0.0
      P(N,2)=  0.0
      P(N,3)=  -elepi
      P(N,5)=  ULMASS(IABS(K(N,2)))
      P(N,4)=  sqrt(P(N,3)**2 + P(N,5)**2)
      V(N,1)=  0.0
      V(N,2)=  0.0
      V(N,3)=  0.0
      V(N,4)=  0.0
      V(N,5)=  0.0

c -> ISR gamma
      N = N+1    ! n=4
      K(N,1)= 21
      K(N,2) = 22
      K(N,3)=  N-1
      K(N,4)=  0
      K(N,5)=  0
      P(N,1)=  0.0
      P(N,2)=  0.0
      P(N,3)=  -EGAMR
      P(N,5)=  0.
      P(N,4)=   EGAMR
      V(N,1)=  0.0
      V(N,2)=  0.0
      V(N,3)=  0.0
      V(N,4)=  0.0
      V(N,5)=  0.0

C --- Final particles

C --------  Outcoming lepton

      N = N+1    ! n=5
      K(N,1)= 21
      K(N,2) = K(2,2)
      V(N,1)=  0.0
      V(N,2)=  0.0
      V(N,3)=  0.0
      V(N,4)=  0.0
      V(N,5)=  0.0
      K(N,3)=  0
      K(N,4)=  0
      K(N,5)=  0
      P(N,1) = plolab(1)
      P(N,2) = plolab(2)
      P(N,3) = plolab(3)
      P(N,4) = plolab(4)
c      write(6,*) P(N,4)**2-P(N,1)**2-P(N,2)**2-P(N,3)**2
      P(N,5)=  ULMASS(IABS(K(N,2)))
c      P(N,5) = sqrt(P(N,4)**2-P(N,1)**2-P(N,2)**2-P(N,3)**2)

C --------  Outcoming hadron

      N = N+1    ! n=6

      K(N,1)= 21

      if (stIELAS.eq.1) then
       K(N,2) = K(1,2)
      else
       call PDGCODE_RESO(ipout)
       if (iabs(ipout).ne.2210) then
        K(N,2) = KC_FREE(ipout)
       else
        K(N,2) = ipout
       endif
      endif

      V(N,1)=  0.0
      V(N,2)=  0.0
      V(N,3)=  0.0
      V(N,4)=  0.0
      V(N,5)=  0.0
      K(N,3)=  0
      K(N,4)=  0
      K(N,5)=  0
      P(N,1) = ppolab(1)
      P(N,2) = ppolab(2)
      P(N,3) = ppolab(3)
      P(N,4) = ppolab(4)
      P(N,5) = sqrt(P(N,4)**2-P(N,1)**2-P(N,2)**2-P(N,3)**2)
c      P(N,5) = ym


C --------  Outcoming real gamma

      N = N+1    ! n=7
      K(N,1)= 21
      K(N,2) = 22
      V(N,1)=  0.0
      V(N,2)=  0.0
      V(N,3)=  0.0
      V(N,4)=  0.0
      V(N,5)=  0.0
      K(N,3)=  0
      K(N,4)=  0
      K(N,5)=  0
      P(N,1) = prglab(1)
      P(N,2) = prglab(2)
      P(N,3) = prglab(3)
      P(N,4) = prglab(4)
c      P(N,5) = sqrt(P(N,4)**2-P(N,1)**2-P(N,2)**2-P(N,3)**2)
      P(N,5) = 0.

      do i=8, 11  ! 8-9-10-11
       n = n+1
       j = i-3
       if (i.eq.11) j=4
       do kk=1,5
        k(i,kk) = k(j,kk)
       enddo
       do kk=1,5
        v(i,kk) = v(j,kk)
        p(i,kk) = p(j,kk)
       enddo
      enddo


      k(8,1)  =  1
      if (stIELAS.eq.1) then
       k(9,1)  =  1
      else
       k(9,1)  =  11
      endif

      k(10,1) =  1
      k(11,1) =  1


      endif ! IRAD=1


      if (stIDEBUG.ge.1)then
       write(6,*) 'FILLLUJET : before LUPREP '
       call LULIST(1)
      endif


      if (stIELAS.eq.0) then

       if (ym.le.1.9) then
        call DECNST(1,IDEC, IDECN)

       else
        call FRAGPJ

       endif
      endif


      call luprep(21)

cep      call lushow(6,8,dble(sqrt(p(6,1)**2+p(6,2)**2) ))

      CALL LUEXEC


      if (stIDEBUG.ge.1) then
       write(6,*) 'fin de FILLLUJET'
       call LULIST(1)
      endif

      return
      end




