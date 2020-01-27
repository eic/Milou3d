************************************************
*          Printing Event Information          *
*        of each event into ASCII file         *
*             S. Fazio, April 2011             *
************************************************
* Input : /LUJETS/ (from pythia comkons), event #(IEvt), LUN for ASCII file
* Output: ASCII file
*
*
      subroutine PRINT_asc(LUN,ievt)

c      implicit NONE
      include 'dvcs.common'
      include 'forpaw.common'

* ------------- PYTHIA common -------------
      integer      N,  K(4000,5)
      real            P(4000,5), V(4000,5)
c      double precision       P(4000,5), V(4000,5)
      COMMON /LUJETS/ N,K,P,V
cvcv      include 'pythia.inc'              ! All PYTHIA commons blocks

C * ------ Argument ------
      integer  LUN,  IEvt

      common /dvcs_VAR/x_main,q_main,phi_main,t_main,ym_main


* ---------- Event record stuff ----------
      integer           K_e(5),K_p(5)
      double precision  P_e(5),P_p(5)

      integer  K6, i, j, IN
      integer  USGEvt_IPA(1:10)
      real*4   Weight, USGEvt_RPA(1:10)

* ---------------- DATA ----------------
      data  K6/0/
* -------------------------------------------
      
        if (Ievt.GE.1.and.Ievt.LE.stNgen) then !!!!!!!!!!!!
         weight=1

* <========= Writing information on each event =========>
         N=5
         if(stirad.eq.1.and.p(4,4).gt.0.0001) N=6
* ------ General information ------
C   This is what we write in the ascii-file

c        write(LUN,*)' PYTHIA EVENT FILE '
c        write(LUN,*)'============================================'
c        write(LUN,30) 

cc      write(LUN,*) IEvt , N , Weight !event num, lines num, weight

      USGEvt_RPA(1)=xntp
      USGEvt_RPA(2)=Qntp
      USGEvt_RPA(3)=yntp
      USGEvt_RPA(4)=-tntp       ! -t  
      USGEvt_RPA(5)=phintp
cc      USGEvt_RPA(5)=sqrt(Qntp*(1/xntp-1)+0.93827231**2)  ! W
c      USGEvt_RPA(6)=-tntp                                ! -t
      USGEvt_RPA(6)=phibelgen                              !phi bel gen
      USGEvt_RPA(7)=resol                              !phi resolution
      USGEvt_RPA(8)=phibelrec                              !phi bel reconstucted
c      USGEvt_RPA(8)=pphi   !pphi = angle between lepton plane and transverse polarization vec (see Belitsky, Mueller and Kirchner)

      USGEVT_IPA(1)=stipro     ! process type generated 
      USGEVT_IPA(2)=stIRAD
      USGEVT_IPA(3)=0
      USGEVT_IPA(4)=0
      USGEVT_IPA(5)=0
      USGEVT_IPA(6)=0
      USGEVT_IPA(7)=0
      USGEVT_IPA(8)=0
      USGEVT_IPA(9)=0
      USGEVT_IPA(10)=0

c write: event num, lines num, weight
       write(LUN,*) 0, IEvt , N , Weight,
c       write(LUN,*)
c     +        (USGEvt_IPA(j),j=1,10) , (USGEvt_RPA(j),j=1,10)
     +        (USGEvt_IPA(j),j=1,2) , (USGEvt_RPA(j),j=1,8)
       write(LUN,*)'============================================'

* -------- SWAP for ZDIS --------
         do i=1,5
           P_e(i) = P(2,i)
           K_e(i) = K(2,i)
           P_p(i) = P(1,i)
           K_p(i) = K(1,i)
         enddo
         do i=1,5
           P(1,i) = P_e(i)
           K(1,i) = K_e(i)
           P(2,i) = P_p(i)
           K(2,i) = K_p(i)
         enddo
         do 500 IN = 3, N
           if (K(IN,3).EQ.1) then
             K(IN,3) = 2
             GOTO 500
           endif
           if (K(IN,3).EQ.2)  K(IN,3) = 1
 500     continue
* -------------------------------

* --------- Event record ---------
         do IN = 1, N
c           write(LUN,*) (real(P(IN,i)),i=1,5), (K(IN,j),j=1,5),K6
         enddo
* --------------------------------
      write(lun,*) 1, 21,  k(1,2), 0, 0, 0,         ! e
     >  p(1,1),p(1,2),p(1,3),p(1,4),p(1,5),
     >  0,0,0
      write(lun,*) 2, 21, k(2,2), 0, 0, 0,          ! p
     >  p(2,1),p(2,2),p(2,3),p(2,4),p(2,5), 
     >  0,0,0

      if (stIRAD.eq.0) then
      write(lun,*) 3, 1,  k(3,2), 1, 2, 0,              ! e'
     >   p(3,1),p(3,2),p(3,3),p(3,4),p(3,5),
     >   0,0,0
      write(lun,*) 4, k(5,1),  k(5,2), 1, 2, 0,              ! virtual gamma
     >   p(5,1),p(5,2),p(5,3),p(5,4),p(5,5),
     >   0,0,0 
      write(lun,*) 5, 1,  k(6,2), 1, 2, 0,              ! gamma
     >     p(6,1),p(6,2),p(6,3),p(6,4),p(6,5),
     >     0,0,0
      write(lun,*) 6, 1, k(4,2), 1, 2, 0,               ! p'
     >     p(4,1),p(4,2),p(4,3),p(4,4),p(4,5), 
     >     0,0,0

      elseif(stIRAD.eq.1) then
      write(lun,*) 3, 1,  k(6,2), 1, 2, 0,               ! e'
     >    p(6,1),p(6,2),p(6,3),p(6,4),p(6,5), 
     >    0,0,0
      write(lun,*) 4, k(4,1),  k(4,2), 1, 2, 0,              ! virtual gamma
     >    p(4,1),p(4,2),p(4,3),p(4,4),p(4,5), 
     >    0,0,0 
      write(lun,*) 5, 1,  k(8,2), 1, 2, 0,              ! gamma
     >    p(8,1),p(8,2),p(8,3),p(8,4),p(8,5), 
     >    0,0,0
      write(lun,*) 6, 1, k(7,2), 1, 2, 0,               ! p'
     >     p(7,1),p(7,2),p(7,3),p(7,4),p(7,5),
     >     0,0,0
      if (p(5,4).gt.0.0001) then  ! only store ISR gamma if energetic enough 
      write(lun,*) 7, 1,   k(5,2), 1, 2, 0,             ! ISR gamma
     >     p(5,1),p(5,2),p(5,3),p(5,4),p(5,5), 
     >     0,0,0
      endif
      endif

        endif !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

** 34      format(2(I6,1x,$),I10,1x,$,3(I8,1x,$),8(f15.6,1x,$),/)
         write(LUN,*)'=============== Event finished ==============='

********* Termination of ASCII Event Storing *********
c      if (Ievt .GE. 1) then
      if (Ievt .GE. stNgen) then
        close(LUN)
      endif
******************************************************

      return
      end
*#############################################
**********************************************
*         Printing header of ASCII file      *
**********************************************
* Input : LUN for ASCII file
* Output: ASCII file
*

      subroutine  Print_asc_head(LUN)
      implicit NONE

      include 'dvcs.common'

* ------ Argument ------
      integer  LUN
* ----------------------

* ------------ Event record stuff for TOZDIS ------------
      character
     +   USGRun_Comment(1:2)*20 , DISGenerator*8 , Institute*8
      integer
     +   EvtCTyp , KRows , LRows ,
     +   K6
      integer   USGRun_IPA(1:10)
      real*4    USGRun_RPA(1:10)
* -------------------------------------------------------
* -------- Local variables --------
      integer  IKF,NKF,KF
       parameter  (NKF=17)
      integer       KFSD(NKF)
      character*16  NASD(NKF)
      integer   i
* --------------------------------- 
* ---------------- DATA ----------------
      data KFSD / 310,  221, 3122, 3222, 3212
     &           ,3112, 3322, 3312, 3334
     &           ,411,   421,  431, 4122
     &           ,13,        211, 321, 130 /

      data  K6/0/
* --------------------------------------
 
 
********************** <<<<<< TOZDIS >>>>>> ********************
* ------ General informations on this generator ------
      USGRun_Comment(1) = 'MILOU32'  !Name and Version of the generator
ccc      USGRun_Comment(2) = 'GPDmodel'      ! model used
ccc      USGRun_Comment(3) = 'DVCS'    ! type of prosess
ccc      USGRun_Comment(4) = 'noF2set'      ! structure function set, if used
c      USGRun_Comment(4) = 'ALLM'      ! structure function set, if used
ccc      USGRun_Comment(5) = 'S. Fazio'  !contact person who produced this file
      USGRun_Comment(2) = 'S. Fazio'  !contact person who produced this file
cm      DISGenerator = 'MILOU'     ! name as a DIS generator
      Institute = 'BNL'          ! institute of the contact person
cm      EvtCTyp = 0                  ! event common type
*                                        1 for LUND common
*                                        2 for HEP common
      KRows = 1       ! # of rows to be filled in the USGRun table
      LRows = 1       ! # of rows to be filled in the USGEvt table
 
* <-------- Run parameters and variables -------->
 
      USGRun_IPA( 1) = stipro   ! process type generated (1=el 2=inel...)
      USGRun_IPA( 2) = stirad  ! does dataset include isr or not?
      USGRun_IPA( 3) = stielas ! IELAS 1=elastic 0=dissociation
      USGRun_IPA( 4) = 0
      USGRun_IPA( 5) = 0
      USGRun_IPA( 6) = 0
      USGRun_IPA( 7) = 0
      USGRun_IPA( 8) = 0
      USGRun_IPA( 9) = 0
      USGRun_IPA(10) = 0
 
      USGRun_RPA( 1) = stELEP      ! electron beam momentum before isr 
      USGRun_RPA( 2) = stetarg      ! proton beam momentum 
      USGRun_RPA( 3) = 0.
      USGRun_RPA( 4) = 0.
      USGRun_RPA( 5) = 0.
      USGRun_RPA( 6) = 0.
      USGRun_RPA( 7) = 0.
      USGRun_RPA( 8) = 0.
      USGRun_RPA( 9) = 0.
      USGRun_RPA(10) = 0.
 
C ------ Opening an output ascii file ------
      open(LUN, FILE='./asc_15x50_dvcs.out', STATUS='unknown', err=80)
           GOTO 90
 80        write(6,*) '!!!Error in PRINT_asc_head!!!'
           write(6,*) '  ---> cannot open ASCII file!'
           write(6,*) '  ---> Good By!'
           STOP
 90   continue
C   This is what we write in the ascii-file
C ------ Writing parameters and variables of the RUN ------
      write(LUN,'(6A12/A8,2I3)') 
     +     USGRun_Comment(1), 
     +     USGRun_Comment(2), 
ccc     +     USGRun_Comment(3), USGRun_Comment(4), 
ccc     +     USGRun_Comment(5),
     +     Institute , KRows , LRows
 
ccc      write(LUN,*)
ccc     +     (USGRun_IPA(i),i=1,3) , (USGRun_RPA(i),i=1,2)

        write(LUN,*)' MILOU EVENT FILE '
        write(LUN,*)'============================================'
        write(LUN,30) 
30      format('I, ievent, linesnum, weight, genprocess, radcorr, 
     & truex, trueQ2, truey, truet, treuphi, phibelgen, phibelres,
     & phibelrec') 
cc        write(LUN,31)
cc31      format('genprocess, radcorr, truex, trueQ2, truey, truet, 
cc     & treuphi')
        write(LUN,*)'============================================'

        write(LUN,*)'I, K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5) 
     &   P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5) V(I,1)  V(I,2)  V(I,3)'
        write(LUN,*)'============================================'
 
 
      return
      end
*#########################################################################
