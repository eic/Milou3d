*CMZ :          05/10/2004  13.26.59  by  H1 Saclay
*-- Author :    Unknown   17/12/2003


******************************************************************
*
      subroutine dvcs_ini(mfail)
*
******************************************************************

      implicit double precision (a-h,o-z)

      include 'dvcs.common'

      external fdvcs
      external dmy

* --- Common for steering, read via FFKEY

      REAL*4 SPACE
      COMMON /CFREAD/SPACE(1000)

      real*4 rndm

      logical Angle_Integ
      common/dvcs_IN/spin,Z,A,s,tchar,xlaml,xlamp,
     >               Angle_Integ,iord,icount,nx,nq,IPRO
      common/dvcs_OUT/res
      common/dvcs_VAR/x_main,q_main,phi_main,t_main,ym_main

      real elepi,EGAMR,sweight,srad,elepin,ehadi
      common /RADGEN/ elepi,EGAMR,sweight,srad,elepin,ehadi

      logical stFIXED_save
      common /FIXED/ stFIXED_save

      PARAMETER (MXDIM = 50 )
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,IG(MXDIM),NCALL
      COMMON /BPARM2/ ACC1,ACC2,ITMX1,ITMX2

      REAL*4 SIGMATOT,DSIGMATOT
      COMMON/XSECTOT/SIGMATOT,DSIGMATOT

      COMMON/STEP/ISTEP

      integer debug
      common/debug/idebug

      integer i_outrange
      common/out_of_range/i_outrange

      real*8 mp 
      real*8  it_form
      real*8  bnew_q,bnew_g,xdel_s,bnew_slopeq,bnew_slopeg,q20
      common/form_factors/it_form,bnew_q,bnew_g,xdel_s,q20,
     +                            bnew_slopeq,bnew_slopeg

      integer isel_twist3
      common /isel_twist3/ isel_twist3

      integer isel_asym
      common /isel_asym/ isel_asym

      integer isel_interf
      common /isel_interf/ isel_interf

      integer isel_realpart
      common /isel_realpart/ isel_realpart

      integer isel_tintin
      common /isel_tintin/ isel_tintin,isel_f2qcd,isel_dipole

      real*8 bti,rti
      common /tintin/ SIGN,bti,rti

      integer igen_exp
      common /genexp/ igen_exp

* => COMMON FOR PAW :
      REAL hmemor
      PARAMETER( NWPAWC=500000 )
      COMMON/pawc/hmemor(NWPAWC)


      common /order/ iiord

      common /counter1/ icountxq,icountt

      common /setting/ nset

      parameter(mx=300)
      dimension f(mx),f1(mx),res(100),f2(mx),f3(mx)
      dimension f4(mx),f5(mx),f6(mx),f7(mx),f8(mx),f9(mx),resut(2,mx)
      dimension f10(mx),f11(mx)
      character*2 temp,temp0
      character*78 temp1,temp2

      mfail = 0


      call readsteer

C --------------------------------------------------------
      if (stIPRO.eq.1) isel_asym   = 0 !stASYM
      if (stIPRO.eq.2) isel_asym   = 0 !stASYM
      if (stIPRO.eq.3) isel_asym   = 0 !stASYM
      if (stIPRO.eq.4) isel_asym   = 1 !stASYM
      if (stIPRO.eq.5) isel_asym   = 1 !stASYM

Cls---------------------------------------new
      isel_realpart = 1 ! 0:from Freund param ; 1:from lambda ; 2:Re=0
Cls---------------------------------------new

      isel_interf = 0
      igen_exp    = stEXP
ccc      print *,'isel_asym,igen_exp=',isel_asym,igen_exp
C --------------------------------------------------------

C --------------------------------------------------------
      if (stFIXED) stETARG = 0.93827

      elepin  = abs(stELEP)
      elepi   = abs(stELEP)
      ehadi   = abs(stETARG)
      EGAMR   = 0.
      sweight = 1.
      srad    = real(stIRAD)
C --------------------------------------------------------

C --------------------------------------------------------
      if (stTINTIN) then
      isel_tintin = 1
      else
      isel_tintin = 0
      endif
      SIGN = -dble(stLCHAR) !dble(stSIGN)
      bti  =  dble(stBTIN)
      rti  =  dble(stRTIN)
      if (stF2QCD) then
      isel_f2qcd = 1
      else
      isel_f2qcd = 0
      endif
      if (stDIPOLE) then
      isel_dipole = 1
      else
      isel_dipole = 0
      endif
C --------------------------------------------------------

C --------------------------------------------------------
      it_form     = dble(real(stITFORM))
      bnew_q      = dble(stBQ)
      bnew_slopeq = dble(stBQ_SLOPE)
      bnew_g      = dble(stBG)
      bnew_slopeg = dble(stBG_SLOPE)
      xdel_s      = dble(stX0)
      q20         = dble(stQ02)

c... it_form = 0 : exp(bt) dependence for DVCS-xsection
c...               with b = b_u+b_slopeu * log(q2/2)
c... it_form = 2 : Pauli/Dirac form factors
c... it_form = 4 : exp(bt) for x<xdel_s, P/D otherwise
C --------------------------------------------------------

      idebug = stIDEBUG

C --------------------------------------------------------
      istep  = stIGEN

c... istep = 0 : 1st run for grid calculation
c... istep = 2 : 2nd step for grid
c... istep = 4 : generation only
C --------------------------------------------------------


*---initialize HBOOK:
      CALL HLIMIT(NWPAWC)

*---Open the HBOOK file:
      LUNNT = 31
      LRECL = 1024
      CALL HROPEN(LUNNT, 'toto','bookhis_form_modif.ntp',
     +            'N', LRECL, ISTAT)

      IF (ISTAT.NE. 0) THEN
       PRINT * ,'*** ERROR ',ISTAT,' DURING OPEN OF FILE 31'
       CLOSE(LUNNT)
       mfail = 1
       GOTO 998
      ENDIF

* -> Initialize the HBOOK file
      CALL NTINIT

      NTID= 1

      LUN = 23

      if (istep.ne.0) then
       open(LUN,file='bases.data',status='old',
     +     form='unformatted')
      else
       open(LUN,file='bases.data',status='unknown',
     +     form='unformatted')
      endif


*=========================================================
*          Initialization of BASES by calling BSINIT
*=========================================================
*         -------------
           CALL BSINIT
*         -------------

       if (istep.ne.0) then
        call bsread(lun)
       endif

*=========================================================
*      Initialization of BASES parameters
*=========================================================

      if (srad.eq.0)then

         NDIM  = 3
         NWILD = 3

      else

         NDIM  = 4
         NWILD = 4
      endif


         NCALL = stNCALL

         ITMX1 = stITMX1
         ACC1  = 0.05
         ITMX2 = stITMX2
         ACC2  = 0.01

* X(1) -> x
      xmin=dble(stXMIN)
      xmax=dble(stXMAX)
      if (igen_exp.eq.1) then
      XL(1) = dlog(xmin)
      XU(1) = dlog(xmax)
      else
      XL(1) = xmin
      XU(1) = xmax
      endif

* X(2) -> q
      qmin=dble(stQMIN)
      qmax=dble(stQMAX)
      if (igen_exp.eq.1) then
      XL(2) = dlog(qmin)
      XU(2) = dlog(qmax)
      else
      xl(2) = qmin
      xu(2) = qmax
      endif

* X(3) -> t
      tmax  = dble(stTMAX)
      tmin  = dble(stTMIN)
      XL(3) = tmax
      XU(3) = tmin

* X(4) -> Egamma
      XL(4) = dble(stEIMIN) !0.00001
      XU(4) = dble(abs(stELEP))-5.d0


*=========================================================
*      Initialization of Histograms
*=========================================================


*      CALL XHINIT(1,XL(1),XU(1),40,'log(x)')
*      CALL XHINIT(2,XL(2),XU(2),40,'log(q)')

      CALL XHINIT(1,XL(1),XU(1),40,'x')
      CALL XHINIT(2,XL(2),XU(2),40,'q')
      CALL XHINIT(3,XL(3),XU(3),40,'t')
      CALL XHINIT(4,XL(4),XU(4),40,'Egam radie')


c ipro = 1 : BH
c ipro = 2 : DVCS
c ipro = 3 : BH + DVCS

      ipro =  stIPRO

C
C     initalize value for dvcs subroutine
C
      spin = dble(stSPIN)
      Z    = dble(stZTAR)
      A    = dble(stATAR)

C
C     # of points in x and Q^2
C
      nx = stNX
      nq = stNQ


C
C     initializing flags
C
      icount   = 1
      icountxq = 1
      icountt  = 1

      nset = stNSET

      Angle_Integ = .false.
      if (nset.eq.2) Angle_Integ=.true.


C------- An EIC setting (electron = 5 GeV, target = 200 GeV)
ccc      s = 4.*5.*200.

C------- HERMES setting
ccc      s = 2.*27.55*0.938272D0

C------- HERA setting with 920 GEV protons
ccc      s = 4.*27.55*920.

C------- CLAS setting: CLAS I = 4.3 GEV, CLAS II = 5.75 GeV
ccc      s = 2.*4.3D0*0.938272D0

       stFIXED_save = stFIXED

       mp     = 0.93827d0

      if (stFIXED) then
       s = 2.*abs(stELEP) * mp
      else
       s = 4.*abs(stELEP) * stETARG
      endif

C
C     dummy settings for the angles at the moment
C
      phi   = 0.0
      pphi  = 0.0
      theta = 0.0
C
C     charge, polarization and order in alpha_s
C
      tchar =  dble(stLCHAR)
      xlaml =  dble(stLPOL)
      xlamp =  dble(stTPOL)

C
C     ORDER
C
      iord  = stIORD
      iiord = iord

      icountt  = 1
      icountxq = 1

C
C     TWIST 3
C
      if (stTWIST3) then
      isel_twist3 = 1
      else
      isel_twist3 = 0
      endif

      if (istep.ne.4) then

********************************************************************
*              Nimerical Integration by BASES V5.1
********************************************************************

      if (stIELAS.eq.0) then
      weight_inel = dgausskeps(dmy,dble(stMYMIN),
     +                         dble(stMYMAX),1.d-8)
      else
      weight_inel = 1.d0
      endif
      WRITE(6,*) 'weight_inel=',weight_inel

      CALL BASES( FDVCS, ESTIM, SIGMA, CTIME, IT1, IT2 )

      lu = 66
      open(66,file='bases.out',status='unknown')


      CALL BSINFO( LU )

      CALL BHPLOT( LU )


      WRITE(66,*) '*************************************'
      WRITE(66,*) '*     INTEGRATION TERMINATED        *'
      WRITE(66,*) '*************************************'
      WRITE(6,*) '                                     '
      SUM = ESTIM*weight_inel
      ERR = SIGMA*weight_inel
      WRITE(66,*) 'Total cross-section (nb) : ',SNGL(SUM)
      WRITE(66,*) 'Error                    : ',SNGL(ERR)

      SIGMATOT  = SNGL(SUM)
      DSIGMATOT = SNGL(ERR)


*==================================================
*     Save the probability information to the file
*==================================================

      CALL BSWRIT( LUN )

      endif

998   continue

      return
      end



