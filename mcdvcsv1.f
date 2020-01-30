*CMZ :          13/07/2004  22.55.59  by  H1 Saclay
*-- Author :    Unknown   17/12/2003


****************************************************************
*
       subroutine dvcsob(spin2,Z1,A1,nx1,nq1,nt1,count,x,q,s,t,
     > phi,pphi,theta,ichar,laml,lamp,iord,results,phi_newcalc)
*
****************************************************************
C
C     subroutine which calculates the DVCS cross section + other observables
C

C
C     Definition block
C
cls      implicit none
         implicit real*8(a-h,o-z)

      integer i_outrange
      common/out_of_range/i_outrange

      real*8  it_form
      real*8  bnew_q,bnew_g,xdel_s,bnew_slopeq,bnew_slopeg,q20
      common/form_factors/it_form,bnew_q,bnew_g,xdel_s,q20,
     +                            bnew_slopeq,bnew_slopeg

      real*8 eps_ls

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

C
C     number of points in x and Q^2 in global grid
C
      integer nx,nq,nt,nx1,nq1,nt1
C
C     external counter variable to avoid repeated readin of amplitudes.
C
        integer count
C
C     flag to indicate whether x,Q^2 or t has changed from the last call.
C     if not do not repeat interpolation for the amplitudes.
C
        integer icountxq,icountt
        common /counter1/ icountxq,icountt
C
C     flag to indicate which setting to choose (with or without angles)
C
        integer nset
        common /setting/ nset
C
C     flag to indicate LO or NLO
C
        integer  iord
C
C     external kinematical variables
C
      real*8 s,t,q,x,A1,Z1,spin1,spin2
C
C     angles for transverse polarization & mixing between trans. and long.
C     proton pol. not yet used!
C
      real*8 phi,pphi,theta
C
C     output variables, cross sections
C
C
      real*8 tot,dvcssq,intdvcsbh,bhsq
C
C     charge: e^+ (ichar = 1), e^- (ichar = -1)
C
      real*8 ichar
C
C     polarization flags (lamp not yet used)
C
      real*8 laml,lamp
C
C     various counters
C
      integer i,j,k,k1,k11
C
C     arrays for zeta=x_bj and q
C
      integer mx,mq,mt
      parameter(mx=100,mq=100,mt=100)
      real*8 xx(mx), qq(mq), tt(mt)
C
C     Definition of the various contributions to the diff. Xsection from DVCS,
C     interference and BH
C     up = unpolarized, lp = long. pol, tp=trans. pol Proton.
C

      real*8  tbhup,tiup,tiup1,tilp,tilp1, tilp11,
     > titp,titp1,titp2,titp3,tiup11,titp4,titp5,
     > tbhlp,tbhtp,tdvcslp,tdvcstp,tdvcstp1,tdvcsup,
     > tbhupc1,tbhupc2,tbhups1,tbhtps1,tbhlpc1,tbhtpc1,
     > tdvcsupcos,tdvcsupsin,tiupcos2,tiupsin2,
     > tdvcslpcos,tdvcslpsin,tilpsin2,tilpcos2,
     > tdvcstpcoscos,tdvcstpcossin,tdvcstpsincos,
     > tdvcstpsinsin,titpsinsin,titpsincos,titpcossin,
     > titpcoscos,titpcos2cos,titpcos2sin,titpsin2sin,
     > titpsin2cos

C
C     Kinematical constants and parameters used throughout the program
C
      real*8 dM,tmax,tmin,eps,pi
      real*8 phimin,phimax
      real*8 ymax,ysim, ep1, kfac
      real*8 y,del,qsq,testt,A,Z

C
C     parameters for parametrization of t-dependence (form factors)
C
      real*8 al,kp,kn,ksea,mp,mn,mv,ma,mpi,ga3,mp1
      data al,kp,kn,ksea,mp1,mn,mv,ma,mpi,ga3/7.29735D-3,
     &        1.79,-1.91,0.0,0.93827,0.939565,
     &        0.840,0.900,0.134977,1.267/
      parameter(pi = 3.141592653589793)

C
C     initiation of phi integration for diff. cross section. First define variables
C     for the various pieces of the phi dependence of Xsection
C
      real*8 corint,corcosi,cor,corcos,result,cor1,corint1
      real*8 corcosi2,corcos2,corsini,corsin,corsin2,corsinsq
      real*8 corcossq,corcos2phi,corcos2phi2,corsin2phi,corsin2phi2
      real*8 corsin2p,corsin2p2,corcos2p2

C
C     set up for expressions for diff. Xsection and asymmetries
C
      real*8 asym1,asym2,asym3,difxuw,difx
      real*8 asym11,asym21,asym31,aaauw,ssauw,resa1
      real*8 cauw,upltuw,cadsfluw,fin,resa2,resa3,resa4,resa5,resa6
      real*8 resa7,resa8,resa9,resa10,resa11,resa12,resa13,resa14,resa15
      real*8 resa16,resa17,resa18,resa19,resa20,resa21,resa22,resa23
      real*8 resa24,resa25,resa26,resa27,resa28,resa29,resa30,resa31,
     > resa32,resa33,resa34,resa35,resa36,resa37,resa38,resa39,resa40

C
C     array to store results of calculations
C

      real*8 results(mx)

C
C     common blocks
C

C
C     x and Q^2 arrays
C

       common /kins/ xx,qq,tt

C
C     kinematical variables
C
      common /req/ y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1

C
C     # of x and Q^2
C
      common /counter/ nx,nq,nt

C
C     external functions
C

C
C     porpagator functions of the type: cos(phi)/(P_1*P_2) see BMK
C

      external corcos,cor1,corcos2,corsin,corsin2,corcos2phi,
     > corsin2phi,corcos2phi2,corsin2phi2

C
C     specify readin function
C

      external readin

C
C     Definition of the various contributions to the diff. Xsection from DVCS,
C     interference and BH
C     up = unpolarized, lp = long. pol, tp=trans. pol Proton.
C


      external  tbhup,tiup,tiup1,tilp,tilp1, tilp11,
     > titp,titp1,titp2,titp3,tiup11,titp4,titp5,
     > tbhlp,tbhtp,tdvcslp,tdvcstp,tdvcstp1,tdvcsup,
     > tbhupc1,tbhupc2,tbhups1,tbhtps1,tbhlpc1,tbhtpc1,
     > tdvcsupcos,tdvcsupsin,tiupcos2,tiupsin2,
     > tdvcslpcos,tdvcslpsin,tilpsin2,tilpcos2,
     > tdvcstpcoscos,tdvcstpcossin,tdvcstpsincos,
     > tdvcstpsinsin,titpsinsin,titpsincos,titpcossin,
     > titpcoscos,titpcos2cos,titpcos2sin,titpsin2sin,
     > titpsin2cos

C
C     set up of integration routine
C

ccc      external D01AKF

      INTEGER LW,IW,LIW,IFAIL
      PARAMETER(LW=4000,LIW=LW/4)
      REAL*8 epsabs,epsrel
      real*8 abserr,w
      DIMENSION w(lw),iw(liw)

C     initiate error boundaries for integral routine

      EPSABS = 0.0000001D0
      EPSREL =1.0D-5
      IFAIL =-1

C
C     SET spin of target
C
      spin1 = spin2

C
C     set # in x and q
C

            nx = nx1
            nq = nq1
            nt = nt1

C
C     adjust del for nuclear number A
C

            x  = x/A1
            mp = mp1*A1
            A  = A1
            Z  = Z1
C
C     set x_bj, Q^2 and now dummy var. i, j
C

            qsq = q
            del = x

            i = 1
            j = 1

C
C   Read in amplitudes and t-dependence, depending whether LO (iord=1) or
C   NLO (iord=2). Readin only one time!
C

      if (count.eq.1) then

      call readin(iord,nx,nq,nt)

ccc      if (A.gt.1D0) then
ccc         call readtdep(A)
ccc      endif

      endif

C
C     Calculating parameters depending on x and Q^2!
C
            testt = t

            ep1 = (2.*del*mp)**2D0/qsq

            tmin = - qsq*(2.*(1.-del)*(1. - dsqrt(1+ep1)) + ep1)/
     >           (4.*del*(1.-del) + ep1)


C
C     ymax = y_coll from BMK in order to avoid the singularity in P_1 at y=y_coll
C

caf            ymax = 2.*(dsqrt(1+ep1) - 1.)/ep1

               ymax = (qsq+testt)/(qsq+(A*del)*testt)

               y    = qsq/((A*del)*(s - (mp/A)**2))

C
C     y,t warning if y,t out of ex. range. Skip the calculation for these x, Q^2 value.
C

            if ((y.gt.1.0d0).or.(y.lt.0.d0)) then
ccc               print *,'yval=',y
               i_outrange = 1
               goto 200
            endif

            if ((y.gt.ymax).or.(testt.gt.tmin)) then
               i_outrange = 1
               goto 200
            endif

C
C     calculation of the k-factor from BMK (Belitsky,Mueller, Kirchner)
C

       kfac = - (testt/qsq)*(1.-del)*(1.-y-y**2*ep1/4.)
     >         *(1.-tmin/testt)*(dsqrt(1.+ep1) + (4.*del*(1.-del)+ep1)/
     >          (4.*(1.-del))*(testt-tmin)/qsq)

            if (kfac.lt.0.) then
               write(6,*) 'kfac negative!'
               i_outrange = 1
               goto 200
            endif

C
C     calculation of LO with OLD formulae
C

      if (isel_tintin.eq.1) then

      F2A =  3.1
      F2B =  0.76
      F2C =  0.124
      F2D = -0.188
      F2E = -2.91
      F2F = -0.043
      F2G =  3.69
      F2H =  1.4

      ALPHA = (1.d0/137.d0)

      tval  = -testt
      q2val =  qsq
      xval  =  del
      yval  =  qsq/(del*(s - mp**2))
      seff  =  s

      TAU   = tval/(4d0*(mp**2))
      GE    = (1.+(TVAL/(0.71d0)))**(-2.d0)
      GM    = 2.7d0*GE

      ETA   = (PI/2.)*(0.176+0.033*dlog(Q2VAL))

cc      F2VAL = ((F2A*(XVAL**F2B))+(F2C*(XVAL**F2D)*
cc     &        (1.D0+F2E*DSQRT(XVAL))*
cc     &        (dlog(Q2VAL)+(F2F*((dlog(Q2VAL))**2))+(F2H/Q2VAL))))*
cc     &        ((1.D0-XVAL)**F2G)

c Corrected April 2011--> allm97 like in the paper: hep-ph/9712415v2
      F2VAL=f2allm(xval,q2val)

cc      ytest =  sf2test(real(XVAL),real(Q2VAL))
cc      write(6,*) xval,q2val,f2val,f2val/ytest

      if (isel_f2qcd.eq.1) F2VAL = dble(sf2test(real(XVAL),real(Q2VAL)))

      SIGBH =(((alpha**3.)*SEFF*(yval**2.)*(1.D0+(1.D0-yval)**2.))/
     &     (PI*(Q2VAL**2.)*tval*(1.D0-yval)))*((GE**2.+TAU*(GM**2.))/
     &     (1.D0+TAU))


      results(4) = 2d0*pi* SIGBH * yval/q2val * 3.8937966E+5
ccc
c S. Fazio -  October 2010 - include the gamma-p xsec
c flux = alpha/(2*pi*q2)*(1+(1-y)**2)/y 
c      FLUX =(ALPHA/(2.D0*PI*Q2VAL))*(1.D0+(1.D0-YVAL)**2)/YVAL

c xsec e-p is:
      SIGDVCS =((PI*(ALPHA**3.)*SEFF)/(4.D0*(RTI**2.)*Q2VAL**3.))*
     &     (1.D0+(1.D0-YVAL)**2.)*(dexp(-1.D0*BTI*TVAL))*(F2VAL**2.)*
     &     (1.D0+ETA**2.)

c      SIGDVCS = SIGDVCS/FLUX
c      SIGDVCS = SIGDVCS * xval/yval
c      SIGDVCS = SIGDVCS * xval

      results(2) = 2d0*pi* SIGDVCS * yval/q2val * 3.8937966E+5

c xsec gamma-p is:
cc      SIGDVCS =PI**3.*(ALPHA**2./RTI**2.
cc     &     *1.D0/Q2VAL**2.
cc     &     *(dexp(-1.D0*BTI*TVAL))*(F2VAL**2.)
cc     &     *(1.D0+ETA**2.)

cc      results(2) = SIGDVCS * 3.8937966E+5

ccc
      if (isel_dipole.eq.1) then

       Q20dip    =   1.00
       W0dip     =  20.00
       const1dip =   0.16
       const2dip =   0.84
       const3dip = 311.00
       const4dip =   5.33
       const5dip =   1.33
       const6dip =   5.37
       const7dip =  49.42
       const8dip =   7.65
       const9dip =   4.94

      WVAL  = dsqrt(seff*yval)

      rsoft = (const4dip
     &        -const5dip*dexp(-4.*(Q2VAL/Q20dip))
     &        +const6dip*(Q2VAL/Q20dip))**(-1.)
      rhard = (const7dip
     &        -const8dip*dexp(-4.*(Q2VAL/Q20dip))
     &        +const9dip*(Q2VAL/Q20dip))**(-1.)
      sigdip=  (1./(16*pi))
     &        *(((rsoft*((WVAL/W0dip)**(const1dip)))
     &          +(rhard*((WVAL/W0dip)**(const2dip))))**2)
     &        *const3dip
     &        *1000.

      sigdip=sigdip*ALPHA*((1.D0+(1.D0-YVAL)**2.)/2.d0)*
     +                     (dexp(-1.D0*BTI*TVAL))
      sigdip=sigdip/(pi*yval*q2val)
      sigdip=sigdip * yval / xval

c S. Fazio -  October 2010 - include the gamma-p xsec
c flux = alpha/(2*pi*q2)*(1+(1-y)**2)/y 
cc      FLUX =(ALPHA/(2.D0*PI*Q2VAL))*(1.D0+(1.D0-YVAL)**2)/YVAL
cc      sigdip=sigdip/FLUX
c      write(6,*) 'flux! :',FLUX

      results(2) = sigdip

      endif ! isel_dipole=1

C
C     PHIRVAL=(had,lep) in HCMS = (rndm) --> 0
C
      SIGINT = -1.D0*SIGN*((ETA*(ALPHA**3.)*SEFF*YVAL*
     &         (1.D0+(1.D0-YVAL)**2.)*
     &          1.D0*(dexp(-1.D0*BTI*TVAL/2.))*F2VAL)/
     &         (2.D0*(Q2VAL**(5.D0/2.D0))*DSQRT(TVAL)*
     &         (DSQRT(1.D0-YVAL))*RTI))*
     &         ((GE+TAU*GM)/(1.D0+TAU)) *DCOS(phi_newcalc)

cls      SIGINT = 0d0 ! in case of integration on angles...
cls                  => no cst int here!!

      results(3) = 2d0*pi* SIGINT * yval/q2val * 3.8937966E+5

ccc      print *,results(2),results(3)

      goto 200

      endif

C
C     compute diff cross section unpol., pol. ,transversely pol.
C

C
C     First option: all the angles are integrated out and the results
C     for the total cross section, DVCS, Interference, BH as well as
C     the results for the SSA and CA with and without twist-3 are stored in
C     the array results. Note here only unpol. target and probe except for the SSA.
C     Second option: most general option with all angles and polarizations free.
C

            if (nset.eq.2) then

C
C     phi integrals
C
      phimin = 0.0
      phimax = 2.*pi

      eps_ls = 10.0d0**(-8.0d0)

C
C     integral over phi from 0 to 2*pi of 1/P_1*P_2
C

ccc      call D01AKF(cor1,phimin,phimax,epsabs,epsrel,result,abserr,w,
ccc     &            lw,iw,liw,ifail)
ccc      corint1 = result

         corint1 = dgausskeps(cor1,phimin,phimax,eps_ls)

ccc      print *,'1/P_1*P_2 all=',corint1
C-----test
ccc      xbi = 2d0*(2d0-y)*dsqrt(-testt*(1d0-y)*
ccc     +          (1d0-del)/qsq)/(1d0-y)
ccc      if (xbi.ge.0.9d0) print *,'WARNING:t,q2,x,A,1/P=>',testt,qsq,del,
ccc     +                  xbi,corint1
C-----test
ccc      term_a = 2d0*3.1415d0*dsqrt((1d0+xbi)/(1d0-xbi))/
ccc     +         (1d0+xbi)
ccc      term_a = term_a * y*2/(1d0-y)
ccc      print *,'1/P1P2 new=',term_a
cls-----------------------
CSF      xv_max = 1000.d0
cls-----------------------
ccc      if (dabs(corint1).gt.100.d0) then
ccc               print *,'t,q2,x,1/P (>100) =',testt,qsq,del,corint1
ccc      endif

       if (dabs(corint1).gt.xv_max) then
CSF               print *,'t,q2,x,1/P (>1000)=',testt,qsq,del,corint1
CSF                i_outrange = 1
CSF                goto 200
       endif

C
C     integral over phi from 0 to 2*pi of cos(phi)/P_1*P_2
C

ccc      call D01AKF(corcos,phimin,phimax,epsabs,epsrel,result,abserr,w,
ccc     &            lw,iw,liw,ifail)
ccc      corcosi = result

CLS      corcosi = dgausskeps(corcos,phimin,phimax,eps_ls)
         corcosi = 0.d0 ! this is a good approximation!
Cls--------------------------------------new
         corcosi = dgausskeps(corcos,phimin,phimax,eps_ls)
Cls--------------------------------------new

C
C     integral over phi from 0 to 2*pi of cos(2phi)/P_1*P_2
C

ccc      call D01AKF(corcos2phi,phimin,phimax,epsabs,epsrel,result,abserr
ccc     &            ,w,lw,iw,liw,ifail)
ccc      corcosi2 = result

CLS      corcosi2 = dgausskeps(corcos2phi,phimin,phimax,eps_ls)
         corcosi2 = 0.d0 ! this is a good approximation!

C
C     integral over phi from 0 to 2*pi of sin(2phi)/P_1*P_2
C

ccc      call D01AKF(corsin2phi,phimin,phimax,epsabs,epsrel,result,abserr
ccc     &            ,w,lw,iw,liw,ifail)
ccc      corsin2p = result
         corsin2p = 0d0

C
C     integral over phi from 0 to 2*pi of sin(phi)/P_1*P_2
C

ccc      call D01AKF(corsin,phimin,phimax,epsabs,epsrel,result,abserr,w,
ccc     &            lw,iw,liw,ifail)
ccc      corsini = result
         corsini = 0d0


C
C     compute dM from phase space (1/dx.dqsq.dt)
C

      dM = al**3*del*y**2/(16.*pi**2*qsq**2*
     >    (1. + 4.*del**2*mp**2/qsq)**0.5*A)


C
C     unpolarized diff. Xsection (unpol. probe/unpol. target)
C

C
C     const. DVCS^2
C
      resa1 = 4.*pi**2*tdvcsup(1D0,1D0,j,i)*dM

C
C    t does not change anymore
C
      icountt = 2


cls>> no interference a O(1/Q)
C
C     cos(phi) interference
C
cint      resa2 = ichar*2.*pi*corcosi*tiup(1D0,1D0,j,i)*dM
          resa2 = 0.d0
Cls--------------------------------------new
          resa2 = ichar*2.*pi*corcosi*tiup(1D0,1D0,j,i)*dM
Cls--------------------------------------new

C
C     const. interference
C
          resa4 = ichar*2.*pi*corint1*tiup11(1D0,1D0,j,i)*dM

C
C     cos(2phi) interference
C
cint      resa23 = ichar*2.*pi*corcosi2*tiupcos2(1D0,1D0,j,i)*dM
          resa23 = 0.d0

C
C     sin(phi) interference
C
cint      resa11 = ichar*2.*pi*corsini*tiup1(laml,1D0,j,i)*dM
          resa11 = 0.d0

C
C     sin(2phi) term interference
C
cint      resa25 = ichar*2.*pi*corsin2p*tiupsin2(laml,1D0,j,i)*dM
          resa25 = 0.d0

C
C     const. BH^2
C
      resa3 = 2.*pi*corint1*tbhup(1D0,1D0,j,i)*dM

C
C     cos(phi) BH^2
C
CLS      resa5 = 2.*pi*corcosi*tbhupc1(1D0,1D0,j,i)*dM
         resa5 = 0.d0
C
C     cos(2phi) BH^2
C
CLS      resa6 = 2.*pi*corcosi2*tbhupc2(1D0,1D0,j,i)*dM
         resa6 = 0.


C
C     long. pol. probe + unpol. target, diff Xsection
C


************************ final results ************************

C
C     collect final result
C

C
C     first DVCS^2
C
      q2val =  qsq
      ETA   = (PI/2.)*(0.176+0.033*dlog(Q2VAL))
ccc      ETA   = 0.d0

      dvcssq     = resa1
      if (isel_realpart.eq.1) then
       dvcssq     = dvcssq * (1.D0+ETA**2.)
      endif

CLS ICI4 
ccc      print *,eta**2

      results(2) = dvcssq*3.8937966E+5

C
C     interference
C

      intdvcsbh  = resa4 + resa2 + resa23 + resa11 + resa25
      if (isel_interf.eq.0) intdvcsbh = 0d0
      results(3) = intdvcsbh*3.8937966E+5

ccc         intdvcsbh  = resa4 ! keep only the cst term
ccc         intdvcsbh  = resa4 + resa2! keep only the cst+cos terms
ccc         results(3) = intdvcsbh*3.8937966E+5 ! S.F. Dec. 2012
ccc         results(3) = intdvcsbh

C
C     BH^2
C

      bhsq       = resa3 + resa5  + resa6 ! only resa3
      results(4) = bhsq*3.8937966E+5

cc>>    OLD > for BH checks...
cc      ALPHA = (1.d0/137.d0)
cc      tval  = -testt
cc      q2val =  qsq
cc      xval  =  del
cc      yval  =  qsq/(del*(s - mp**2))
cc      seff  =  s
cc      TAU   = tval/(4d0*(mp**2))
cc      GE    = (1.+(TVAL/(0.71d0)))**(-2.d0)
cc      GM    =  2.7d0*GE
cc      SIGBH = (((alpha**3.)*SEFF*(yval**2.)*(1.D0+(1.D0-yval)**2.))/
cc     &       (PI*(Q2VAL**2.)*tval*(1.D0-yval)))*((GE**2.+TAU*(GM**2.))/
cc     &       (1.D0+TAU))
cc      SIGBH = 2d0*pi* SIGBH * yval/q2val


C
C     total Xsection
C

      tot = dvcssq + intdvcsbh + bhsq
      results(1) = tot*3.8937966E+5
C
C     total Xsection (twist-3)
C

ccc      results(13) = (resa1+resa4+resa2+resa3+resa5+resa6)
ccc     > *3.8937966E+5

      if (isel_asym.eq.0) goto 200

************************ final results ************************


C
C     integral over phi from 0 to 2*pi of cos(phi)^2/P_1*P_2
C

ccc      call D01AKF(corcos2,phimin,phimax,epsabs,epsrel,result,abserr,w,
ccc     &            lw,iw,liw,ifail)
ccc      corcossq = result
         corcossq = dgausskeps(corcos2,phimin,phimax,eps_ls)

C
C     integral over phi from 0 to 2*pi of sin(phi)^2/P_1*P_2
C

ccc      call D01AKF(corsin2,phimin,phimax,epsabs,epsrel,result,abserr,w,
ccc     &            lw,iw,liw,ifail)
ccc      corsinsq = result
         corsinsq = dgausskeps(corsin2,phimin,phimax,eps_ls)

C
C     integral over phi from 0 to 2*pi of cos(2phi)*cos(phi)/P_1*P_2
C

ccc      call D01AKF(corcos2phi2,phimin,phimax,epsabs,epsrel,result,abserr
ccc     &            ,w,lw,iw,liw,ifail)
ccc      corcos2p2 = result
         corcos2p2 = dgausskeps(corcos2phi2,phimin,phimax,eps_ls)

C
C     integral over phi from 0 to 2*pi of sin(2phi)*sin(phi)/P_1*P_2
C

ccc      call D01AKF(corsin2phi2,phimin,phimax,epsabs,epsrel,result,abserr
ccc     &            ,w,lw,iw,liw,ifail)
ccc      corsin2p2 = result
         corsin2p2 = dgausskeps(corsin2phi2,phimin,phimax,eps_ls)

C
C     cos^2(phi) term in CA
C

      resa20 = ichar*2.*pi*corcossq*tiup(1D0,1D0,j,i)*dM

C
C     cos(phi) term in CA
C

CLS      resa21 = ichar*2.*pi*corcosi*tiup11(1D0,1D0,j,i)*dM
         resa21 = 0d0

C
C     cos(phi)*cos(2phi) term in CA
C
      resa24 = ichar*2.*pi*corcos2p2*tiupcos2(1D0,1D0,j,i)*dM

C
C     sin^2(phi) term in SSA
C

      resa22 = ichar*2.*pi*corsinsq*tiup1(laml,1D0,j,i)*dM

C
C     sin(phi)*sins(2phi) term in SSA
C
      resa26 = ichar*2.*pi*corsin2p2*tiupsin2(laml,1D0,j,i)*dM

C
C     full CA
C

      results(5) = 2.*(resa20+(resa21+resa24))*3.8937966E+5
      results(6) = (dvcssq+bhsq)*3.8937966E+5

C
C     CA without WW Twist-3
C

      results(7) = 2.*(resa20+resa21)*3.8937966E+5
      results(8) = (dvcssq+bhsq)*3.8937966E+5
      if (rndm(1).lt.0.001) print *,'results(7)=',results(7)

CLS   note :  resa21~0
CLS   note :  once integrated over x/t =>
CLS   note :  *  ca= results(7)/results(8) *
CLS   => generate with results(7)!!

C
C     full SSA
C

      results(9) = 2.*(resa22+resa26)*3.8937966E+5
      results(10)= (dvcssq+bhsq+resa2+resa4+resa23)*3.8937966E+5

C
C     SSA without WW Twist-3
C

      results(11)= 2.*(resa22)*3.8937966E+5
      results(12) = (dvcssq+bhsq+resa2+resa4)*3.8937966E+5

CLS   note :  resa2~0 ; resa4 = cst interf.
CLS   note :  once integrated over x/t =>
CLS   note :  * ssa = results(11)/results(12) *
CLS   => generate with results(11)!!

      else ! nset=1 >>

C
C     unpol. cos(phi) interference
C
      results(14) = ichar*corcos(phi)*tiup(1D0,1D0,j,i)*dM


C
C    t does not change anymore
C

      icountt = 2

C
C     unpol. const. interference
C
      results(13) = ichar*cor1(phi)*tiup11(1D0,1D0,j,i)*dM

C
C     unpol. cos(2phi) interference
C
      results(15) = ichar*corcos2(phi)*tiupcos2(1D0,1D0,j,i)*dM

C
C     const. BH^2
C
      results(33) = cor1(phi)*tbhup(1D0,1D0,j,i)*dM

C
C     cos(phi) BH^2
C
      results(34) = corcos(phi)*tbhupc1(1D0,1D0,j,i)*dM

C
C     cos(2phi) BH^2
C
      results(35) = corcos2(phi)*tbhupc2(1D0,1D0,j,i)*dM

C
C     const. DVCS^2
C
      results(1) = tdvcsup(1D0,1D0,j,i)*dM

C
C     cos(phi) DVCS^2
C
      results(2) = dcos(phi)*tdvcsupcos(1D0,1D0,j,i)*dM

C
C     long. pol. probe + unpol. target, diff Xsection
C

C
C     sin(phi) interference
C
      results(16) = ichar*corsin(phi)*tiup1(laml,1D0,j,i)*dM

C
C     sin(2phi) interference
C
      results(17) = ichar*corsin2(phi)*tiupsin2(laml,1D0,j,i)*dM

C
C     sin(phi) DVCS^2
C
      results(3) = dsin(phi)*tdvcsupsin(laml,1D0,j,i)*dM

C
C     long.pol. probe + long.pol. target
C

C
C     const. interference
C
      results(18) = ichar*cor1(phi)*tilp11(laml,lamp,j,i)*dM

C
C     cos(phi) interference
C
      results(19) = ichar*corcos(phi)*tilp(laml,lamp,j,i)*dM

C
C     sin(phi) interference
C
      results(20) = ichar*corsin(phi)*tilp1(1D0,lamp,j,i)*dM

C
C     cos(2phi) interference
C
      results(21) = ichar*corcos2(phi)*tilpcos2(laml,lamp,j,i)*dM

C
C     sin(2phi) interference
C
      results(22) = ichar*corsin2(phi)*tilpsin2(1D0,lamp,j,i)*dM

C
C     const. BH^2
C
      results(36) = cor1(phi)*tbhlp(laml,lamp,j,i)*dM

C
C     cos(phi) BH^2
C
      results(37) = corcos(phi)*tbhlpc1(laml,lamp,j,i)*dM

C
C     const. DVCS^2
C
      results(4) = tdvcslp(laml,lamp,j,i)*dM

C
C     cos(phi) DVCS^2
C
      results(5) = dcos(phi)*tdvcslpcos(laml,lamp,j,i)*dM

C
C     sin(phi) DVCS^2
C
      results(6) = dsin(phi)*tdvcslpsin(1D0,lamp,j,i)*dM

C
C     unpol./pol probe + trans. pol. target
C

C
C     const.*sin(pphi) interference
C
      results(23) = ichar*dsin(pphi)*cor1(phi)*titp1(1D0,1D0,j,i)*dM

C
C     const.*cos(pphi) interference
C
      results(24) = ichar*dcos(pphi)*cor1(phi)*titp(laml,1D0,j,i)*dM

C
C     cos(phi)*cos(pphi) interference
C
      results(25) = ichar*dcos(pphi)*corcos(phi)*titp2(laml,1D0,j,i)*dM

C
C     cos(phi)*sin(pphi) interference
C
      results(26) = ichar*dsin(pphi)*corcos(phi)*titp3(1D0,1D0,j,i)*dM

C
C     sin(phi)*cos(pphi) interference
C
      results(27) = ichar*dcos(pphi)*corsin(phi)*titp4(1D0,1D0,j,i)*dM

C
C     sin(phi)*sin(pphi) interference
C
      results(28) = ichar*dsin(pphi)*corsin(phi)*titp5(laml,1D0,j,i)*dM

C
C     cos(2phi)*cos(pphi) interference
C
      results(29) = ichar*dcos(pphi)*corcos2(phi)*
     >        titpcos2cos(laml,1D0,j,i)*dM

C
C     cos(2phi)*sin(pphi) interference
C
      results(30) = ichar*dsin(pphi)*corcos2(phi)*
     >        titpcos2sin(1D0,1D0,j,i)*dM

C
C     sin(2phi)*cos(pphi) interference
C
      results(31) = ichar*dcos(pphi)*corsin2(phi)*
     >        titpsin2cos(1D0,1D0,j,i)*dM

C
C    sin(2phi)*sin(pphi) interference
C
      results(32) = ichar*dsin(pphi)*corsin2(phi)*
     >        titpsin2sin(laml,1D0,j,i)*dM

C
C     const*cos(pphi) BH^2
C
      results(38) = dcos(pphi)*cor1(phi)*tbhtp(laml,1D0,j,i)*dM

C
C     cos(phi)*cos(pphi) BH^2
C
      results(39) = dcos(pphi)*corcos(phi)*tbhtpc1(laml,1D0,j,i)*dM

C
C     sin(phi)*sin(pphi) BH^2
C
      results(40) = dsin(pphi)*corsin(phi)*tbhtps1(laml,1D0,j,i)*dM

C
C     const.*cos(pphi) DVCS^2
C
      results(7) = dcos(pphi)*tdvcstp(laml,1D0,j,i)*dM

C
C     const.*sin(pphi) DVCS^2
C
      results(8) = dsin(pphi)*tdvcstp1(laml,1D0,j,i)*dM

C
C     cos(phi)*cos(pphi) DVCS^2
C
      results(9) = dcos(phi)*dcos(pphi)*tdvcstpcoscos(laml,1D0,j,i)*dM

C
C     cos(phi)*sin(pphi) DVCS^2
C
      results(10) = dcos(phi)*dsin(pphi)*tdvcstpcossin(1D0,1D0,j,i)*dM

C
C     sin(phi)*cos(pphi) DVCS^2
C
      results(11) = dsin(phi)*dcos(pphi)*tdvcstpsincos(1D0,1D0,j,i)*dM

C
C     sin(phi)*sin(pphi) DVCS^2
C
      results(12) = dsin(phi)*dsin(pphi)*tdvcstpsinsin(laml,1D0,j,i)*dM

C
C     collect final result
C

C
C     first DVCS^2
C

      results(41) = (results(1) + results(2) + results(3)
     > + dcos(theta)*(results(4) + results(5) + results(6)
     > ) + dsin(theta)*(results(7) + results(8) + results(9)
     >+ results(10) + results(11) + results(12)))*3.8937966E+5
      print *,results(41)


C
C     interference
C

      results(42) = (results(13) + results(14) + results(15)
     > + results(16) + results(17)
     > + dcos(theta)*(results(18) + results(19) + results(20)
     > + results(21) + results(22))
     > + dsin(theta)*(results(23) + results(24) + results(25)
     > + results(26) + results(27) + results(28)+ results(29)
     > + results(30) + results(31) + results(32)))*3.8937966E+5


C
C     BH^2
C

      results(43) = (results(33) + results(34) + results(35)
     > + dcos(theta)*(results(36) + results(37))
     > + dsin(theta)*(results(38) + results(39) + results(40)))
     > *3.8937966E+5



C
C     total Xsection
C

      results(44) = results(41) + results(42) + results(43)

         endif ! nset=1,2


         return


 200     continue

ccc         print *,"OUT OF RANGE!"," y=",y," ,ymax=",ymax," ,t=",testt,
ccc     >   " ,tmin=",tmin


      return

      end

CMM   Input data files which contain re and imag parts of DVCS amplitudes
CMM   for the appropriate flavour and order
CMM   zeta(j=1)
CMM   Q(i) ReA(j,i) ImA(j.i)
CMM   ...
CMM   zeta(j=nx)
CMM   ...
CMM   Q(nq) ReA(nx,nq) ImA(nx.nq)
C
C     at the moment only H! Negelect others!
C
****************************************************************
      subroutine readin(nord,nx,nq,nt)
****************************************************************

      implicit none

      integer nord
      integer nx,nq,nt,mx,mq,mt
      parameter(mx=100,mq=100,mt=100)
      real*8 xx(mx),qq(mq),tt(mt)
      real*8 dum1,dum2
      integer i,j,k
      real*8 x,q,t
      real*8 y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1
      common /kins/ xx,qq,tt
      common /req/ y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1

      interface
        function readinsingle(nx,nq,nt,path,isim) result(arr)
            integer nx,nq,nt
            character(len=*) path
            logical isim
            real*8 arr(nx,nq,nt)
        end function
      end interface

      real*8 reu(mx,mq,mt),imu(mx,mq,mt),red(mx,mq,mt),imd(mx,mq,mt)
      real*8 res(mx,mq,mt),ims(mx,mq,mt),reg(mx,mq,mt),img(mx,mq,mt)
      real*8 reup(mx,mq,mt),imup(mx,mq,mt),redp(mx,mq,mt),imdp(mx,mq,mt)
      real*8 resp(mx,mq,mt),imsp(mx,mq,mt),regp(mx,mq,mt),imgp(mx,mq,mt)
      real*8 reue(mx,mq,mt),imue(mx,mq,mt),rede(mx,mq,mt),imde(mx,mq,mt)
      real*8 rese(mx,mq,mt),imse(mx,mq,mt),rege(mx,mq,mt),imge(mx,mq,mt)
      real*8 reuep(mx,mq,mt),imuep(mx,mq,mt),redep(mx,mq,mt),
     > imdep(mx,mq,mt)
      real*8 resep(mx,mq,mt),imsep(mx,mq,mt),regep(mx,mq,mt),
     > imgep(mx,mq,mt)

      real*8 reut3(mx,mq,mt),imut3(mx,mq,mt),redt3(mx,mq,mt)
      real*8 rest3(mx,mq,mt),imst3(mx,mq,mt),imdt3(mx,mq,mt)
      real*8 reupt3(mx,mq,mt),imupt3(mx,mq,mt),redpt3(mx,mq,mt)
      real*8 respt3(mx,mq,mt),imspt3(mx,mq,mt),imdpt3(mx,mq,mt)
      real*8 reuet3(mx,mq,mt),imuet3(mx,mq,mt),redet3(mx,mq,mt)
      real*8 reset3(mx,mq,mt),imset3(mx,mq,mt),imdet3(mx,mq,mt)
      real*8 reuept3(mx,mq,mt),imuept3(mx,mq,mt),redept3(mx,mq,mt)
      real*8 resept3(mx,mq,mt),imsept3(mx,mq,mt),imdept3(mx,mq,mt)

      real*8 reut3d(mx,mq,mt),imut3d(mx,mq,mt),redt3d(mx,mq,mt)
      real*8 rest3d(mx,mq,mt),imst3d(mx,mq,mt),imdt3d(mx,mq,mt)
      real*8 reupt3d(mx,mq,mt),imupt3d(mx,mq,mt),redpt3d(mx,mq,mt)
      real*8 respt3d(mx,mq,mt),imspt3d(mx,mq,mt),imdpt3d(mx,mq,mt)
      real*8 reuet3d(mx,mq,mt),imuet3d(mx,mq,mt),redet3d(mx,mq,mt)
      real*8 reset3d(mx,mq,mt),imset3d(mx,mq,mt),imdet3d(mx,mq,mt)
      real*8 reuept3d(mx,mq,mt),imuept3d(mx,mq,mt),redept3d(mx,mq,mt)
      real*8 resept3d(mx,mq,mt),imsept3d(mx,mq,mt),imdept3d(mx,mq,mt)

      real*8 reul(mx,mq,mt),imul(mx,mq,mt),redl(mx,mq,mt),imdl(mx,mq,mt)
      real*8 resl(mx,mq,mt),imsl(mx,mq,mt),regl(mx,mq,mt),imgl(mx,mq,mt)
      real*8 reupl(mx,mq,mt),imupl(mx,mq,mt),redpl(mx,mq,mt),
     > imdpl(mx,mq,mt)
      real*8 respl(mx,mq,mt),imspl(mx,mq,mt),regpl(mx,mq,mt),
     > imgpl(mx,mq,mt)
      real*8 reuel(mx,mq,mt),imuel(mx,mq,mt),redel(mx,mq,mt),
     > imdel(mx,mq,mt)
      real*8 resel(mx,mq,mt),imsel(mx,mq,mt),regel(mx,mq,mt),
     > imgel(mx,mq,mt)
      real*8 reuepl(mx,mq,mt),imuepl(mx,mq,mt),redepl(mx,mq,mt),
     > imdepl(mx,mq,mt)
      real*8 resepl(mx,mq,mt),imsepl(mx,mq,mt),regepl(mx,mq,mt),
     > imgepl(mx,mq,mt)

      common /ampsn/reul,imul,redl,imdl,resl,imsl,regl,imgl,reupl,imupl,
     > redpl,imdpl,respl,imspl,regpl,imgpl,reuel,imuel,redel,imdel,resel
     > ,imsel,regel,imgel,reuepl,imuepl,redepl,imdepl,resepl,imsepl
     > ,regepl,imgepl

      common /amps/ reu,imu,red,imd,res,ims,reg,img,reup,imup,
     > redp,imdp,resp,imsp,regp,imgp,reue,imue,rede,imde,rese,imse
     > ,rege,imge,reuep,imuep,redep,imdep,resep,imsep,regep,imgep

      common /ampst3d/ reut3d,imut3d,redt3d,imdt3d,rest3d,imst3d
     >  ,reupt3d,imupt3d,redpt3d,imdpt3d,respt3d,imspt3d
     >  ,reuet3d,imuet3d,redet3d,imdet3d,reset3d,imset3d
     >  ,reuept3d,imuept3d,redept3d,imdept3d,resept3d,imsept3d

      common /ampst3/ reut3,imut3,redt3,imdt3,rest3,imst3
     >  ,reupt3,imupt3,redpt3,imdpt3,respt3,imspt3
     >  ,reuet3,imuet3,redet3,imdet3,reset3,imset3
     >  ,reuept3,imuept3,redept3,imdept3,resept3,imsept3

C     READ KINEMATIC VARIABLES

      open(unit=11,file='luamp.dat',status='unknown')

      do j = 1, nx
            read(11,*) x

            xx(j) = x

         do i = 1, nq
            read(11,*) q
            
            qq(i) = q
            
            do k = 1, nt
                read(11,*) t,dum1,dum2

                tt(k) = t

            enddo
         enddo
      enddo

      close(unit=11)


      if (spin1.eq.0D0) then

        if (nord.eq.1) then

           reu = readinsingle(nx,nq,nt,'luamp.dat',.false.)
           imu = readinsingle(nx,nq,nt,'luamp.dat',.true.)
           reul = readinsingle(nx,nq,nt,'luamp.dat',.false.)
           imul = readinsingle(nx,nq,nt,'luamp.dat',.true.)

           red = readinsingle(nx,nq,nt,'ldamp.dat',.false.)
           imd = readinsingle(nx,nq,nt,'ldamp.dat',.true.)
           redl = readinsingle(nx,nq,nt,'ldamp.dat',.false.)
           imdl = readinsingle(nx,nq,nt,'ldamp.dat',.true.)

           res = readinsingle(nx,nq,nt,'lsamp.dat',.false.)
           ims = readinsingle(nx,nq,nt,'lsamp.dat',.true.)
           resl = readinsingle(nx,nq,nt,'lsamp.dat',.false.)
           imsl = readinsingle(nx,nq,nt,'lsamp.dat',.true.)

        else

           reu = readinsingle(nx,nq,nt,'nlouamp.dat',.false.)
           imu = readinsingle(nx,nq,nt,'nlouamp.dat',.true.)
           reul = readinsingle(nx,nq,nt,'luamp.dat',.false.)
           imul = readinsingle(nx,nq,nt,'luamp.dat',.true.)

           red = readinsingle(nx,nq,nt,'nlodamp.dat',.false.)
           imd = readinsingle(nx,nq,nt,'nlodamp.dat',.true.)
           redl = readinsingle(nx,nq,nt,'ldamp.dat',.false.)
           imdl = readinsingle(nx,nq,nt,'ldamp.dat',.true.)

           res = readinsingle(nx,nq,nt,'nlosamp.dat',.false.)
           ims = readinsingle(nx,nq,nt,'nlosamp.dat',.true.)
           resl = readinsingle(nx,nq,nt,'lsamp.dat',.false.)
           imsl = readinsingle(nx,nq,nt,'lsamp.dat',.true.)

      endif

           reut3 = readinsingle(nx,nq,nt,'luamptw3.dat',.false.)
           imut3 = readinsingle(nx,nq,nt,'luamptw3.dat',.true.)
           reut3d = readinsingle(nx,nq,nt,'luamptw3d.dat',.false.)
           imut3d = readinsingle(nx,nq,nt,'luamptw3d.dat',.true.)

           redt3 = readinsingle(nx,nq,nt,'ldamptw3.dat',.false.)
           imdt3 = readinsingle(nx,nq,nt,'ldamptw3.dat',.true.)
           redt3d = readinsingle(nx,nq,nt,'ldamptw3d.dat',.false.)
           imdt3d = readinsingle(nx,nq,nt,'ldamptw3d.dat',.true.)

           rest3 = readinsingle(nx,nq,nt,'lsamptw3.dat',.false.)
           imst3 = readinsingle(nx,nq,nt,'lsamptw3.dat',.true.)
           rest3d = readinsingle(nx,nq,nt,'lsamptw3d.dat',.false.)
           imst3d = readinsingle(nx,nq,nt,'lsamptw3d.dat',.true.)

      elseif(spin1.eq.1D0) then

         if (nord.eq.1) then

           reu = readinsingle(nx,nq,nt,'luamp.dat',.false.)
           imu = readinsingle(nx,nq,nt,'luamp.dat',.true.)
           reul = readinsingle(nx,nq,nt,'luamp.dat',.false.)
           imul = readinsingle(nx,nq,nt,'luamp.dat',.true.)

           red = readinsingle(nx,nq,nt,'ldamp.dat',.false.)
           imd = readinsingle(nx,nq,nt,'ldamp.dat',.true.)
           redl = readinsingle(nx,nq,nt,'ldamp.dat',.false.)
           imdl = readinsingle(nx,nq,nt,'ldamp.dat',.true.)

           res = readinsingle(nx,nq,nt,'lsamp.dat',.false.)
           ims = readinsingle(nx,nq,nt,'lsamp.dat',.true.)
           resl = readinsingle(nx,nq,nt,'lsamp.dat',.false.)
           imsl = readinsingle(nx,nq,nt,'lsamp.dat',.true.)


           reup = readinsingle(nx,nq,nt,'luamppol.dat',.false.)
           imup = readinsingle(nx,nq,nt,'luamppol.dat',.true.)
           reupl = readinsingle(nx,nq,nt,'luamppol.dat',.false.)
           imupl = readinsingle(nx,nq,nt,'luamppol.dat',.true.)

           redp = readinsingle(nx,nq,nt,'ldamppol.dat',.false.)
           imdp = readinsingle(nx,nq,nt,'ldamppol.dat',.true.)
           redpl = readinsingle(nx,nq,nt,'ldamppol.dat',.false.)
           imdpl = readinsingle(nx,nq,nt,'ldamppol.dat',.true.)

           resp = readinsingle(nx,nq,nt,'lsamppol.dat',.false.)
           imsp = readinsingle(nx,nq,nt,'lsamppol.dat',.true.)
           respl = readinsingle(nx,nq,nt,'lsamppol.dat',.false.)
           imspl = readinsingle(nx,nq,nt,'lsamppol.dat',.true.)


           reue = readinsingle(nx,nq,nt,'luampe.dat',.false.)
           imue = readinsingle(nx,nq,nt,'luampe.dat',.true.)
           reuel = readinsingle(nx,nq,nt,'luampe.dat',.false.)
           imuel = readinsingle(nx,nq,nt,'luampe.dat',.true.)

           rede = readinsingle(nx,nq,nt,'ldampe.dat',.false.)
           imde = readinsingle(nx,nq,nt,'ldampe.dat',.true.)
           redel = readinsingle(nx,nq,nt,'ldampe.dat',.false.)
           imdel = readinsingle(nx,nq,nt,'ldampe.dat',.true.)

           rese = readinsingle(nx,nq,nt,'lsampe.dat',.false.)
           imse = readinsingle(nx,nq,nt,'lsampe.dat',.true.)
           resel = readinsingle(nx,nq,nt,'lsampe.dat',.false.)
           imsel = readinsingle(nx,nq,nt,'lsampe.dat',.true.)


           reuep = readinsingle(nx,nq,nt,'luamppole.dat',.false.)
           imuep = readinsingle(nx,nq,nt,'luamppole.dat',.true.)
           reuepl = readinsingle(nx,nq,nt,'luamppole.dat',.false.)
           imuepl = readinsingle(nx,nq,nt,'luamppole.dat',.true.)

           redep = readinsingle(nx,nq,nt,'ldamppole.dat',.false.)
           imdep = readinsingle(nx,nq,nt,'ldamppole.dat',.true.)
           redepl = readinsingle(nx,nq,nt,'ldamppole.dat',.false.)
           imdepl = readinsingle(nx,nq,nt,'ldamppole.dat',.true.)

           resep = readinsingle(nx,nq,nt,'lsamppole.dat',.false.)
           imsep = readinsingle(nx,nq,nt,'lsamppole.dat',.true.)
           resepl = readinsingle(nx,nq,nt,'lsamppole.dat',.false.)
           imsepl = readinsingle(nx,nq,nt,'lsamppole.dat',.true.)

        else

           reu = readinsingle(nx,nq,nt,'nlouamp.dat',.false.)
           imu = readinsingle(nx,nq,nt,'nlouamp.dat',.true.)
           reul = readinsingle(nx,nq,nt,'luamp.dat',.false.)
           imul = readinsingle(nx,nq,nt,'luamp.dat',.true.)

           red = readinsingle(nx,nq,nt,'nlodamp.dat',.false.)
           imd = readinsingle(nx,nq,nt,'nlodamp.dat',.true.)
           redl = readinsingle(nx,nq,nt,'ldamp.dat',.false.)
           imdl = readinsingle(nx,nq,nt,'ldamp.dat',.true.)

           res = readinsingle(nx,nq,nt,'nlosamp.dat',.false.)
           ims = readinsingle(nx,nq,nt,'nlosamp.dat',.true.)
           resl = readinsingle(nx,nq,nt,'lsamp.dat',.false.)
           imsl = readinsingle(nx,nq,nt,'lsamp.dat',.true.)

           reg = readinsingle(nx,nq,nt,'nlogamp.dat',.false.)
           img = readinsingle(nx,nq,nt,'nlogamp.dat',.true.)


           reup = readinsingle(nx,nq,nt,'nlouamppol.dat',.false.)
           imup = readinsingle(nx,nq,nt,'nlouamppol.dat',.true.)
           reupl = readinsingle(nx,nq,nt,'luamppol.dat',.false.)
           imupl = readinsingle(nx,nq,nt,'luamppol.dat',.true.)

           redp = readinsingle(nx,nq,nt,'nlodamppol.dat',.false.)
           imdp = readinsingle(nx,nq,nt,'nlodamppol.dat',.true.)
           redpl = readinsingle(nx,nq,nt,'ldamppol.dat',.false.)
           imdpl = readinsingle(nx,nq,nt,'ldamppol.dat',.true.)

           resp = readinsingle(nx,nq,nt,'nlosamppol.dat',.false.)
           imsp = readinsingle(nx,nq,nt,'nlosamppol.dat',.true.)
           respl = readinsingle(nx,nq,nt,'lsamppol.dat',.false.)
           imspl = readinsingle(nx,nq,nt,'lsamppol.dat',.true.)

           regp = readinsingle(nx,nq,nt,'nlogamppol.dat',.false.)
           imgp = readinsingle(nx,nq,nt,'nlogamppol.dat',.true.)


           reue = readinsingle(nx,nq,nt,'nlouampe.dat',.false.)
           imue = readinsingle(nx,nq,nt,'nlouampe.dat',.true.)
           reuel = readinsingle(nx,nq,nt,'luampe.dat',.false.)
           imuel = readinsingle(nx,nq,nt,'luampe.dat',.true.)

           rede = readinsingle(nx,nq,nt,'nlodampe.dat',.false.)
           imde = readinsingle(nx,nq,nt,'nlodampe.dat',.true.)
           redel = readinsingle(nx,nq,nt,'ldampe.dat',.false.)
           imdel = readinsingle(nx,nq,nt,'ldampe.dat',.true.)

           rese = readinsingle(nx,nq,nt,'nlosampe.dat',.false.)
           imse = readinsingle(nx,nq,nt,'nlosampe.dat',.true.)
           resel = readinsingle(nx,nq,nt,'lsampe.dat',.false.)
           imsel = readinsingle(nx,nq,nt,'lsampe.dat',.true.)

           rege = readinsingle(nx,nq,nt,'nlogampe.dat',.false.)
           imge = readinsingle(nx,nq,nt,'nlogampe.dat',.true.)


           reuep = readinsingle(nx,nq,nt,'nlouamppole.dat',.false.)
           imuep = readinsingle(nx,nq,nt,'nlouamppole.dat',.true.)
           reuepl = readinsingle(nx,nq,nt,'luamppole.dat',.false.)
           imuepl = readinsingle(nx,nq,nt,'luamppole.dat',.true.)

           redep = readinsingle(nx,nq,nt,'nlodamppole.dat',.false.)
           imdep = readinsingle(nx,nq,nt,'nlodamppole.dat',.true.)
           redepl = readinsingle(nx,nq,nt,'ldamppole.dat',.false.)
           imdepl = readinsingle(nx,nq,nt,'ldamppole.dat',.true.)

           resep = readinsingle(nx,nq,nt,'nlosamppole.dat',.false.)
           imsep = readinsingle(nx,nq,nt,'nlosamppole.dat',.true.)
           resepl = readinsingle(nx,nq,nt,'lsamppole.dat',.false.)
           imsepl = readinsingle(nx,nq,nt,'lsamppole.dat',.true.)

           regep = readinsingle(nx,nq,nt,'nlogamppole.dat',.false.)
           imgep = readinsingle(nx,nq,nt,'nlogamppole.dat',.true.)

        endif

           reut3 = readinsingle(nx,nq,nt,'luamptw3.dat',.false.)
           imut3 = readinsingle(nx,nq,nt,'luamptw3.dat',.true.)
           reut3d = readinsingle(nx,nq,nt,'luamptw3d.dat',.false.)
           imut3d = readinsingle(nx,nq,nt,'luamptw3d.dat',.true.)

           redt3 = readinsingle(nx,nq,nt,'ldamptw3.dat',.false.)
           imdt3 = readinsingle(nx,nq,nt,'ldamptw3.dat',.true.)
           redt3d = readinsingle(nx,nq,nt,'ldamptw3d.dat',.false.)
           imdt3d = readinsingle(nx,nq,nt,'ldamptw3d.dat',.true.)

           rest3 = readinsingle(nx,nq,nt,'lsamptw3.dat',.false.)
           imst3 = readinsingle(nx,nq,nt,'lsamptw3.dat',.true.)
           rest3d = readinsingle(nx,nq,nt,'lsamptw3d.dat',.false.)
           imst3d = readinsingle(nx,nq,nt,'lsamptw3d.dat',.true.)


           reupt3 = readinsingle(nx,nq,nt,'luamppoltw3.dat',.false.)
           imupt3 = readinsingle(nx,nq,nt,'luamppoltw3.dat',.true.)
           reupt3d = readinsingle(nx,nq,nt,'luamppoltw3d.dat',.false.)
           imupt3d = readinsingle(nx,nq,nt,'luamppoltw3d.dat',.true.)

           redpt3 = readinsingle(nx,nq,nt,'ldamppoltw3.dat',.false.)
           imdpt3 = readinsingle(nx,nq,nt,'ldamppoltw3.dat',.true.)
           redpt3d = readinsingle(nx,nq,nt,'ldamppoltw3d.dat',.false.)
           imdpt3d = readinsingle(nx,nq,nt,'ldamppoltw3d.dat',.true.)

           respt3 = readinsingle(nx,nq,nt,'lsamppoltw3.dat',.false.)
           imspt3 = readinsingle(nx,nq,nt,'lsamppoltw3.dat',.true.)
           respt3d = readinsingle(nx,nq,nt,'lsamppoltw3d.dat',.false.)
           imspt3d = readinsingle(nx,nq,nt,'lsamppoltw3d.dat',.true.)


           reuet3 = readinsingle(nx,nq,nt,'luampetw3.dat',.false.)
           imuet3 = readinsingle(nx,nq,nt,'luampetw3.dat',.true.)
           reuet3d = readinsingle(nx,nq,nt,'luampetw3d.dat',.false.)
           imuet3d = readinsingle(nx,nq,nt,'luampetw3d.dat',.true.)

           redet3 = readinsingle(nx,nq,nt,'ldampetw3.dat',.false.)
           imdet3 = readinsingle(nx,nq,nt,'ldampetw3.dat',.true.)
           redet3d = readinsingle(nx,nq,nt,'ldampetw3d.dat',.false.)
           imdet3d = readinsingle(nx,nq,nt,'ldampetw3d.dat',.true.)

           reset3 = readinsingle(nx,nq,nt,'lsampetw3.dat',.false.)
           imset3 = readinsingle(nx,nq,nt,'lsampetw3.dat',.true.)
           reset3d = readinsingle(nx,nq,nt,'lsampetw3d.dat',.false.)
           imset3d = readinsingle(nx,nq,nt,'lsampetw3d.dat',.true.)


           reuept3 = readinsingle(nx,nq,nt,'luamppoletw3.dat',.false.)
           imuept3 = readinsingle(nx,nq,nt,'luamppoletw3.dat',.true.)
           reuept3d = readinsingle(nx,nq,nt,'luamppoletw3d.dat',.false.)
           imuept3d = readinsingle(nx,nq,nt,'luamppoletw3d.dat',.true.)

           redept3 = readinsingle(nx,nq,nt,'ldamppoletw3.dat',.false.)
           imdept3 = readinsingle(nx,nq,nt,'ldamppoletw3.dat',.true.)
           redept3d = readinsingle(nx,nq,nt,'ldamppoletw3d.dat',.false.)
           imdept3d = readinsingle(nx,nq,nt,'ldamppoletw3d.dat',.true.)

           resept3 = readinsingle(nx,nq,nt,'lsamppoletw3.dat',.false.)
           imsept3 = readinsingle(nx,nq,nt,'lsamppoletw3.dat',.true.)
           resept3d = readinsingle(nx,nq,nt,'lsamppoletw3d.dat',.false.)
           imsept3d = readinsingle(nx,nq,nt,'lsamppoletw3d.dat',.true.)

      endif

      return
      end

C
C     subroutine to readin the t-dependence for various nuclei
C
****************************************************************
      subroutine readtdep(A1)
****************************************************************

      character*78 name
      integer marray
      real*8 A1,dum1,dum2
      parameter(marray = 40)

      real*8 tarray(marray),tvalarray(marray)

      common /tdependence/ tarray,tvalarray

         if (A1.eq.16D0) then

         name = 'tdep_o-16.dat'

         elseif(A1.eq.40D0) then

         name = 'tdep_ca-40.dat'

         elseif(A1.eq.110D0) then

         name = 'tdep_pd-110.dat'

         elseif(A1.eq.206D0) then

         name = 'tdep_pb-206.dat'

         endif

      open(unit =11,file=name,status='unknown')

       do i = 1,marray
          read(11,*) dum1,dum2
          tarray(i) = dum1
          tvalarray(i) = dum2
       enddo

      close(unit=11)

      return

      end


C
C    subroutine which computes the DVCS amplitudes h1 etc.
C
****************************************************************
      subroutine ampli(t,i,j)
****************************************************************

      integer i,j,n1,n2
      real*8 t
      integer iord
      parameter(mx=100,mq=100,mt=100)
      real*8 resu,resd,ress,resg,resu1,resd1,ress1,resg1,erru
      real*8 resue,resde,resse,resge,resu1e,resd1e,ress1e,resg1e
      real*8 resup,resdp,ressp,resgp,resu1p,resd1p,ress1p,resg1p
      real*8 resupe,resdpe,resspe,resgpe,resu1pe,resd1pe,ress1pe
     > ,resg1pe
      common/tempval/resu,resd,ress,resg,resu1,resd1,ress1,resg1,
     > resue,resde,resse,resge,resu1e,resd1e,ress1e,resg1e,
     > resup,resdp,ressp,resgp,resu1p,resd1p,ress1p,resg1p,
     > resupe,resdpe,resspe,resgpe,resu1pe,resd1pe,ress1pe
     > ,resg1pe
      real*8 f1u,f2u,f1d,f2d
      real*8 f1s,f2s
      real*8 g1,g1sea,gpi,exg
      common /form1/ f1u,f1d,f1s,exg,f2u,f2d,f2s,g1,g1sea,gpi
      integer nx,nq,nt
      common /counter/ nx,nq,nt
      integer icountxq,icountt
      common /counter1/ icountxq,icountt
      real*8 xx(mx),qq(mq),tt(mt)
      common /kins/ xx,qq,tt
      real*8 reu(mx,mq,mt),imu(mx,mq,mt),red(mx,mq,mt),imd(mx,mq,mt)
      real*8 res(mx,mq,mt),ims(mx,mq,mt),reg(mx,mq,mt),img(mx,mq,mt)
      real*8 reup(mx,mq,mt),imup(mx,mq,mt),redp(mx,mq,mt),imdp(mx,mq,mt)
      real*8 resp(mx,mq,mt),imsp(mx,mq,mt),regp(mx,mq,mt),imgp(mx,mq,mt)
      real*8 reue(mx,mq,mt),imue(mx,mq,mt),rede(mx,mq,mt),imde(mx,mq,mt)
      real*8 rese(mx,mq,mt),imse(mx,mq,mt),rege(mx,mq,mt),imge(mx,mq,mt)
      real*8 reuep(mx,mq,mt),imuep(mx,mq,mt),redep(mx,mq,mt),
     > imdep(mx,mq,mt)
      real*8 resep(mx,mq,mt),imsep(mx,mq,mt),regep(mx,mq,mt),
     > imgep(mx,mq,mt)
      real*8 y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1

      common /req/ y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1

      real*8 delt

      common /amps/ reu,imu,red,imd,res,ims,reg,img,reup,imup,
     > redp,imdp,resp,imsp,regp,imgp,reue,imue,rede,imde,rese,imse
     > ,rege,imge,reuep,imuep,redep,imdep,resep,imsep,regep,imgep

      real*8 reutw3,imutw3,redtw3,imdtw3,restw3,imstw3
     >  ,reuptw3,imuptw3,redptw3,imdptw3,resptw3,imsptw3
     >  ,reuetw3,imuetw3,redetw3,imdetw3,resetw3,imsetw3
     >  ,reueptw3,imueptw3,redeptw3,imdeptw3,reseptw3,imseptw3

      common /ampstw3/ reutw3,imutw3,redtw3,imdtw3,restw3,imstw3
     >  ,reuptw3,imuptw3,redptw3,imdptw3,resptw3,imsptw3
     >  ,reuetw3,imuetw3,redetw3,imdetw3,resetw3,imsetw3
     >  ,reueptw3,imueptw3,redeptw3,imdeptw3,reseptw3,imseptw3


      real*8 h1re,h1im,e1re,e1im,h1tre,h1tim,e1tre,e1tim,
     > h1retw3,h1imtw3,e1retw3,e1imtw3,h1tretw3,h1timtw3,e1tretw3
     > ,e1timtw3

      common /amps1/h1re,h1im,e1re,e1im,h1tre,h1tim,e1tre,e1tim

      common /amps2/h1retw3,h1imtw3,e1retw3,e1imtw3,h1tretw3
     > ,h1timtw3,e1tretw3,e1timtw3

      common /order/ iord

      integer isel_twist3
      common /isel_twist3/ isel_twist3

      external forms, amptw3

      interface
        function interpolate3D(nx,nq,nt,xx,qq,tt,x,q2,t,arr) result(amp)
            integer mx,mq,mt
            parameter(mx=100,mq=100,mt=100)
            integer nx,nq,nt
            real*8 xx(mx)
            real*8 qq(mq)
            real*8 tt(mt)
            real*8 x,q2,t
            real*8 arr(mx,mq,mt)

            real*8 amp
        end function
      end interface

      if (iord.eq.1) then

            sw = 0.0

         else

            sw = 1D0

       endif

      call forms(t)

      call amptw3(t,i,j)

C
C     switch back from x_A to x_bj since amplitudes are given in terms
C     of x_bj per nucleon.
C

      delt = del*A

C
C     twist-2 amplitudes
C

      if (spin1.eq.0D0) then

C
C     do 2-dim interpolation in the global grid to the values of x and Q^2
C

         if (icountxq.eq.1) then

            resu = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reu)
            resd = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,red)
            ress = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,res)
            resu1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imu)
            resd1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imd)
            ress1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,ims)

            if (iord.eq.2) then
             resg = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reg)
             resg1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,img)
            endif

         endif

       h1re  =A**2D0*(resu+resd+ress+sw*resg)*f1u
       h1im  =A**2D0*(resu1+resd1+ress1+sw*resg1)*f1u
       e1re  = 0D0
       e1im  = 0D0
       h1tre = 0D0
       h1tim = 0D0
       e1tre = 0D0
       e1tim = 0D0

         else ! spin1=1. (spin 1/2)
CLS ICI1

        if (icountxq.eq.1) then

            resu = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reu)
            resd = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,red)
            ress = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,res)
            resu1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imu)
            resd1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imd)
            ress1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,ims)

            if (iord.eq.2) then
             resg = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reg)
             resg1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,img)
            endif

         endif



      h1re=(resu*f1u+resd*f1d+ress*f1s
     > +sw*resg*exg)
      h1im=(resu1*f1u+resd1*f1d+ress1*f1s
     > +sw*resg1*exg)



      if (icountxq.eq.1) then

            resue = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reue)
            resde = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,rede)
            resse = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,rese)
            resu1e = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imue)
            resd1e = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imde)
            ress1e = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imse)

            if (iord.eq.2) then
             resge = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,rege)
             resg1e = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imge)
            endif

         endif

      e1re=(resue*f2u+resde*f2d+resse*f2s
     > +sw*resge*exg)
      e1im=(resu1e*f2u+resd1e*f2d+ress1e*f2s
     > +sw*resg1e*exg)

      if (icountxq.eq.1) then

            resup = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reup)
            resdp = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redp)
            ressp = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,resp)
            resu1p = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imup)
            resd1p = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdp)
            ress1p = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imsp)

            if (iord.eq.2) then
             resgp = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,regp)
             resg1p = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imgp)
            endif

         endif

      h1tre=(resup*g1+resdp*g1+ressp*g1sea
     > +sw*resgp*exg)
      h1tim=(resu1p*g1+resd1p*g1+ress1p*g1sea
     > +sw*resg1p*exg)

      if (icountxq.eq.1) then

            resupe = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reuep)
            resdpe = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redep)
            resspe = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,resep)
            resu1pe = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imuep)
            resd1pe = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdep)
            ress1pe = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imsep)

            if (iord.eq.2) then
             resgpe = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,regep)
             resg1pe = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imgep)
            endif

         endif

      e1tre=(resupe*gpi+resdpe*gpi+resspe*gpi
     > +sw*resgpe*exg)
      e1tim=(resu1pe*gpi+resd1pe*gpi+ress1pe*gpi
     > +sw*resg1pe*exg)

c      print *,"twist-2",t,delt,qsq,h1re,h1im,e1re,e1im,h1tre,h1tim
c     > ,e1tre,e1tim

            endif

C
C     twist-3 amplitudes
C

      if (spin1.eq.0D0) then

      h1retw3=A**2D0*(reutw3+redtw3+restw3)*f1u

      h1imtw3=A**2D0*(imutw3+imdtw3+imstw3)*f1u

      e1retw3=reuetw3*f2u+redetw3*f2d+resetw3*f2s

      e1imtw3=imuetw3*f2u+imdetw3*f2d+imsetw3*f2s

      h1tretw3=reuptw3*g1+redptw3*g1+resptw3*g1sea

      h1timtw3=imuptw3*g1+imdptw3*g1+imsptw3*g1sea

      e1tretw3=reueptw3*gpi+redeptw3*gpi+reseptw3*gpi

      e1timtw3=imueptw3*gpi+imdeptw3*gpi+imseptw3*gpi

      else

      h1retw3= reutw3*f1u+redtw3*f1d+restw3*f1s

      h1imtw3= imutw3*f1u+imdtw3*f1d+imstw3*f1s

      e1retw3=reuetw3*f2u+redetw3*f2d+resetw3*f2s

      e1imtw3=imuetw3*f2u+imdetw3*f2d+imsetw3*f2s

      h1tretw3=reuptw3*g1+redptw3*g1+resptw3*g1sea

      h1timtw3=imuptw3*g1+imdptw3*g1+imsptw3*g1sea

      e1tretw3=reueptw3*gpi+redeptw3*gpi+reseptw3*gpi

      e1timtw3=imueptw3*gpi+imdeptw3*gpi+imseptw3*gpi

c      print *,"twist-3",t,delt,qsq,h1retw3,h1imtw3,e1retw3,e1imtw3
c     > ,h1tretw3,h1timtw3,e1tretw3,e1timtw3
c      read (*,*)

      endif

      return

      end

C
C     subroutine to compute the effective twist-3 amplitudes
C
****************************************************************
      subroutine amptw3(t,i,j)
****************************************************************

      integer i,j,mx,mq
      parameter(mx=100,mq=100,mt=100)
      real*8 mp,delt
      real*8 y,tmin,testt,del,qsq,ep1,kfac,t,A,Z,spin1
      common /req/ y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1
      integer nx,nq,nt
      common /counter/ nx,nq,nt
      integer icountxq,icountt
      common /counter1/ icountxq,icountt
      real*8 xx(mx),qq(mq),tt(mt)
      common /kins/ xx,qq,tt
      real*8 resuh,resdh,ressh,resu1h,resd1h,ress1h,erru
      real*8 resu2h,resd2h,ress2h,resu3h,resd3h,ress3h
      real*8 resu4,resd4,ress4,resu5,resd5,ress5
      real*8 resu6,resd6,ress6,resu7,resd7,ress7
      real*8 resu8,resd8,ress8,resu9,resd9,ress9
      real*8 resuht,resdht,ressht,resu1ht,resd1ht,ress1ht
      real*8 resu2ht,resd2ht,ress2ht,resu3ht,resd3ht,ress3ht
      real*8 resue1,resde1,resse1,resu1e1,resd1e1,ress1e1
      real*8 resu2e1,resd2e1,ress2e1,resu3e1,resd3e1,ress3e1
      real*8 resuet,resdet,resset,resu1et,resd1et,ress1et
      real*8 resu2et,resd2et,ress2et,resu3et,resd3et,ress3et
      common/tempval2/resuh,resdh,ressh,resu1h,resd1h,ress1h,
     > resu2h,resd2h,ress2h,resu3h,resd3h,ress3h,
     > resu4,resd4,ress4,resu5,resd5,ress5,resu6,resd6,ress6,
     > resu7,resd7,ress7,resu8,resd8,ress8,resu9,resd9,ress9,
     > resuht,resdht,ressht,resu1ht,resd1ht,ress1ht,
     > resu2ht,resd2ht,ress2ht,resu3ht,resd3ht,ress3ht,
     > resue1,resde1,resse1,resu1e1,resd1e1,ress1e1,
     > resu2e1,resd2e1,ress2e1,resu3e1,resd3e1,ress3e1,
     > resuet,resdet,resset,resu1et,resd1et,ress1et,
     > resu2et,resd2et,ress2et,resu3et,resd3et,ress3et
      real*8 reut3(mx,mq,mt),imut3(mx,mq,mt),redt3(mx,mq,mt)
      real*8 rest3(mx,mq,mt),imst3(mx,mq,mt),imdt3(mx,mq,mt)
      real*8 reupt3(mx,mq,mt),imupt3(mx,mq,mt),redpt3(mx,mq,mt)
      real*8 respt3(mx,mq,mt),imspt3(mx,mq,mt),imdpt3(mx,mq,mt)
      real*8 reuet3(mx,mq,mt),imuet3(mx,mq,mt),redet3(mx,mq,mt)
      real*8 reset3(mx,mq,mt),imset3(mx,mq,mt),imdet3(mx,mq,mt)
      real*8 reuept3(mx,mq,mt),imuept3(mx,mq,mt),redept3(mx,mq,mt)
      real*8 resept3(mx,mq,mt),imsept3(mx,mq,mt),imdept3(mx,mq,mt)

      real*8 reut3d(mx,mq,mt),imut3d(mx,mq,mt),redt3d(mx,mq,mt)
      real*8 rest3d(mx,mq,mt),imst3d(mx,mq,mt),imdt3d(mx,mq,mt)
      real*8 reupt3d(mx,mq,mt),imupt3d(mx,mq,mt),redpt3d(mx,mq,mt)
      real*8 respt3d(mx,mq,mt),imspt3d(mx,mq,mt),imdpt3d(mx,mq,mt)
      real*8 reuet3d(mx,mq,mt),imuet3d(mx,mq,mt),redet3d(mx,mq,mt)
      real*8 reset3d(mx,mq,mt),imset3d(mx,mq,mt),imdet3d(mx,mq,mt)
      real*8 reuept3d(mx,mq,mt),imuept3d(mx,mq,mt),redept3d(mx,mq,mt)
      real*8 resept3d(mx,mq,mt),imsept3d(mx,mq,mt),imdept3d(mx,mq,mt)

      real*8 reutw3,imutw3,redtw3,imdtw3,restw3,imstw3
     >  ,reuptw3,imuptw3,redptw3,imdptw3,resptw3,imsptw3
     >  ,reuetw3,imuetw3,redetw3,imdetw3,resetw3,imsetw3
     >  ,reueptw3,imueptw3,redeptw3,imdeptw3,reseptw3,imseptw3

      real*8 reul(mx,mq,mt),imul(mx,mq,mt),redl(mx,mq,mt),imdl(mx,mq,mt)
      real*8 resl(mx,mq,mt),imsl(mx,mq,mt),regl(mx,mq,mt),imgl(mx,mq,mt)
      real*8 reupl(mx,mq,mt),imupl(mx,mq,mt),redpl(mx,mq,mt),
     > imdpl(mx,mq,mt)
      real*8 respl(mx,mq,mt),imspl(mx,mq,mt),regpl(mx,mq,mt),
     > imgpl(mx,mq,mt)
      real*8 reuel(mx,mq,mt),imuel(mx,mq,mt),redel(mx,mq,mt),
     > imdel(mx,mq,mt)
      real*8 resel(mx,mq,mt),imsel(mx,mq,mt),regel(mx,mq,mt),
     > imgel(mx,mq,mt)
      real*8 reuepl(mx,mq,mt),imuepl(mx,mq,mt),redepl(mx,mq,mt),
     > imdepl(mx,mq,mt)
      real*8 resepl(mx,mq,mt),imsepl(mx,mq,mt),regepl(mx,mq,mt),
     > imgepl(mx,mq,mt)

      common /ampsn/reul,imul,redl,imdl,resl,imsl,regl,imgl,reupl,imupl,
     > redpl,imdpl,respl,imspl,regpl,imgpl,reuel,imuel,redel,imdel,resel
     > ,imsel,regel,imgel,reuepl,imuepl,redepl,imdepl,resepl,imsepl
     > ,regepl,imgepl

      common /ampst3d/ reut3d,imut3d,redt3d,imdt3d,rest3d,imst3d
     >  ,reupt3d,imupt3d,redpt3d,imdpt3d,respt3d,imspt3d
     >  ,reuet3d,imuet3d,redet3d,imdet3d,reset3d,imset3d
     >  ,reuept3d,imuept3d,redept3d,imdept3d,resept3d,imsept3d

      common /ampst3/ reut3,imut3,redt3,imdt3,rest3,imst3
     >  ,reupt3,imupt3,redpt3,imdpt3,respt3,imspt3
     >  ,reuet3,imuet3,redet3,imdet3,reset3,imset3
     >  ,reuept3,imuept3,redept3,imdept3,resept3,imsept3

      common /ampstw3/ reutw3,imutw3,redtw3,imdtw3,restw3,imstw3
     >  ,reuptw3,imuptw3,redptw3,imdptw3,resptw3,imsptw3
     >  ,reuetw3,imuetw3,redetw3,imdetw3,resetw3,imsetw3
     >  ,reueptw3,imueptw3,redeptw3,imdeptw3,reseptw3,imseptw3

      integer isel_twist3
      common /isel_twist3/ isel_twist3

      if (isel_twist3.eq.0) then

      reutw3   = 0D0
      imutw3   = 0D0
      redtw3   = 0D0
      imdtw3   = 0D0
      restw3   = 0D0
      imstw3   = 0D0
      reuetw3  = 0D0
      redetw3  = 0D0
      resetw3  = 0D0
      imuetw3  = 0D0
      imdetw3  = 0D0
      imsetw3  = 0D0
      reuptw3  = 0D0
      redptw3  = 0D0
      resptw3  = 0D0
      imuptw3  = 0D0
      imdptw3  = 0D0
      imsptw3  = 0D0
      reueptw3 = 0D0
      redeptw3 = 0D0
      reseptw3 = 0D0
      imueptw3 = 0D0
      imdeptw3 = 0D0
      imseptw3 = 0D0

      else ! isel_twist3

C
C     effective tw-3 amplitude H
C

C
C     first change x_A (del) to x_bj per nucleon since the amplitudes are
C     given for x_bj per nulceon
C

      delt = del*A

      if (spin1.eq.0D0) then

         if (icountxq.eq.1) then

            resuh = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reul)
            resdh = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redl)
            ressh = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,resl)
            resu1h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imul)
            resd1h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdl)
            ress1h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imsl)

            resu2h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reut3d)
            resd2h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redt3d)
            ress2h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,rest3d)
            resu3h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imut3d)
            resd3h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdt3d)
            ress3h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imst3d)

         endif

      reutw3 = resuh  + resu2h
      imutw3 = resu1h + resu3h
      redtw3 = resdh  + resd2h
      imdtw3 = resd1h + resd3h
      restw3 = ressh  + ress2h
      imstw3 = ress1h + ress3h

      reuetw3 = 0D0
      redetw3 = 0D0
      resetw3 = 0D0
      imuetw3 = 0D0
      imdetw3 = 0D0
      imsetw3 = 0D0
      reuptw3 = 0D0
      redptw3 = 0D0
      resptw3 = 0D0
      imuptw3 = 0D0
      imdptw3 = 0D0
      imsptw3 = 0D0

      reueptw3 = 0D0
      redeptw3 = 0D0
      reseptw3 = 0D0
      imueptw3 = 0D0
      imdeptw3 = 0D0
      imseptw3 = 0D0

      else

         if (icountxq.eq.1) then

            resuh = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reul)
            resdh = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redl)
            ressh = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,resl)
            resu1h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imul)
            resd1h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdl)
            ress1h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imsl)

            resu2h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reut3d)
            resd2h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redt3d)
            ress2h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,rest3d)
            resu3h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imut3d)
            resd3h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdt3d)
            ress3h = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imst3d)

            resu4 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reut3)
            resd4 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redt3)
            ress4 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,rest3)
            resu5 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imut3)
            resd5 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdt3)
            ress5 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imst3)

            resu6 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reuet3)
            resd6 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redet3)
            ress6 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reset3)
            resu7 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imuet3)
            resd7 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdet3)
            ress7 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imset3)

            resu8 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reupt3)
            resd8 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redpt3)
            ress8 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,respt3)
            resu9 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imupt3)
            resd9 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdpt3)
            ress9 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imspt3)

         endif


       reutw3 = resuh + resu2h +  ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*((-t*del)/(4D0*mp**2*(2D0-del))
     > *(resu4 + resu6 - (2D0-del)*resu8/del))

       imutw3 = resu1h + resu3h +  ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*((-t*del)/(4D0*mp**2*(2D0-del))
     > *(resu5 + resu7 - (2D0-del)*resu9/del))

       redtw3 = resdh + resd2h +  ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*((-t*del)/(4D0*mp**2*(2D0-del))
     > *(resd4 + resd6 - (2D0-del)*resd8/del))

       imdtw3 = resd1h + resd3h +  ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*((-t*del)/(4D0*mp**2*(2D0-del))
     > *(resd5 + resd7 - (2D0-del)*resd9/del))

       restw3 = ressh + ress2h +  ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*((-t*del)/(4D0*mp**2*(2D0-del))
     > *(ress4 + ress6 - (2D0-del)*ress8/del))

       imstw3 = ress1h + ress3h +  ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*((-t*del)/(4D0*mp**2*(2D0-del))
     > *(ress5 + ress7 - (2D0-del)*ress9/del))

C
C     effective tw-3 amplitude E
C

       if (icountxq.eq.1) then

           resue1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reuel)
           resde1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redel)
           resse1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,resel)
           resu1e1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imuel)
           resd1e1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdel)
           ress1e1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imsel)


           resu2e1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reuet3d)
           resd2e1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redet3d)
           ress2e1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reset3d)
           resu3e1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imuet3d)
           resd3e1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdet3d)
           ress3e1 = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imset3d)

       endif

      reuetw3 = resue1 + resu2e1 + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(del/(2D0-del)*
     > (resu4 + resu6 - (2D0-del)*resu8/del))
c      print *,resue1,resu2e1, resu4,resu6,resu8
      imuetw3 = resu1e1 + resu3e1 + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(del/(2D0-del)*
     > (resu5 + resu7 - (2D0-del)*resu9/del))
c      print *,resu1e1,resu3e1,resu5,resu7,resu9
      redetw3 = resde1 + resd2e1 + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(del/(2D0-del)*
     > (resd4 + resd6 - (2D0-del)*resd8/del))

      imdetw3 = resd1e1 + resd3e1 + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(del/(2D0-del)*
     > (resd5 + resd7 - (2D0-del)*resd9/del))

      resetw3 = resse1 + ress2e1 + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(del/(2D0-del)*
     > (ress4 + ress6 - (2D0-del)*ress8/del))

      imsetw3 = ress1e1 + ress3e1 + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(del/(2D0-del)*
     > (ress5 + ress7 - (2D0-del)*ress9/del))



C
C    effective tw-3 amplitude \tilde H
C

      if (icountxq.eq.1) then

           resuht = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reupl)
           resdht = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redpl)
           ressht = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,respl)
           resu1ht = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imupl)
           resd1ht = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdpl)
           ress1ht = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imspl)

           resu2ht = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reupt3d)
           resd2ht = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redpt3d)
           ress2ht = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,respt3d)
           resu3ht = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imupt3d)
           resd3ht = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdpt3d)
           ress3ht = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imspt3d)

       endif

      reuptw3 = resuht + resu2ht + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(del*(1-t/(4D0*mp**2))*resu8
     > /(2D0-del) +(resu4 + resu6)*(t/(4D0*mp**2)))
c      print *,reuptw3,resuht,resu2ht,resu8,resu4,resu6
      imuptw3 = resu1ht + resu3ht + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(del*(1-t/(4D0*mp**2))*resu9
     > /(2D0-del) +(resu5 + resu7)*(t/(4D0*mp**2)))
c      print *,imuptw3
      redptw3 = resdht + resd2ht + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(del*(1-t/(4D0*mp**2))*resd8
     > /(2D0-del) +(resd4 + resd6)*(t/(4D0*mp**2)))
c      print *,redptw3
      imdptw3 = resd1ht + resd3ht + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(del*(1-t/(4D0*mp**2))*resd9
     > /(2D0-del) +(resd5 + resd7)*(t/(4D0*mp**2)))
c      print *,imdptw3
      resptw3 = ressht + ress2ht + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(del*(1-t/(4D0*mp**2))*ress8
     > /(2D0-del) +(ress4 + ress6)*(t/(4D0*mp**2)))
c      print* ,resptw3
      imsptw3 = ress1ht + ress3ht + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(del*(1-t/(4D0*mp**2))*ress9
     > /(2D0-del) +(ress5 + ress7)*(t/(4D0*mp**2)))
c      print* ,imsptw3


C
C     effective tw-3 amplitude \tilde E
C

      if (icountxq.eq.1) then

          resuet = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reuepl)
          resdet = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redepl)
          resset = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,resepl)
          resu1et = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imuepl)
          resd1et = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdepl)
          ress1et = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imsepl)

          resu2et = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,reuept3d)
          resd2et = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,redept3d)
          ress2et = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,resept3d)
          resu3et = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imuept3d)
          resd3et = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imdept3d)
          ress3et = interpolate3D(nx,nq,nt,xx,qq,tt,delt,qsq,t,imsept3d)
       endif

      reueptw3 = resuet + resu2et + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(
     > -(resu4 + resu6) + (2D0-del)*resu8/del)

      imueptw3 = resu1et + resu3et + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(
     > -(resu5 + resu7) + (2D0-del)*resu9/del)

      redeptw3 = resdet + resd2et + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(
     > -(resd4 + resd6) + (2D0-del)*resd8/del)

      imdeptw3 = resd1et + resd3et + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(
     > -(resd5 + resd7) + (2D0-del)*resd9/del)

      reseptw3 = resset + ress2et + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(
     > -(ress4 + ress6) + (2D0-del)*ress8/del)

      imseptw3 = ress1et + ress3et + ((2D0*mp**2*del)/
     > ((1D0-del)*(t-tmin)))*(
     > -(ress5 + ress7) + (2D0-del)*ress9/del)


      endif

      endif ! isel_twist3

      return

      end

C
C     functions to compute the various contributions of DVCS, interference and BH
C     to the UNP, LP and TP Xsections, separated by angular dependence. Also used to form
C     asymmetries
C
****************************************************************
      real*8 function tdvcsup(laml,lamp,jj,ii)
****************************************************************

C
C     definition block
C
      implicit none
      integer i,j
      integer ii,jj
      real*8 t,A,Z,spin1
      integer icountxq,icountt
      common /counter1/ icountxq,icountt
      real*8 h1re,h1im,h1tre,h1tim
      real*8 e1re,e1im,e1tre,e1tim
      real*8 h1retw3,h1imtw3,e1retw3
      real*8 e1imtw3,h1tretw3,h1timtw3
      real*8 e1tretw3,e1timtw3
      common /amps1/h1re,h1im,e1re,e1im,h1tre,h1tim,e1tre,e1tim
      common /amps2/h1retw3,h1imtw3,e1retw3,e1imtw3,h1tretw3
     > ,h1timtw3,e1tretw3,e1timtw3
      real*8 fact, rem,rem1
      real*8 fact2,fact3,rem2
      real*8 lamp, laml
      real*8 mp
      real*8 pi
      parameter(pi = 3.141592653589793)
      real*8 y,tmin,testt,del,qsq,ep1,kfac
      common /req/ y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1
      real*8 tbhup,tiup,tiup1,tilp,tilp1, tilp11,
     > titp,titp1,titp2,titp3,tiup11,titp4,titp5,
     > tbhlp,tbhtp,tdvcslp,tdvcstp,tdvcstp1,
     > tbhupc1,tbhupc2,tbhups1,tbhtps1,tbhlpc1,tbhtpc1,
     > tdvcsupcos,tdvcsupsin,tiupcos2,tiupsin2,
     > tdvcslpcos,tdvcslpsin,tilpsin2,tilpcos2,
     > tdvcstpcoscos,tdvcstpcossin,tdvcstpsincos,
     > tdvcstpsinsin,titpsinsin,titpsincos,titpcossin,
     > titpcoscos,titpcos2cos,titpcos2sin,titpsin2sin,
     > titpsin2cos

      real*8 f1p,f2p,f1n,f2n
      common /form2/ f1p,f2p,f1n,f2n

      integer isel_realpart
      common /isel_realpart/ isel_realpart

      external ampli

C
C     contributions for an unpolarized target and unpol/pol probe
C

C
C     DVCS^2 contributions
C


C
C     constant piece in the DVCS^2 term others are higher twist i.e. cos(phi)
C     etc. Eq. (43,66)
C

      i = jj
      j = ii

      t = testt

C
C     icountt determines whether values of x and Q^2 have changed (new interpolation)
C     or not.
C

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

C
C     spin-0 or spin-1/2
C

CLS ---> ICI2 modifs 
      if ((qsq.gt.2).and.(qsq.lt.80.)) then
ccc      write(92,*) qsq,del,(h1re/h1im)**2
      endif

CLS   HERE : real part=0

      if (isel_realpart.eq.1) then
       h1re  = 0.
       e1re  = 0.
       h1tre = 0.
       e1tre = 0.
      endif
      if (isel_realpart.eq.2) then
       h1re  = 0.
       e1re  = 0.
       h1tre = 0.
       e1tre = 0.
      endif

c      e1im  = 0.
c      e1tim = 0.
c      h1tim = 0.
CLS
   
      if (spin1.eq.0D0) then

      tdvcsup = 2.*(2.-2.*y + y**2)/(y**2*qsq)*(h1re**2 + h1im**2)

      else

       rem = 2.*(2.-2.*y + y**2)/(y**2*(2.-del)**2*qsq)*(
     >-del**2*2.*(h1re*e1re+h1im*e1im+h1tre*e1tre+h1tim*e1tim)
     >)

       rem1 = 2.*(2.-2.*y + y**2)/(y**2*(2.-del)**2*qsq)*
     > 4.*(1.-del)*(h1re**2 + h1im**2 + h1tre**2 + h1tim**2)

       rem2 = 2.*(2.-2.*y + y**2)/(y**2*(2.-del)**2*qsq)*(
     >-(del**2 + (2.-del)**2*t/(4.*mp**2))*(e1re**2 + e1im**2)
     >-del**2*t/(4.*mp**2)*(e1tre**2 + e1tim**2)
     >)

      tdvcsup = rem + rem1 + rem2

      endif

      return

C
C     cos(phi) piece in DVCS^2 (twist-3) Eq.(44,66)
C

      entry tdvcsupcos(laml,lamp,jj,ii)

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      tdvcsupcos = 8.*dsqrt(kfac)*(2.-y)/(2.-del)**2*
     > (4.*(1.-del)*(h1re*h1retw3+h1im*h1imtw3 + h1tre*h1tretw3
     > + h1tim*h1timtw3) -(del**2 + (2.-del)**2*t/(4.*mp**2))*
     > (e1re*e1retw3 + e1im*e1imtw3) -del**2*t/(4.*mp**2)*
     > (e1tre*e1tretw3 + e1tim*e1timtw3) -del**2*
     > (h1retw3*e1re + h1imtw3*e1im + e1retw3*h1re + e1imtw3*h1im
     > + h1tretw3*e1tre + h1timtw3*e1tim + e1tretw3*h1tre
     > + e1timtw3*h1tim))

      return

C
C     sin(phi) piece in DVCS^2 (twist-3) Eq.(44,66)
C

      entry tdvcsupsin(laml,lamp,jj,ii)

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      tdvcsupsin = -laml*8.*dsqrt(kfac)*y/(2.-del)**2*
     > (4.*(1.-del)*(h1imtw3*h1re - h1im*h1retw3 +
     > h1timtw3*h1tre - h1tim*h1tretw3) -(del**2 + (2.-del)**2*t
     > /(4.*mp**2))*(e1imtw3*e1re - e1im*e1retw3) -del**2*t/(4.*mp**2)
     > *(e1timtw3*e1tre - e1tim*e1tretw3) -del**2*(h1imtw3*e1re
     > - h1retw3*e1im + e1imtw3*h1re - h1im*e1retw3 + h1timtw3*e1tre
     > - h1tretw3*e1tim + e1timtw3*h1tre - h1tim*e1tretw3))

      return

C
C     BH^2 contributions
C


      entry tbhup(laml,lamp,jj,ii)
C
C     constant term in BH^2
C     see eq.(35) of BMK
C
      t = testt

      call forms(t)


      if (spin1.eq.0D0) then

      fact = -(Z)**2/(del**2*y**2*(1.+ep1)**2*t)

      tbhup = fact*(((2.-y)**2+y**2*(1+ep1**2)**2)*(4.*del**2*
     > mp**2/t + 4.*(1.-del) + (4.*del+ep1**2)*t/qsq) + 32.*
     > del**2*kfac*mp**2/t + 2.*ep1**2*(4.*(1.-y)*(3.*2.*ep1**2)
     > + y**2*(2.-ep1**4)) - 4.*del**2*(2.-y)**2*(2.+ep1**2)*
     > t/qsq)*f2p**2

      else

       fact = -1./(del**2*y**2*(1.+ep1)**2*t)

      tbhup =  fact*(8.*kfac*((2.+3.*ep1)*(qsq/t)*(f1p**2 - f2p**2
     > *t/(4.*mp**2)) + 2.*del**2*(f1p + f2p)**2) + (2.-y)**2*(
     > (2.+ep1)*((4.*del**2*mp**2/t)*(1+t/qsq)**2 + 4.*(1.-del)*(
     > 1. + del*t/qsq))*(f1p**2-t/(4.*mp**2)*f2p**2) +
     > 4.*del**2*(del + (1.-del+ep1/2.)
     > *(1-t/qsq)**2 - del*(1.-2.*del)*t**2/qsq**2)*(f1p+f2p)**2)
     > + 8.*(1+ep1)*(1.-y-ep1*y**2/4)*(2.*ep1*(1-t/(4.*mp**2))*(f1p**2
     > -t/(4.*mp**2)*f2p**2) - del**2*(1.-t/qsq)**2*(f1p+f2p)**2))

      endif

      return

      entry tbhupc1(laml,lamp,jj,ii)
C
C     cos(phi) term in BH^2 eq.(36)
C
      t = testt

      call forms(t)

      if (spin1.eq.0D0) then

      fact = -(Z)**2/(del**2*y**2*(1.+ep1)**2*t)

      tbhupc1 = fact*dsqrt(kfac)*8.*(2.-y)*(2.*del+ep1**2
     >          -4.*(del*mp)**2/t)*f2p**2

      else

         fact = -1./(del**2*y**2*(1.+ep1)**2*t)

      tbhupc1 =  fact*8.*dsqrt(kfac)*(((4.*del**2*mp**2/t)-2.*del-ep1)
     > *(f1p**2 - f2p**2*t/(4.*mp**2)) + 2.*del**2*(1.- (1.-2.*del)
     > *t/qsq)*(f1p+f2p)**2)*(2.-y)

      endif

      return

      entry tbhupc2(laml,lamp,jj,ii)
C
C     cos(2phi) term in BH^2 eq.(37)
C
      t = testt

      call forms(t)

      if (spin1.eq.0D0) then

      fact = -(Z)**2/(del**2*y**2*(1.+ep1)**2*t)

      tbhupc2 = -f2p**2*fact*kfac*32.*(del*mp)**2/t

      else

         fact = -1./(del**2*y**2*(1.+ep1)**2*t)

      tbhupc2 = fact*8.*kfac*((4.*mp**2/t)*(f1p**2-f2p**2*t/(4.*mp**2))
     > + 2.*(f1p+f2p)**2)*del**2

      endif

      return

C
C     interference contributions
C

      entry tiup11(laml,lamp,jj,ii)

C
C     constant factor in the unpol. interference term (twist-3), eq.(53,69,72)
C

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      if (spin1.eq.0D0) then

      fact = Z*8.*(2.-y)/(del*y**3*t)

      tiup11 = -fact*((2.-del)*(1.-y)-(1.-del)*(2.-y)**2
     > *(1.-tmin/t))*h1re*f2p*(t/qsq)

      else

        fact = -8.*(2.-y)/(del*y**3*t)

      tiup11 = fact*(((2.-y)**2/(1.-y))*kfac
     & * (f1p*h1re + del/(2.-del)*
     & (f1p + f2p) * h1tre - t/4./mp**2 * f2p * e1re)
     & + (t/qsq)*(1.-y)*(2.-del)*(f1p*h1re + del/(2.-del)*
     & (f1p + f2p) * h1tre - t/4./mp**2 * f2p * e1re
     & - del/(2.-del)*(f1p+f2p)*(del*(h1re+e1re)/(2.-del)
     & + h1tre)))

      endif

      return

      entry tiup(laml,lamp,jj,ii)

C
C     cos(phi) factor in interference term eq.(54,69)
C

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      if (spin1.eq.0D0) then

      fact = Z*8.*dsqrt(kfac)/(del*y**3*t)

      tiup = - (2.-2.*y+y**2) * fact * h1re * f2p

      else

         fact = -8.*dsqrt(kfac)/(del*y**3*t)

      tiup = (2.-2.*y+y**2) * fact *
     & (f1p*h1re + del/(2.-del)*
     & (f1p + f2p) * h1tre - t/(4.*mp**2) * f2p * e1re)

      endif

      return



      entry  tiup1(laml,lamp,jj,ii)
C
C     here sin(phi) term in the interference term eq.(54,69)
C
      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      if (spin1.eq.0D0) then

      fact = Z*8.*dsqrt(kfac)/(del*y**2*t)

      tiup1 = laml*(2.-y)*fact*f2p*h1im
c      print *,tiup1,h1im,fact,f2p,(2.-y)

      else

       fact = -8.*dsqrt(kfac)/(del*y**2*t)

      tiup1 = -laml*(2.-y)*fact*(f1p*h1im + del/(2.-del)
     & *(f1p + f2p) * h1tim - t/(4.*mp**2) * f2p * e1im)

      endif

      return

      entry tiupcos2(laml,lamp,jj,ii)

C
C     cos(2phi) factor in interference term eq.(55,69)
C


      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      if (spin1.eq.0D0) then

         fact = -16.*kfac/(del*y**3*t)

        tiupcos2 =  Z*(2.-y)*h1retw3*fact*f2p

         else

      fact = -16.*kfac/(del*y**3*t)

      tiupcos2 = (2.- y) * fact *
     & (f1p*h1retw3 + del/(2.-del)*
     & (f1p + f2p) * h1tretw3 - t/(4.*mp**2) * f2p * e1retw3)

       endif

      return

      entry  tiupsin2(laml,lamp,jj,ii)
C
C     here sin(2phi) term in the interference term eq.(55,69)
C
      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      if (spin1.eq.0D0) then

         fact = 16.*kfac/(del*y**2*t)

         tiupsin2 = -laml*Z*fact*h1imtw3*f2p

      else

      fact = -16.*kfac/(del*y**2*t)

      tiupsin2 = -laml*fact*(f1p*h1imtw3 + del/(2.-del)
     & *(f1p + f2p) * h1timtw3 - t/(4.*mp**2) * f2p * e1imtw3)

      endif

      return

C
C     unpol/pol probe and long. polarized target
C

C
C     DVCS^2
C

C
C     constant term in DVCS^2, eq.(46,67)
C

      entry tdvcslp(laml,lamp,jj,ii)

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      tdvcslp = 2.*laml * lamp * (2.-y)/
     & (y * (2.-del)**2 * qsq) *(
     & 4.*(1.-del)* 2.*(h1re * h1tre + h1im * h1tim)
     & - del**2 *2.*(h1re*e1tre +h1im*e1tim +h1tre*e1re +h1tim
     & *e1im) - del*(del**2/2. + (2.-del)*t/(4.*mp**2))*2.*
     & (e1re*e1tre +e1im*e1tim))

      return

C
C     cos(phi) term in DVCS^2, eq.(47,67)
C

      entry tdvcslpcos(laml,lamp,jj,ii)

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      tdvcslpcos = 8.*laml * lamp * dsqrt(kfac)/
     & (y * (2.-del)**2 * qsq) *(
     & 4.*(1.-del)*(h1retw3*h1tre + h1imtw3*h1tim + h1tretw3*h1re
     & + h1timtw3*h1im) - del**2*(h1retw3*e1tre + h1imtw3*e1tim
     & + e1tretw3*h1re + e1timtw3*h1im + h1tretw3*e1re
     & + h1timtw3*e1im + e1retw3*h1tre + e1imtw3*h1tim)
     & - del*(del**2/2. + (2.-del)*t/(4.*mp**2))*
     & (e1retw3*e1tre + e1imtw3*e1tim + e1tretw3*e1re
     & + e1timtw3*e1im))

      return

C
C     sin(phi) term in DVCS^2, eq.(47,67)
C

      entry tdvcslpsin(laml,lamp,jj,ii)

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      tdvcslpsin = -8.*laml * lamp * dsqrt(kfac)*(2.-y)/
     & (y**2 * (2.-del)**2 * qsq) *(
     & 4.*(1.-del)*(h1imtw3*h1tre - h1tim*h1retw3 +
     & h1timtw3*h1re - h1im*h1tretw3)
     & - del**2*(h1imtw3*e1tre - h1retw3*e1tim + e1timtw3*h1re
     & - h1im*e1tretw3 + h1timtw3*e1re - h1tretw3*e1im
     & + e1imtw3*h1tre - h1tim*e1retw3)
     & - del*(del**2/2. + (2.-del)*t/(4.*mp**2))*
     & (e1imtw3*e1tre - e1tim*e1retw3 + e1timtw3*e1re
     & - e1im*e1tretw3))

      return

C
C     BH^2
C

      entry  tbhlp(laml,lamp,jj,ii)
C
C     here const. term for long pol. eq.(38)
C
      t = testt
      call forms(t)
      fact = -8.*laml*lamp*dsqrt(1+ep1)*(2.-y)*(f1p+f2p)/
     & ((1-t/(4.*mp**2))*t*del*y*(1+ep1)**2)

      tbhlp = fact*(1./2.*(del/2.*(1-t/qsq) - t/(4.*mp**2))*(2.-del
     & - 2*(1.-del)**2*t/qsq + ep1*(1-t/qsq) - del*(1.-2.*del)
     & *(t/qsq)**2)*(f1p+f2p) + (1. - (1.-del)*t/qsq)*(del**2*mp**2
     & *(1.+t/qsq)**2/t + (1.-del)*(1.+del*t/qsq))
     & *(f1p+f2p*t/(4.*mp**2)))

      return


      entry  tbhlpc1(laml,lamp,jj,ii)
C
C     here cos(phi) long pol. eq.(39)
C
      t = testt

      call forms(t)

      fact = 8.*laml*lamp*dsqrt(1+ep1)*dsqrt(kfac)*(f1p+f2p)/
     & ((1-t/(4.*mp**2))*t*del*y*(1+ep1)**2)

      tbhlpc1 = fact*((t/(2*mp**2) - del*(1-t/qsq))*(1.-del*(1.-t/qsq))
     & *(f1p+f2p) + (1. + del - (3.-2*del)*(1+del*t/qsq)
     &-4.*del**2*mp**2*(1.+(t/qsq)**2)/t)*(f1p+f2p*t/(4.*mp**2)))

      return

C
C     interference term
C

      entry tilp(laml,lamp,jj,ii)
C
C     here sin(phi) term for long pol. eq.(58,70)
C
      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact = 8.*dsqrt(kfac)/(t*del*y**3)

      tilp = fact * lamp* (2.-2.*y + y**2)*
     & (del*(f1p + f2p)*(h1im + del/2.*e1im)/(2.-del)
     & + f1p*h1tim -  del*e1tim*(del*f1p/2.
     & + t*f2p/(4.*mp**2))/(2-del))


      return

      entry tilp1(laml,lamp,jj,ii)
C
C     here cos(phi) term for long pol. eq.(58,70)
C
      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact = - 8.*dsqrt(kfac)/(t*del)

      tilp1 =  (fact/y**2)*lamp * laml * (2.-y)*(
     &(del*(f1p+f2p))/(2.-del) *(h1re + e1re*del/2.)
     & +f1p*h1tre
     & - del*e1tre/(2.-del) *(del/2. *f1p
     & + t*f2p/(4.*mp**2)))

      return

      entry tilp11(laml,lamp,jj,ii)
C
C     here const in long pol. eq.(57,70,73)
C
      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact = - 8.*laml*lamp/(del*t*y**2)

      tilp11 = fact*(((2.-y)**2/(1.-y)+2)*kfac
     & *(del*(f1p+f2p)/(2.-del) *(h1re + e1re*del/2.)
     & +f1p*h1tre
     & - del*e1tre/(2.-del) *(del/2. *f1p
     & + t*f2p/(4.*mp**2))) + t/qsq*(1.-y)*(2.-del)
     & *(del*(f1p+f2p)/(2.-del)*(h1re + e1re*del/2.)
     & + f1p*h1tre - del*e1tre/(2.-del)*(del*f1p/2.
     & + t*f2p/(4.*mp**2))  - del*(f1p+f2p)/(2.-del)
     & *(h1re + del*e1re/2. + del*(h1tre + del*e1tre/2.)
     & /(2.-del))))

      return

C
C     cos(2phi) term, eq.(59,70)
C

      entry tilpcos2(laml,lamp,jj,ii)

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact = -16.*kfac/(t*del)

      tilpcos2 =  (fact/y**2)*lamp * laml *(
     &(del*(f1p+f2p))/(2.-del) *(h1retw3 + e1retw3*del/2.)
     & +f1p*h1tretw3
     & - del*e1tretw3/(2.-del) *(del/2. *f1p
     & + t*f2p/(4.*mp**2)))

      return

C
C     sin(2phi) term, eq.(59,70)
C

      entry tilpsin2(laml,lamp,jj,ii)

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact = 16.*kfac/(t*del)

      tilpsin2 = fact * lamp* ((2.-y)/y**3) *
     & (del/(2.-del)*(f1p + f2p)*(h1imtw3 + del/2.*e1imtw3)
     & + f1p*h1timtw3 -
     & del*e1timtw3/(2.-del) *(del*f1p/2. +
     & t*f2p/(4.*mp**2)))

      return

C
C     transversely polarized target and unpol/pol probe
C
C

C
C     DVCS^2
C

      entry tdvcstp(laml,lamp,jj,ii)
C
C     here cos(pphi) of const. DVCS^2 trans. pol. eq.(49,68)
C
      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact3 =  dsqrt(qsq*kfac)*(2.-y)/(y*qsq*mp*dsqrt(1.-y))

      tdvcstp = laml*fact3*
     & (4.*del*(h1re*e1tre + h1im*e1tim) -
     & 4.*(2.-del)*(h1tre*e1re + h1tim*e1im)+2.*del**2*
     & (e1re*e1tre+e1im*e1tim))/(2.-del)**2

      return

      entry tdvcstp1(laml,lamp,jj,ii)
C
C     here sin(pphi) of const. DVCS^2 trans. pol. eq.(49,68)
C
      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact3=-dsqrt(qsq*kfac)*(2.-2.*y+y**2)/(y**2*qsq*mp*dsqrt(1.-y))

      tdvcstp1 = 2.*fact3*(2.*(2.-del)*(h1im*e1re - h1re*e1im)
     & - 2.*del*(h1tim*e1tre - h1tre*e1tim))/(2.-del)**2

      return

C
C     cos(pphi) of cos(phi) term in DVCS^2, eq.(50,68)
C

      entry tdvcstpcoscos(laml,lamp,jj,ii)

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact3=4.*dsqrt(qsq)*kfac/(y*qsq*mp*dsqrt(1.-y))

      tdvcstpcoscos = laml*fact3*
     & (2.*del*(h1retw3*e1tre + h1imtw3*e1tim
     & + e1tretw3*h1re + e1timtw3*h1im) -
     & 2.*(2.-del)*(h1tretw3*e1re
     & + h1timtw3*e1im + e1retw3*h1tre + e1imtw3*h1tim)
     & +del**2*(e1retw3*e1tre + e1imtw3*e1tim + e1tretw3*e1re
     & + e1timtw3*e1im))/(2.-del)**2

      return

C
C    sin(pphi) of cos(phi) term in DVCS^2, eq.(50,68)
C

      entry tdvcstpsincos(laml,lamp,jj,ii)

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact3=-4.*dsqrt(qsq)*kfac*(2.-y)/(y**2*qsq*mp*dsqrt(1.-y))

      tdvcstpsincos = 2.*fact3*((2.-del)*(h1imtw3*e1re
     & - h1retw3*e1im - e1imtw3*h1re + h1im*e1retw3)
     & - del*(h1timtw3*e1tre
     & - h1tretw3*e1tim - e1timtw3*h1tre + h1tim*e1tretw3))
     & /(2.-del)**2

      return

C
C     cos(pphi) of sin(phi) term in DVCS^2, eq.(50,68)
C

      entry tdvcstpcossin(laml,lamp,jj,ii)

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact3=-4.*dsqrt(qsq)*kfac*(2.-y)/(y**2*qsq*mp*dsqrt(1.-y))

      tdvcstpcossin = fact3*
     & (2.*del*(h1imtw3*e1tre - h1retw3*e1tim + e1timtw3*h1re
     & - h1im*e1tretw3) -2.*(2.-del)*(h1timtw3*e1re - h1tretw3*e1im
     & + e1imtw3*h1tre - h1tim*e1retw3) + del**2*(e1imtw3*e1tre
     & - e1tim*e1retw3 + e1timtw3*e1re
     & - e1im*e1tretw3))/(2.-del)**2

      return

C
C     sin(pphi) of sin(phi) term in DVCS^2, eq.(50,68)
C

      entry tdvcstpsinsin(laml,lamp,jj,ii)

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact3=-4.*dsqrt(qsq)*kfac/(y*qsq*mp*dsqrt(1.-y))

      tdvcstpsinsin = 2.*laml*fact3*((2.-del)*(h1retw3*e1re
     & + h1imtw3*e1im - e1retw3*h1re - e1imtw3*h1im)
     & - del*(h1tretw3*e1tre + h1timtw3*e1tim - e1tretw3*h1tre
     & - e1timtw3*h1tim))/(2.-del)**2

      return

C
C     BH^2 term
C


      entry tbhtp(laml,lamp,jj,ii)
c
C     here const.*cos(pphi) BH^2 for trans. pol. eq.(40)
C
      t = testt

      call forms(t)

      fact3 = 8.*laml*(2.-y)*dsqrt(qsq*(1.+ep1)*kfac)
     & *(f1p+f2p)/(del**2*y*(1+ep1)**2*t*mp
     & *dsqrt(1.-y-ep1*y**2/4.))

      tbhtp = fact3*(del**3*mp**2*(1.-t/qsq)*(f1p+f2p)/qsq
     & + (1.-(1.-del)*t/qsq)*(del**2*mp**2*(1.-t/qsq)*f1p/t
     & + del*f2p/2.))

      return

      entry tbhtpc1(laml,lamp,jj,ii)
c
C     here cos(phi)*cos(pphi) BH^2 for trans. pol. eq.(41)
C
      t = testt

      call forms(t)

      fact3 = 16.*laml*dsqrt((1.-y-ep1*y**2/4.)*(1+ep1))*mp
     & *(f1p+f2p)/(del*y*(1.+ep1)**2*t*dsqrt(qsq))

      tbhtpc1 = fact3*(2*kfac*qsq/(t*(1.-y-ep1*y**2/4.))*(del
     & *(1-t/qsq)*f1p + t*f2p/(4.*mp**2)) + (1.+ep1)*del
     & *(1.-t/qsq)*(f1p+t*f2p/(4.*mp**2)))

      return

      entry tbhtps1(laml,lamp,jj,ii)
c
C     here sin(phi)*sin(pphi)  BH^2 for trans pol. eq.(42)
C
      t = testt

      call forms(t)

      fact3 = -16.*laml*dsqrt((1.-y-ep1*y**2/4.)*(1+ep1)**3)
     & *mp*(1.-t/qsq)*(f1p+f2p)/(y*(1.+ep1)**2*t*dsqrt(qsq))

      tbhtps1 = fact3*(f1p+t*f2p/(4.*mp**2))

      return

C
C     interference term
C


      entry titp(laml,lamp,jj,ii)

C
C     here cos(pphi) interference for trans. pol. (actually twist-3) eq.(61,71,74)
C

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact2 = - 8.*mp*dsqrt((1.-y)*kfac)/(del*y**2*t*dsqrt(qsq))

      titp = laml*(((2.-y)**2/(1.-y)+2)*((f1p+f2p)*(del**2*(h1re
     & + del*e1re/2.)/(2.-del) + del*t*e1re/(4.*mp**2)) - del**2
     & *f1p*(h1tre + del*e1tre/2.)/(2.-del) + t/(4.*mp**2)*(4.
     & *(1.-del)*f2p*h1tre/(2.-del) - (del*f1p + del**2*f2p
     & /(2.-del))*e1tre)) - t/mp**2*(f2p*h1tre - del*
     & (f1p+del*f2p/2.)*e1tre/(2.-del)))*fact2

      return

      entry titp1(laml,lamp,jj,ii)
C
C     here sin(pphi) interference for trans. pol. (actually twist-3) eq.(61,71,74)
C
      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact2= - 8.*mp*dsqrt((1.-y)*kfac)*(2.-y)
     & /(del*y**3*t*dsqrt(qsq))

      titp1 = fact2*(((2.-y)**2/(1.-y))*((del**2*f1p
     & - (1.-del)*t*f2p/mp**2)*h1im/(2.-del) + ((t/(4.*mp**2))
     & *((2.-del)*f1p + del**2*f2p/(2.-del)) + del**2*f1p
     & /(2.-del))*e1im - del**2*(f1p+f2p)*(h1tim + t*e1tim
     & /(4.*mp**2))) + t/mp**2*(f2p*h1im - f1p*e1im))

      return

       entry titp2(laml,lamp,jj,ii)
C
C     here  cos(phi)*cos(pphi) in interference for trans. pol. eq.(62,71,74)
C
       t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact2 = -8.*mp*dsqrt(1.-y)*(2.-y)/(del*y**2*t*dsqrt(qsq))

      titp2 =  laml*fact2*((f1p+f2p)*(del**2*(h1re
     & + del*e1re/2.)/(2.-del) + del*t*e1re/(4.*mp**2)) - del**2
     & *f1p*(h1tre + del*e1tre/2.)/(2.-del) + t/(4.*mp**2)*(4.
     & *(1.-del)*f2p*h1tre/(2.-del) - (del*f1p + del**2*f2p
     & /(2.-del))*e1tre))

      return


      entry titp3(laml,lamp,jj,ii)
C
C     here cos(phi)*sin(pphi) term eq.(62,71,74)
C
      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact2= 8.*mp*dsqrt(1.-y)*(2.-2.*y+y**2)
     & /(del*y**3*t*dsqrt(qsq))

      titp3 = fact2*((del**2*f1p
     & - (1.-del)*t*f2p/mp**2)*h1im/(2.-del) + ((t/(4.*mp**2))
     & *((2.-del)*f1p + del**2*f2p/(2.-del)) + del**2*f1p
     & /(2.-del))*e1im - del**2*(f1p+f2p)*(h1tim + t*e1tim
     & /(4.*mp**2)))

      return

      entry titp4(laml,lamp,jj,ii)
C
C     here  sin(phi)*cos(pphi) in interference for trans. pol. eq.(62,71,74)
C
       t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact2 = 8.*mp*dsqrt(1.-y)*(2.-2.*y+y**2)
     & /(del*y**3*t*dsqrt(qsq))

      titp4 =  fact2*((f1p+f2p)*(del**2*(h1im
     & + del*e1im/2.)/(2.-del) + del*t*e1im/(4.*mp**2)) - del**2
     & *f1p*(h1tim + del*e1tim/2.)/(2.-del) + t/(4.*mp**2)*(4.
     & *(1.-del)*f2p*h1tim/(2.-del) - (del*f1p + del**2*f2p
     & /(2.-del))*e1tim))

      return



      entry titp5(laml,lamp,jj,ii)
C
C     here sin(phi)*sin(pphi) term eq.(62,71,74)
C
      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact2= 8.*mp*dsqrt(1.-y)*(2.-y)
     & /(del*y**2*t*dsqrt(qsq))

      titp5 = laml*fact2*((del**2*f1p
     & - (1.-del)*t*f2p/mp**2)*h1re/(2.-del) + ((t/(4.*mp**2))
     & *((2.-del)*f1p + del**2*f2p/(2.-del)) + del**2*f1p
     & /(2.-del))*e1re - del**2*(f1p+f2p)*(h1tre + t*e1tre
     & /(4.*mp**2)))

      return

C
C     cos(pphi)*cos(2phi) term, eq.(63,71,74) start here on Monday!
C

      entry titpcos2cos(laml,lamp,jj,ii)

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact2 = -16.*mp*dsqrt((1.-y)*kfac)/(del*y**2*t*dsqrt(qsq))

      titpcos2cos =  laml*fact2*((f1p+f2p)*(del**2*(h1retw3
     & + del*e1retw3/2.)/(2.-del)+del*t*e1retw3/(4.*mp**2))-del**2
     & *f1p*(h1tretw3 + del*e1tretw3/2.)/(2.-del) + t/(4.*mp**2)*
     & (4.*(1.-del)*f2p*h1tretw3/(2.-del) - (del*f1p + del**2*f2p
     & /(2.-del))*e1tretw3))

      return

C
C     cos(2phi)*sin(pphi) term, eq.(63,71,74)
C

      entry titpcos2sin(laml,lamp,jj,ii)

      t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact2= 16.*mp*dsqrt((1.-y)*kfac)*(2.-y)
     & /(del*y**3*t*dsqrt(qsq))

      titpcos2sin = fact2*((del**2*f1p
     & - (1.-del)*t*f2p/mp**2)*h1imtw3/(2.-del) + ((t/(4.*mp**2))
     & *((2.-del)*f1p + del**2*f2p/(2.-del)) + del**2*f1p
     & /(2.-del))*e1imtw3 - del**2*(f1p+f2p)*(h1timtw3 + t*e1timtw3
     & /(4.*mp**2)))

      return

C
C     sin(2phi)*cos(pphi) term, eq.(63,71,74)
C

      entry titpsin2cos(laml,lamp,jj,ii)

       t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact2 = 16.*mp*dsqrt((1.-y)*kfac)*(2.-y)
     & /(del*y**3*t*dsqrt(qsq))

      titpsin2cos =  fact2*((f1p+f2p)*(del**2*(h1imtw3
     & + del*e1imtw3/2.)/(2.-del)+del*t*e1imtw3/(4.*mp**2))-del**2
     & *f1p*(h1timtw3 + del*e1timtw3/2.)/(2.-del) + t/(4.*mp**2)
     & *(4.*(1.-del)*f2p*h1timtw3/(2.-del) - (del*f1p + del**2*f2p
     & /(2.-del))*e1timtw3))

      return

C
C     sin(2phi)*sin(pphi) term, eq.(63,71,74)
C

      entry titpsin2sin(laml,lamp,jj,ii)

       t = testt

      i = jj
      j = ii

      if (icountt.eq.1) then

      call ampli(t,i,j)

      else

         call forms(t)

      endif

      fact2= 16.*mp*dsqrt((1.-y)*kfac)
     & /(del*y**2*t*dsqrt(qsq))

      titpsin2sin = laml*fact2*((del**2*f1p
     & - (1.-del)*t*f2p/mp**2)*h1retw3/(2.-del)+((t/(4.*mp**2))
     & *((2.-del)*f1p + del**2*f2p/(2.-del)) + del**2*f1p
     & /(2.-del))*e1retw3 - del**2*(f1p+f2p)*(h1tretw3
     & + t*e1tretw3/(4.*mp**2)))


      end

C
C     functions that store the different angular depnedences of the
C     various DVCS,
C     interference and BH terms
C

****************************************************************
      real*8 function corsin(phi)
      real*8 phi,cor1,arg
      external cor1
      arg = phi
      corsin = DSIN(arg)*cor1(phi)
      return
      end

      real*8 function corsin2phi(phi)
      real*8 phi,cor1,arg
      external cor1
      arg = phi
      corsin2phi = DSIN(2.*arg)*cor1(phi)
      return
      end

      real*8 function corsin2(phi)
      real*8 phi,cor1,arg
      external cor1
      arg = phi
      corsin2 = DSIN(arg)**2*cor1(phi)
      return
      end

      real*8 function corsin2phi2(phi)
      real*8 phi,cor1,arg
      external cor1
      arg = phi
      corsin2phi2 = DSIN(2.*arg)*DSIN(arg)*cor1(phi)
      return
      end

      real*8 function corcos(phi)
      real*8 phi,cor1,arg
      external cor1
      arg = phi
      corcos = DCOS(arg)*cor1(phi)
      return
      end

      real*8 function corcos2phi(phi)
      real*8 phi,cor1,arg
      external cor1
      arg = phi
      corcos2phi = DCOS(2.*arg)*cor1(phi)
      return
      end

      real*8 function corcos2(phi)
      real*8 phi,cor1,arg
      external cor1
      arg = phi
      corcos2 = DCOS(arg)**2*cor1(phi)
      return
      end

      real*8 function corcos2phi2(phi)
      real*8 phi,cor1,arg
      external cor1
      arg = phi
      corcos2phi2 = DCOS(2.*arg)*DCOS(arg)*cor1(phi)
      return
      end
****************************************************************

C
C   function that computes minus
C   the value of the product of the two BH propagators!
C   i.e. cor = - P_1 * P_2
C

****************************************************************
      real*8 function cor(phi)
      real*8 phi,brac,twophi,t,spin1
      real*8 y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp
      real*8 c1,c2,c3,c4
      common /req/ y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1
      t = testt
      twophi = 2.*phi
      c0 = kfac
      if(c0.lt.0.0) then
         write(6,*) 'c0 negative...for'
         write(6,*) 'phi,twophi,y,tmin,t,del,qsq'
         write(6,*)  phi,twophi,y,tmin,t,del,qsq
         write(6,*) 'STOP'
         stop
      else
      c1 = 2.*(c0)**0.5
      c2 = -t/qsq*(1.-del*(2.-y) + y*ep1/2.)
      c3 = -t/qsq*(1.-y*(1.+ep1)-del*(2.-y)+y*ep1/2.)
      c4 = y*ep1/2.
      brac = (1.-y+c1*dcos(phi)+c2-c4)*(1.+c1*dcos(phi)+c3+c4)
      cor = brac/(y**2*(1+ep1)**2)
      endif
      return
      end

C
C     function that computes minus the reciprocal of the product of the two
C     BH propagators!
C     cor1 = -1/ (P_1 P_2)
C     now for the unweigted cross section!
C
      real*8 function cor1(phi)
      real*8 phi,brac,twophi,t,spin1
      real*8 y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp
      real*8 c1,c2,c3,c4
      common /req/ y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1
      t = testt
      twophi = 2.*phi
      c0 = kfac
      if(c0.lt.0.0) then
         write(6,*) 'c0 negative...for'
         write(6,*) 'phi,twophi,y,tmin,t,del,qsq'
         write(6,*)  phi,twophi,y,tmin,t,del,qsq
         write(6,*) 'STOP'
         stop
      else
      c1 = 2.*(c0)**0.5
      c2 = -t/qsq*(1.-del*(2.-y) + y*ep1/2.)
      c3 = -t/qsq*(1.-y*(1.+ep1)-del*(2.-y)+y*ep1/2.)
      c4 = y*ep1/2.
      brac = (1.-y+c1*dcos(phi)+c2-c4)*(1.+c1*dcos(phi)+c3+c4)
      cor1 = y**2*(1+ep1)**2/brac
      endif
      return
      end
****************************************************************

C
C     function that sets the t dependence of the various amplitudes
C
****************************************************************
      subroutine forms (t)
****************************************************************

ccc   modif_forms

      implicit none

      integer marray
      parameter(marray = 40)
      real*8 t,spin1,dy,tin,temp
      real*8 kp,kn,ksea,mp,mn,mv,ma,mpi,ga3,mp1
      data kp,kn,ksea,mp1,mn,mv,ma,mpi,ga3/
     &       1.79,-1.91,-2.0,0.938272,0.939565,
     &       0.840,1.100,0.134977,1.267/
      real*8 gep,gmp,gmn,gen
      real*8 f1u,f2u,f1d,f2d
      real*8 ges,gms,f1s,f2s
      real*8 g1,g1sea,gpi,exg
      common /form1/ f1u,f1d,f1s,exg,f2u,f2d,f2s,g1,g1sea,gpi
      real*8 f1p,f2p,f1n,f2n
      common /form2/ f1p,f2p,f1n,f2n
      real*8 y,tmin,testt,del,qsq,ep1,kfac,A,Z
      common /req/ y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1
      real*8 tarray(marray),tvalarray(marray)
      common /tdependence/ tarray,tvalarray

      real*8  it_form
      real*8  bnew_q,bnew_g,xdel_s,bnew_slopeq,bnew_slopeg,q20
      common/form_factors/it_form,bnew_q,bnew_g,xdel_s,q20,
     +                            bnew_slopeq,bnew_slopeg
      real*8  bnew
      real*8  alpha_g,alpha_q

      gep = 1./(1.-t/mv**2)**2
      gmp = gep * (1. + kp)
      gmn = gep * kn
      gen = 0.0
      ges = 1./(1.-t/mv**2)**3
      gms = (1. + ksea) * ges
      f1n = -t*gmn/(4.*mn**2)/(1.-t/(4.*mn**2))
      f2n = gmn/(1.-t/(4.*mn**2))

      if (spin1.eq.0D0) then

         if (A.gt.1D0) then
            call ratint(tarray,tvalarray,marray,t,f2p,dy)
            f1p = f2p
            f1u = f2p
         endif

      else ! spin1.eq.0D0

       f1p = (gep - t*gmp/(4.*mp1**2))/(1.-t/(4.*mp1**2))
       f2p = (gmp - gep)/(1.-t/(4.*mp**2))

C
C     allow for different t-dependences for collider settings
C     Q^2 dependent slopes used for colliders
C

       if (it_form.eq.0d0) then
        bnew = bnew_q+bnew_slopeq*dlog(qsq/q20)
        bnew = bnew/2.d0
        f1u  = dexp(t*bnew)
       endif
       if (it_form.eq.2d0) then
        f1u  = (2.*f1p + f1n)/2.
       endif
       if (it_form.eq.4d0) then
        bnew = bnew_q+bnew_slopeq*dlog(qsq/q20)
        bnew = bnew/2.d0
        if (del.lt.xdel_s) f1u  = dexp(t*bnew)
        if (del.ge.xdel_s) f1u  = (2.*f1p + f1n)/2.
       endif

      endif ! spin1.eq.0D0

       alpha_q = 1.
       f1u = f1u * alpha_q
 
       if (it_form.eq.0d0) then
        f2u = f1u
        f1d = f1u
        f2d = f1u
        f1s = f1u
        f2s = f1u
       endif
       if (it_form.eq.2d0) then
        f2u = (2.*f2p + f2n)/2.
        f1d = f1p + 2.*f1n
        f2d = f2p + 2.*f2n
        f1s = (ges - t/(4.*mv**2)*gms)/(1-t/(4.*mv**2))
        f2s = (gms - ges)/(1.-t/(4.*mv**2))
       endif
       if (it_form.eq.4d0) then
        if (del.lt.xdel_s) then
        f2u = f1u
        f1d = f1u
        f2d = f1u
        f1s = f1u
        f2s = f1u
        endif
        if (del.ge.xdel_s) then
        f2u = (2.*f2p + f2n)/2.
        f1d = f1p + 2.*f1n
        f2d = f2p + 2.*f2n
        f1s = (ges - t/(4.*mv**2)*gms)/(1-t/(4.*mv**2))
        f2s = (gms - ges)/(1.-t/(4.*mv**2))
        endif
       endif

      g1    = 1./(1.-t/ma**2)**2
      g1sea = 1./(1.-t/ma**2)**3
      gpi   = 4 * ga3 * mp1**2/(mpi**2 -t)

      alpha_g = 1.

      if (it_form.eq.0d0) then
       bnew = bnew_g+bnew_slopeg*dlog(qsq/q20)
       bnew = bnew/2.d0
       exg  = alpha_g*dexp(t*bnew)
      endif
      if (it_form.eq.2d0) then
       exg = 1./(1.-t/ma**2)**3
      endif
      if (it_form.eq.4d0) then
       bnew = bnew_g+bnew_slopeg*dlog(qsq/q20)
       bnew = bnew/2.d0
       if (del.lt.xdel_s)  exg  = dexp(t*bnew)
       if (del.ge.xdel_s)  exg  = 1./(1.-t/ma**2)**3
      endif

      return
      end

C-----------------------------------------------
C
C     2-dim interpolation routine "Spline"
C
C-----------------------------------------------

****************************************************************
      subroutine polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
****************************************************************

      implicit none

      integer m,n,nmax,mmax
      parameter (nmax=100,mmax=100)
      real*8 dy,x1,x2,y,x1a(mmax),x2a(mmax),ya(mmax,mmax)
      real*8 b_n(nmax),c_n(nmax),d_n(nmax)
      real*8 b_m(mmax),c_m(mmax),d_m(mmax)
      integer j,k

      real*8 ymtmp(mmax), yntmp(nmax)

      real*8 seval_ls

      do 12 j=1,m

         do 11 k = 1,n

            yntmp(k) = ya(j,k)

 11      continue


            call spline_ls(n,x2a,yntmp,b_n,c_n,d_n)
	    ymtmp(j) = seval_ls(n,x2,x2a,yntmp,b_n,c_n,d_n)

 12      continue


            call spline_ls(m,x1a,ymtmp,b_m,c_m,d_m)
	    y  = seval_ls(m,x1,x1a,ymtmp,b_m,c_m,d_m)
	
	    dy = 0.d0

            return

      end

C-----------------------------------------------
C
C     2-dim interpolation routine "Andreas"
C
C-----------------------------------------------

****************************************************************
      subroutine polin2_and(x1a,x2a,ya,m,n,x1,x2,y,dy)
****************************************************************

      implicit none

      integer m,n,nmax,mmax
      parameter (nmax=100,mmax=100)
      real*8 dy,x1,x2,y,x1a(mmax),x2a(mmax),ya(mmax,mmax)
      integer j,k

      real*8 ymtmp(mmax), yntmp(nmax),temp(mmax),ybis
      real*8 b_n(nmax),c_n(nmax),d_n(nmax)
      real*8 b_m(mmax),c_m(mmax),d_m(mmax)
      real*8 seval_ls
      integer jsave,ksave

      ksave = 1
      jsave = 1

      do 12 j=1,m

         do 11 k = 1,n

            yntmp(k) = ya(j,k)

 11      continue

            call ratint(x2a,yntmp,n,x2,ymtmp(j),dy)
	
 12      continue

            call ratint(x1a,ymtmp,m,x1,y,dy)
	
            return

      end

C-----------------------------------------------
C
C     2-dim interpolation routine "EP/LS (1)"
C
C-----------------------------------------------

****************************************************************
      subroutine polin2_mod_1(x1a,x2a,ya,m,n,x1,x2,y,dy)
****************************************************************

C Modif 29.09.03, LS/EP
C because of unstabilities with initial routine

      implicit none

      integer m,n,nmax,mmax
      parameter (nmax=100,mmax=100)
      real*8 dy,x1,x2,y,x1a(mmax),x2a(mmax),ya(mmax,mmax)
      integer i,j,k
      real*8 ymtmp(nmax), yntmp(nmax)

      real*4 ymtmpr(58), yntmpr(40)
      real*4 x2ar(40),x1ar(58)

      real*4 r,rr,divdif

      do i=1,40
      x2ar(i) = sngl(x2a(i))
      enddo
      do i=1,58
      x1ar(i) = sngl(x1a(i))
      enddo

      do 12 j=1,58!m

         do 11 k = 1,40!n

            yntmp(k)  = ya(j,k)
            yntmpr(k) = sngl(ya(j,k))

 11      continue

            call ratint(x2a,yntmp,40,x2,ymtmp(j),dy)
            ymtmpr(j)=sngl(ymtmp(j))

c            call polint(yntmpr,x2ar,40,sngl(x2),r)
c            ymtmpr(j)=r

c            print *,'First ratint',sngl(ymtmp(j)),
c     +              '  polint',r

 12      continue

c            call ratint(x1a,ymtmp,58,x1,y,dy)

            rr =  divdif(ymtmpr,x1ar,58,sngl(x1),2)

c            print *,'Second ratint',sngl(y),'  polint',rr

            y=dble(rr)

            dy =0.d0

            return

      end

C-----------------------------------------------
C
C     2-dim interpolation routine "EP/LS (2)"
C
C-----------------------------------------------

****************************************************************
      subroutine polin2_mod_2(x1a,x2a,ya,m,n,x1,x2,y,dy)
****************************************************************

      implicit none

      integer m,n,nmax,mmax
      parameter (nmax=100,mmax=100)
      real*8 dy,x1,x2,y,x1a(mmax),x2a(mmax),ya(mmax,mmax)
      integer i,j,k
      real*8 ymtmp(nmax), yntmp(nmax)

      real*4 ymtmpr(58), yntmpr(40)
      real*4 x2ar(40),x1ar(58)

      real*4 r,rr,divdif

      do i=1,40
      x2ar(i) = sngl(x2a(i))
      enddo
      do i=1,58
      x1ar(i) = sngl(x1a(i))
      enddo

      do 12 j=1,58!m

         do 11 k = 1,40!n

            yntmp(k)  = ya(j,k)
            yntmpr(k) = sngl(ya(j,k))

 11      continue

cc            call ratint(x2a,yntmp,40,x2,ymtmp(j),dy)

            call polint(yntmpr,x2ar,40,sngl(x2),r)
            ymtmpr(j)=r

c            print *,'ratint',sngl(ymtmp(j)),
c     +              '  polint',r

 12      continue

cc            call ratint(x1a,ymtmp,58,x1,y,dy)

c            do i=1,58
c            print *,x1ar(i),ymtmpr(i)-sngl(ymtmp(i))
c            enddo
            rr =  divdif(ymtmpr,x1ar,58,sngl(x1),2)


            y=dble(rr)

cc            print *,x1,'ratint',sngl(y),'  polint',rr

            dy =0.d0

            return

      end

C-----------------------------------------------
C
C     2-dim interpolation routine "Linear"
C
C-----------------------------------------------

****************************************************************
      subroutine polin2_lin(x1a,x2a,ya,m,n,x1,x2,y,dy)
****************************************************************

      implicit none

      integer m,n,nmax,mmax
      integer na(2)
      parameter (nmax=100,mmax=100)
      real*8 dy,x1,x2,y,x1a(mmax),x2a(nmax),ya(mmax,nmax)
      integer j,k,ii

      logical fprint
      common/fordebug/fprint

      real xtmp(2)
      real ftmp(58,40)
      real atmp(58+40)
      real fint

      logical first
      data first /.true./

      NA(1) = m
      NA(2) = n

      XTMP(1) = sngl(x1)
      XTMP(2) = sngl(x2)

      do j=1,m
       do k=1,n
        ftmp(j,k) = sngl(ya(j,k))
       enddo
      enddo

      do j=1,m
       atmp(j) = sngl(x1a(j))
      enddo
      do k=1,n
       atmp(m+k) = sngl(x2a(k))
      enddo


      y = dble(FINT(2,XTMP,NA,ATMP,FTMP))

      if (fprint) then
       write(6,*) 'x q2 res = ',x1,x2,y
       write(6,*) 'xtmp ',xtmp(1),xtmp(2)
       write(6,*) 'valeurs de f quand x=xmin : '
       do ii=1,n
        write(6,*) ii,ftmp(1,ii)
       enddo
       write(6,*) 'valeurs de yk quand x=xmin : '
       do ii=1,n
        write(6,*) ii,ya(1,ii)
       enddo
      endif

cc      print *,'test2',sngl(y)

      dy = 0.0d0

      return
      end

C-----------------------------------------------
C
C     2-dim interpolation routine "Andreas OLD"
C
C-----------------------------------------------

****************************************************************
      SUBROUTINE POLINT_andreas_OLD (XA,YA,N,X,Y,DY)
****************************************************************

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (NMAX=100)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
c          IF(DEN.EQ.0.) PAUSE
c          THEN
c             print *,X
c             PAUSE
c             endif
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

      SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=100,TINY=1.E-25,MXX=1050)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
	
      NS=1
      HH=ABS(X-XA(1))
      DO 11 I=1,N
        H=ABS(X-XA(I))
        IF (H.EQ.0.)THEN
          Y=YA(I)
          DY=0.0
          RETURN
        ELSE IF (H.LT.HH) THEN
          NS=I
          HH=H
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)+TINY
C	print *, 99, YA(I), I
11    CONTINUE
      Y=YA(NS)
C	print *, 100, Y, NS
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          W=C(I+1)-D(I)
          H=XA(I+M)-X
          T=(XA(I)-X)*D(I)/H
          DD=T-C(I+1)
          IF(DD.EQ.0.) GOTO 12
          DD=W/DD
          D(I)=C(I+1)*DD
          C(I)=T*DD
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
C	print *, 101, Y
13    CONTINUE
      RETURN
	
      END

c ---------------------------------------------------------------------
      subroutine spline_ls(n,x,y,b,c,d)
c ---------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension x(100),y(100),b(100),c(100),d(100)

      nm1=n-1
      if(n.lt.2) return
      if(n.lt.3) go to 250
      d(1)=x(2)-x(1)
      c(2)=(y(2)-y(1))/d(1)
      do 210 i=2,nm1
        d(i)=x(i+1)-x(i)
        b(i)=2.0d0*(d(i-1)+d(i))
        c(i+1)=(y(i+1)-y(i))/d(i)
        c(i)=c(i+1)-c(i)
  210 continue
      b(1)=-d(1)
      b(n)=-d(n-1)
      c(1)=01.0d0
      c(n)=0.0d0
      if(n.eq.3) go to 215
      c(1)=c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
      c(n)=c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
      c(1)=c(1)*d(1)**2.0d0/(x(4)-x(1))
      c(n)=-c(n)*d(n-1)**2.0d0/(x(n)-x(n-3))
  215 continue
      do 220 i=2,n
        t=d(i-1)/b(i-1)
        b(i)=b(i)-t*d(i-1)
        c(i)=c(i)-t*c(i-1)
  220 continue
      c(n)=c(n)/b(n)
      do 230 ib=1,nm1
        i=n-ib
        c(i)=(c(i)-d(i)*c(i+1))/b(i)
  230 continue
      b(n)=(y(n)-y(nm1))/d(nm1)+d(nm1)*(c(nm1)+2.0d0*c(n))
      do 240 i=1,nm1
        b(i)=(y(i+1)-y(i))/d(i)-d(i)*(c(i+1)+2.0d0*c(i))
        d(i)=(c(i+1)-c(i))/d(i)
        c(i)=3.0d0*c(i)
  240 continue
      c(n)=3.0d0*c(n)
      d(n)=d(n-1)
      return
  250 continue
      b(1)=(y(2)-y(1))/(x(2)-x(1))
      c(1)=0.0d0
      d(1)=0.0d0
      b(2)=b(1)
      c(2)=0.0d0
      d(2)=0.0d0

      return
      end

c ---------------------------------------------------------------------
      real*8 function seval_ls(n,xx,x,y,b,c,d)
c ---------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension x(100),y(100),b(100),c(100),d(100)
      data i/1/

      if(i.ge.n) i=1
      if(xx.lt.x(i)) go to 310
      if(xx.le.x(i+1)) go to 330
  310 continue
      i=1
      j=n+1
  320 continue
      k=(i+j)/2
      if(xx.lt.x(k)) j=k
      if(xx.ge.x(k)) i=k
      if(j.gt.i+1) go to 320
  330 continue
      dx=xx-x(i)
      seval_ls=y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))

      return
      end


c ---------------------------------------------------------------------
      real*8 function dgausskeps(extern,a,b,eps)
c ---------------------------------------------------------------------
      implicit real*8(a-h,o-z)

      del=eps
      gauss2=0.0d0
      aa=a
    5 y=b-aa
      if(del-dabs(y))2,2,1
    2 bb=aa+y
      c1=0.5d0*(aa+bb)
      c2=c1-aa
      c3=c2
      n=1
      s8=0.d0
      s16=0.d0
    8 s8=s8+
     1   0.101228536290376d0*extern(c1+c2*0.960289856497536d0)+
     2   0.222381034453374d0*extern(c1+c2*0.796666477413627d0)+
     3   0.313706645877887d0*extern(c1+c2*0.525532409916329d0)+
     4   0.362683783378362d0*extern(c1+c2*0.183434642495650d0)
      s16=s16+
     1   0.027152459411754d0*extern(c1+c2*0.989400934991650d0)+
     2   0.062253523938648d0*extern(c1+c2*0.944575023073233d0)+
     3   0.095158511682493d0*extern(c1+c2*0.865631202387832d0)+
     4   0.124628971255534d0*extern(c1+c2*0.755404408355003d0)+
     5   0.149595988816577d0*extern(c1+c2*0.617876244402644d0)+
     6   0.169156519395003d0*extern(c1+c2*0.458016777657227d0)+
     7   0.182603415044924d0*extern(c1+c2*0.281603550779259d0)+
     8   0.189450610455069d0*extern(c1+c2*0.095012509837637d0)
      go to (11,12),n
   11 c2=-c2
      n=2
      go to 8
   12 s8=s8*c3
      s16=s16*c3
      if(dabs(s16-s8)-eps*(dabs(s16)+1.d0))3,3,4
    3 gauss2=gauss2+s16
      aa=bb
      go to 5
    4 y=0.5d0*y
      if(del-dabs(y))2,6,6
    6 write(6,*) ' too high accuracy required '
    1 dgausskeps=gauss2

      return
      end


***********************************************************
*
       function  sf2test(xin,q2in)
*
***********************************************************

***********************************************************
c  parametrization of f2-ep  from a fit to dis data
c  available for :  0.971085 >= x >=  0.000005  and
c
c                    50000.000 >= q2 >=       0.500
c*
c   for q2 or x outside limits, the closest
c  limit is assumed : f2(x>xmax,q2)=f2(xmax,q2)
c                     f2(x<xmin,q2)=f2(xmin,q2)
c                     f2(x,q2>q2max)=f2(x,q2max)
c                     f2(x,q2<q2min)=f2(x,q2min)
c*
c  input : x and q2
c  library : cernlib (mathlib) for tchebytchev poly. routine
c  comments, etc... to c. pascaud or f. zomer
***********************************************************
       real x,q2,q2min,q2max,xin,q2in
       double precision xmax,t,dchsum,coef(0:100),f(2)
       parameter(nq= 77,xmax= 0.971085)
       parameter(q2min=      0.500,q2max=  50000.000)
       double precision qbin( 77),xmin( 77)
       integer nqbin( 77)
       double precision param(0: 64,100)
       data   ilsave/0/

       data   qbin /      0.200
     &,     0.300,     0.400,     0.500,     0.600,     0.700,     0.800
     &,     0.900,     1.000,     1.100,     1.200,     1.300,     1.400
     &,     1.500,     1.600,     1.700,     1.800,     1.900,     1.960
     &,     2.000,     2.100,     2.200,     2.300,     2.400,     2.500
     &,     3.000,     3.500,     4.000,     4.500,     5.000,     6.000
     &,     7.000,     8.000,     9.000,    10.000,    15.000,    20.000
     &,    20.250,    25.000,    30.000,    35.000,    40.000,    45.000
     &,    50.000,    55.000,    60.000,    65.000,    70.000,    75.000
     &,    80.000,    85.000,    90.000,    95.000,   100.000,   150.000
     &,   200.000,   250.000,   300.000,   350.000,   400.000,   450.000
     &,   500.000,   550.000,   600.000,   650.000,   700.000,   750.000
     &,   800.000,   850.000,   900.000,   950.000,  1000.000,  1500.000
     &,  2000.000,  2500.000,  3000.000,  3500.000
     & /
       data   xmin /  0.0000050
     &, 0.0000050, 0.0000050, 0.0000050, 0.0000050, 0.0000050, 0.0000050
     &, 0.0000050, 0.0000050, 0.0000050, 0.0000050, 0.0000050, 0.0000050
     &, 0.0000050, 0.0000050, 0.0000050, 0.0000050, 0.0000050, 0.0000050
     &, 0.0000050, 0.0000050, 0.0000050, 0.0000050, 0.0000050, 0.0000050
     &, 0.0000050, 0.0000050, 0.0000050, 0.0000050, 0.0000050, 0.0000060
     &, 0.0000070, 0.0000080, 0.0000090, 0.0000100, 0.0000150, 0.0000200
     &, 0.0000202, 0.0000250, 0.0000300, 0.0000350, 0.0000400, 0.0000449
     &, 0.0000499, 0.0000549, 0.0000599, 0.0000649, 0.0000699, 0.0000749
     &, 0.0000799, 0.0000849, 0.0000899, 0.0000949, 0.0000999, 0.0001498
     &, 0.0001998, 0.0002497, 0.0002996, 0.0003496, 0.0003995, 0.0004495
     &, 0.0004994, 0.0005493, 0.0005993, 0.0006492, 0.0006992, 0.0007491
     &, 0.0007990, 0.0008490, 0.0008989, 0.0009489, 0.0009988, 0.0014982
     &, 0.0019976, 0.0024970, 0.0029964, 0.0034958
     & /


       if(ilsave.EQ.0) THEN
          ilsave=1
 280      format(6d12.5)
          open(unit=28,file='sf2test.dat',status='old')
ccc          open(unit=28,file='sf2test.dat',status='unknown')
          read(28,280)(param(k,  1),k=0, 64)
          read(28,280)(param(k,  2),k=0, 64)
          read(28,280)(param(k,  3),k=0, 64)
          read(28,280)(param(k,  4),k=0, 64)
          read(28,280)(param(k,  5),k=0, 64)
          read(28,280)(param(k,  6),k=0, 64)
          read(28,280)(param(k,  7),k=0, 64)
          read(28,280)(param(k,  8),k=0, 64)
          read(28,280)(param(k,  9),k=0, 64)
          read(28,280)(param(k, 10),k=0, 64)
          read(28,280)(param(k, 11),k=0, 64)
          read(28,280)(param(k, 12),k=0, 64)
          read(28,280)(param(k, 13),k=0, 64)
          read(28,280)(param(k, 14),k=0, 64)
          read(28,280)(param(k, 15),k=0, 64)
          read(28,280)(param(k, 16),k=0, 64)
          read(28,280)(param(k, 17),k=0, 64)
          read(28,280)(param(k, 18),k=0, 64)
          read(28,280)(param(k, 19),k=0, 64)
          read(28,280)(param(k, 20),k=0, 64)
          read(28,280)(param(k, 21),k=0, 64)
          read(28,280)(param(k, 22),k=0, 64)
          read(28,280)(param(k, 23),k=0, 64)
          read(28,280)(param(k, 24),k=0, 64)
          read(28,280)(param(k, 25),k=0, 64)
          read(28,280)(param(k, 26),k=0, 64)
          read(28,280)(param(k, 27),k=0, 64)
          read(28,280)(param(k, 28),k=0, 64)
          read(28,280)(param(k, 29),k=0, 64)
          read(28,280)(param(k, 30),k=0, 64)
          read(28,280)(param(k, 31),k=0, 64)
          read(28,280)(param(k, 32),k=0, 64)
          read(28,280)(param(k, 33),k=0, 64)
          read(28,280)(param(k, 34),k=0, 64)
          read(28,280)(param(k, 35),k=0, 64)
          read(28,280)(param(k, 36),k=0, 64)
          read(28,280)(param(k, 37),k=0, 64)
          read(28,280)(param(k, 38),k=0, 64)
          read(28,280)(param(k, 39),k=0, 64)
          read(28,280)(param(k, 40),k=0, 64)
          read(28,280)(param(k, 41),k=0, 64)
          read(28,280)(param(k, 42),k=0, 56)
          read(28,280)(param(k, 43),k=0, 42)
          read(28,280)(param(k, 44),k=0, 42)
          read(28,280)(param(k, 45),k=0, 42)
          read(28,280)(param(k, 46),k=0, 56)
          read(28,280)(param(k, 47),k=0, 42)
          read(28,280)(param(k, 48),k=0, 42)
          read(28,280)(param(k, 49),k=0, 42)
          read(28,280)(param(k, 50),k=0, 56)
          read(28,280)(param(k, 51),k=0, 56)
          read(28,280)(param(k, 52),k=0, 56)
          read(28,280)(param(k, 53),k=0, 30)
          read(28,280)(param(k, 54),k=0, 30)
          read(28,280)(param(k, 55),k=0, 29)
          read(28,280)(param(k, 56),k=0, 36)
          read(28,280)(param(k, 57),k=0, 43)
          read(28,280)(param(k, 58),k=0, 38)
          read(28,280)(param(k, 59),k=0, 40)
          read(28,280)(param(k, 60),k=0, 40)
          read(28,280)(param(k, 61),k=0, 40)
          read(28,280)(param(k, 62),k=0, 46)
          read(28,280)(param(k, 63),k=0, 40)
          read(28,280)(param(k, 64),k=0, 46)
          read(28,280)(param(k, 65),k=0, 31)
          read(28,280)(param(k, 66),k=0, 26)
          read(28,280)(param(k, 67),k=0, 24)
          read(28,280)(param(k, 68),k=0, 21)
          read(28,280)(param(k, 69),k=0, 28)
          read(28,280)(param(k, 70),k=0, 22)
          read(28,280)(param(k, 71),k=0, 51)
          read(28,280)(param(k, 72),k=0, 52)
          read(28,280)(param(k, 73),k=0, 34)
          read(28,280)(param(k, 74),k=0, 21)
          read(28,280)(param(k, 75),k=0, 63)
          read(28,280)(param(k, 76),k=0, 35)
          read(28,280)(param(k, 77),k=0, 33)
 281      format(25i3)
          read(28,281)(nqbin(k),k=1, 77)
          rewind(28)
          close(28)
C--                  precision is better than   0.229E-02
       ENDIF
       sf2test=0.
       x=xin
       q2=q2in
       if(x.GT.xmax)x=xmax*0.999999
       if(q2.GT.q2max)q2=q2max*0.999999
       if(q2.LE.q2min)q2=q2min*1.000001

       i = 0
       do il = 1,nq-1
          if(q2.GT.qbin(il).and.q2.LE.qbin(il+1)) i = il
       ENDDO
       if (i.eq.0) return

 1     do 3 k=0,1
          if(x.lt.xmin(i+k)) f(k+1)=-35.

          t=xmin(i+k)
          t=(2*log10(x)-log10(t)-log10(xmax))/(log10(xmax)-log10(t))

          do 2 j=0,nqbin(i+k)
 2        coef(j)=param(j,i+k)

          f(k+1) = dchsum(1,coef,nqbin(i+k),t)
          if(f(k+1).GT.6.305116e-16) THEN
             f(k+1)=log(f(k+1))
          ELSEIF(-f(k+1).GT.6.305116e-16)THEN
             f(k+1)=-70.-log(-f(k+1))
          ELSE
             f(k+1)=-35.
          ENDIF

 3     continue
       a=(f(2)-f(1))/log(qbin(i+1)/qbin(i))
       b=f(2)-a*log(qbin(i+1))
       a=a*log(q2)+b
       if(a.LT.-35.) THEN
          a=-a-70.
          a=-exp(a)
       ELSE
          a=exp(a)
       ENDIF
       sf2test=a-a*x
       return
       end
*
* DHCSUM
*
      FUNCTION DCHSUM(MODE,C,N,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(0:N)

      IF(MODE .EQ. 1) THEN
       H=X
       F=H
       V=1
      ELSE IF(MODE .EQ. 2) THEN
       H=2*X**2-1
       F=H
       V=1
      ELSE IF(MODE .EQ. 3) THEN
       H=2*X**2-1
       F=1
       V=X
      ELSE IF(MODE .EQ. 4) THEN
       H=2*X-1
       F=H
       V=1
      END IF

      ALFA=H+H
      B1=0
      B2=0
      DO 1 I = N,0,-1
      B0=C(I)+ALFA*B1-B2
      B2=B1
    1 B1=B0
      DCHSUM=V*(B0-F*B2)
      RETURN
      END
