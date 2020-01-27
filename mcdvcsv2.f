C
C     subroutine which calculates the DVCS cross section + other observables
C

       subroutine dvcsob(spin2,Z1,A1,nx1,nq1,count,x,q,s,t,
     > phi,pphi,theta,ichar,laml,lamp,iord,results) 

C 
C     Definition block 
C 

C
C     Definition block
C
cls      implicit none
         implicit real*8(a-h,o-z)

C
C     number of points in x and Q^2 in global grid
C

      integer nx,nq,nx1,nq1

C
C     external counter variable to avoid repeated readin of amplitudes. 
C     count has to be set to 1 whenever one changes from LO to NLO or vice versa!
C
      real*8 eps_ls

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
      integer mx,mq 
      parameter(mx=100,mq=100) 
      real*8 xx(mx), qq(mq)  

 
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
     &        1.79,-1.91,0.0,0.938272,0.939565, 
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

       common /kins/ xx,qq

C
C     kinematical variables
C
      common /req/ y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1 

C
C     # of x and Q^2
C      
      common /counter/ nx,nq

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

C
C     adjust del for nuclear number A 
C

            x = x/A1
            mp = mp1*A1
            A = A1
            Z = Z1
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
        
      call readin(iord,nx,nq) 
      
      if (A.gt.1D0) then
         call readtdep(A)
      endif

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

c            ymax = 2.*(dsqrt(1+ep1) - 1.)/ep1

            ymax = (qsq+testt)/(qsq+(A*del)*testt)

            y = qsq/((A*del)*(s - (mp/A)**2)) 
           
 
C 
C     y,t warning if y,t out of ex. range. Skip the calculation for these x, Q^2 value. 
C  
            if (y.gt.ymax.or.testt.gt.tmin) then 

               goto 200
                 
            endif 

C
C     calculation of the k-factor from BMK (Belitsky,Mueller, Kirchner)
C

            kfac = - (testt/qsq)*(1.-del)*(1.-y-y**2*ep1/4.)
     > *(1.-tmin/testt)*(dsqrt(1.+ep1) + (4.*del*(1.-del)+ep1)/
     > (4.*(1.-del))*(testt-tmin)/qsq)
 
C 
C     compute diff cross section unpol., pol. ,transversely pol. 
c 

C
C     First option: all the angles are integrated out and the results 
C     for the total cross section, DVCS, Interference, BH as well as
C     the results for the SSA and CA with and without twist-3 are stored in
C     the array results. Note here only unpol. target and probe except for the SSA.
C     Second option: most general option with all angles and polarizations free.

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
c        call D01AKF(cor1,phimin,phimax,epsabs,epsrel,result,abserr,w, 
c      &            lw,iw,liw,ifail) 
c       corint1 = result 
      
         corint1 = dgausskeps(cor1,phimin,phimax,eps_ls)

C
C     integral over phi from 0 to 2*pi of cos(phi)/P_1*P_2
C

c       call D01AKF(corcos,phimin,phimax,epsabs,epsrel,result,abserr,w, 
c      &            lw,iw,liw,ifail) 
c       corcosi = result

      corcosi   = dgausskeps(corcos,phimin,phimax,eps_ls)
     
C
C     integral over phi from 0 to 2*pi of cos(2phi)/P_1*P_2
C

c      call D01AKF(corcos2phi,phimin,phimax,epsabs,epsrel,result,abserr 
c     &            ,w,lw,iw,liw,ifail) 
c      corcosi2 = result
      
      corcosi2  = dgausskeps(corcos2phi,phimin,phimax,eps_ls)
      
C
C     integral over phi from 0 to 2*pi of sin(2phi)/P_1*P_2
C

c      call D01AKF(corsin2phi,phimin,phimax,epsabs,epsrel,result,abserr
c     &            ,w,lw,iw,liw,ifail)
c      corsin2p = result
      
      corsin2p  = dgausskeps(corsin2phi,phimin,phimax,eps_ls)

C
C     integral over phi from 0 to 2*pi of sin(phi)/P_1*P_2
C

c       call D01AKF(corsin,phimin,phimax,epsabs,epsrel,result,abserr,w, 
c      &            lw,iw,liw,ifail)
c       corsini = result
      corsini  = dgausskeps(corsin,phimin,phimax,eps_ls)
      
C
C     integral over phi from 0 to 2*pi of cos(phi)^2/P_1*P_2
C

c       call D01AKF(corcos2,phimin,phimax,epsabs,epsrel,result,abserr,w, 
c      &            lw,iw,liw,ifail) 
c       corcossq = result
      corcossq  = dgausskeps(corcos2,phimin,phimax,eps_ls)
      
C
C     integral over phi from 0 to 2*pi of sin(phi)^2/P_1*P_2
C

c      call D01AKF(corsin2,phimin,phimax,epsabs,epsrel,result,abserr,w, 
c     &            lw,iw,liw,ifail)
c      corsinsq = result
      corsinsq = dgausskeps(corsin2,phimin,phimax,eps_ls)
      
C
C     integral over phi from 0 to 2*pi of cos(2phi)*cos(phi)/P_1*P_2
C

c       call D01AKF(corcos2phi2,phimin,phimax,epsabs,epsrel,result,abserr 
c     &            ,w,lw,iw,liw,ifail) 
c      corcos2p2 = result
      corcos2p2 = dgausskeps(corcos2phi2,phimin,phimax,eps_ls)
ccc      print *,corcos2p2
      
C
C     integral over phi from 0 to 2*pi of sin(2phi)*sin(phi)/P_1*P_2
C

c      call D01AKF(corsin2phi2,phimin,phimax,epsabs,epsrel,result,abserr 
c     &            ,w,lw,iw,liw,ifail)
c      corsin2p2 = result
      
      corsin2p2 = dgausskeps(corsin2phi2,phimin,phimax,eps_ls)
 
C 
C     compute dM from phase space 
C 
 
      dM = al**3*del*y**2/(16.*pi**2*qsq**2* 
     >   (1. + 4.*del**2*mp**2/qsq)**0.5*A) 


C 
C     unpolarized diff. Xsection (unpol. probe/unpol. target)
C 

C
C     cos(phi) interference
C
      resa2 = ichar*2.*pi*corcosi*tiup(1D0,1D0,j,i)*dM

C
C    t does not change anymore 
C

      icountt = 2

C
C     const. interference
C
      resa4 = ichar*2.*pi*corint1*tiup11(1D0,1D0,j,i)*dM

C
C     cos^2(phi) term in CA
C

      resa20 = ichar*2.*pi*corcossq*tiup(1D0,1D0,j,i)*dM
cc      resa20lolo = ichar*2.*pi*corint1*tiup(1D0,1D0,j,i)*dM
      resa20lolo = ichar*2.*pi*1.*tiup(1D0,1D0,j,i)*dM
cc      print *,'cos^2(phi) term in CA',resa20,corcossq

C
C     cos(phi) term in CA
C

      resa21 = ichar*2.*pi*corcosi*tiup11(1D0,1D0,j,i)*dM
cc      resa21lolo = ichar*2.*pi*corint1*tiup11(1D0,1D0,j,i)*dM
      resa21lolo = ichar*2.*pi*1.*tiup11(1D0,1D0,j,i)*dM

C
C     cos(2phi) in interference
C

      resa23 = ichar*2.*pi*corcosi2*tiupcos2(1D0,1D0,j,i)*dM

C
C     cos(phi)*cos(2phi) term in CA
C
      resa24 = ichar*2.*pi*corcos2p2*tiupcos2(1D0,1D0,j,i)*dM
cc      resa24lolo = ichar*2.*pi*corint1*tiupcos2(1D0,1D0,j,i)*dM
      resa24lolo = ichar*2.*pi*1.*tiupcos2(1D0,1D0,j,i)*dM
cc      print *,'cos(phi)*cos(2phi) term in CA',resa24,corcos2p2

C
C     const. BH^2
C
cc      resa3 = 2.*pi*corint1*tbhup(1D0,1D0,j,i)*dM
      resa3 = 2.*pi*1.*tbhup(1D0,1D0,j,i)*dM

C
C     cos(phi) BH^2
C
      resa5 = 0.!2.*pi*corcosi*tbhupc1(1D0,1D0,j,i)*dM

C
C     cos(2phi) BH^2
C
      resa6 = 0.!2.*pi*corcosi2*tbhupc2(1D0,1D0,j,i)*dM

C
C     const. DVCS^2
C
      resa1 = 4.*pi**2*tdvcsup(1D0,1D0,j,i)*dM


C
C     long. pol. probe + unpol. target, diff Xsection
C      

C
C     sin(phi) interference
C
      resa11 = ichar*2.*pi*corsini*tiup1(laml,1D0,j,i)*dM

C
C     sin^2(phi) term in SSA
C

      resa22 = ichar*2.*pi*corsinsq*tiup1(laml,1D0,j,i)*dM

C
C     sin(2phi) term interference
C
      resa25 = ichar*2.*pi*corsin2p*tiupsin2(laml,1D0,j,i)*dM


C
C     sin(phi)*sins(2phi) term in SSA 
C
      resa26 = ichar*2.*pi*corsin2p2*tiupsin2(laml,1D0,j,i)*dM

C
C     collect final result
C

C
C     first DVCS^2
C

      dvcssq = resa1 
      results(2) = dvcssq*3.8937966E+5

C
C     interference
C

      intdvcsbh = resa4 + resa2 + resa23 + resa11 + resa25

      results(3) = intdvcsbh*3.8937966E+5

C
C     BH^2
C

      bhsq = resa3 + resa5  + resa6 

      results(4) = bhsq*3.8937966E+5
C
C     total Xsection w twist-3
C

      tot = dvcssq + intdvcsbh + bhsq
      results(1) = tot*3.8937966E+5

C
C     total Xsection wo twist-3
C

      results(13) = (resa1+resa4+resa2+resa3+resa5+resa6)
     > *3.8937966E+5 

C
C     full CA
C

ccc ICILOLO
ccc      results(5) = 2.*(resa20+(resa21+resa24))*3.8937966E+5
      results(5) = (0.*(resa21lolo)+
     +              1.*(resa20lolo)+
     +              0.*(resa24lolo))           *3.8937966E+5
      results(6) = (bhsq)*3.8937966E+5!(dvcssq+bhsq)*3.8937966E+5
C
C     CA without WW Twist-3
C

      results(7) = 2.*(resa20+resa21)*3.8937966E+5
      results(8) = (dvcssq+bhsq)*3.8937966E+5

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

      else

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

         endif


         return

 
 200     print *,"OUT OF RANGE!"," y=",y," ,ymax=",ymax," ,t=",testt,
     >   " ,tmin=",tmin 


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
      subroutine readin(nord,nx,nq) 
      implicit none 
      integer nord 
      integer nx,nq,mx,mq 
      parameter(mx=100,mq=100) 
      real*8 xx(mq),qq(mq) 
      real*8 dum1,dum2 
      integer i,j 
      real*8 x,q 
      real*8 y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1
      common /kins/ xx,qq 
      common /req/ y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1

      real*8 reu(mx,mq),imu(mx,mq),red(mx,mq),imd(mx,mq) 
      real*8 res(mx,mq),ims(mx,mq),reg(mx,mq),img(mx,mq) 
      real*8 reup(mx,mq),imup(mx,mq),redp(mx,mq),imdp(mx,mq) 
      real*8 resp(mx,mq),imsp(mx,mq),regp(mx,mq),imgp(mx,mq)
      real*8 reue(mx,mq),imue(mx,mq),rede(mx,mq),imde(mx,mq) 
      real*8 rese(mx,mq),imse(mx,mq),rege(mx,mq),imge(mx,mq) 
      real*8 reuep(mx,mq),imuep(mx,mq),redep(mx,mq),imdep(mx,mq) 
      real*8 resep(mx,mq),imsep(mx,mq),regep(mx,mq),imgep(mx,mq) 

      real*8 reut3(mx,mq),imut3(mx,mq),redt3(mx,mq) 
      real*8 rest3(mx,mq),imst3(mx,mq),imdt3(mx,mq) 
      real*8 reupt3(mx,mq),imupt3(mx,mq),redpt3(mx,mq) 
      real*8 respt3(mx,mq),imspt3(mx,mq),imdpt3(mx,mq)
      real*8 reuet3(mx,mq),imuet3(mx,mq),redet3(mx,mq)
      real*8 reset3(mx,mq),imset3(mx,mq),imdet3(mx,mq) 
      real*8 reuept3(mx,mq),imuept3(mx,mq),redept3(mx,mq) 
      real*8 resept3(mx,mq),imsept3(mx,mq),imdept3(mx,mq) 

      real*8 reut3d(mx,mq),imut3d(mx,mq),redt3d(mx,mq) 
      real*8 rest3d(mx,mq),imst3d(mx,mq),imdt3d(mx,mq) 
      real*8 reupt3d(mx,mq),imupt3d(mx,mq),redpt3d(mx,mq) 
      real*8 respt3d(mx,mq),imspt3d(mx,mq),imdpt3d(mx,mq)
      real*8 reuet3d(mx,mq),imuet3d(mx,mq),redet3d(mx,mq) 
      real*8 reset3d(mx,mq),imset3d(mx,mq),imdet3d(mx,mq) 
      real*8 reuept3d(mx,mq),imuept3d(mx,mq),redept3d(mx,mq) 
      real*8 resept3d(mx,mq),imsept3d(mx,mq),imdept3d(mx,mq) 
      
      real*8 reul(mx,mq),imul(mx,mq),redl(mx,mq),imdl(mx,mq) 
      real*8 resl(mx,mq),imsl(mx,mq),regl(mx,mq),imgl(mx,mq) 
      real*8 reupl(mx,mq),imupl(mx,mq),redpl(mx,mq),imdpl(mx,mq) 
      real*8 respl(mx,mq),imspl(mx,mq),regpl(mx,mq),imgpl(mx,mq)
      real*8 reuel(mx,mq),imuel(mx,mq),redel(mx,mq),imdel(mx,mq) 
      real*8 resel(mx,mq),imsel(mx,mq),regel(mx,mq),imgel(mx,mq) 
      real*8 reuepl(mx,mq),imuepl(mx,mq),redepl(mx,mq),imdepl(mx,mq) 
      real*8 resepl(mx,mq),imsepl(mx,mq),regepl(mx,mq),imgepl(mx,mq)

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
       
      
CMM   Read in data, depending whether LO (tord=1) or NLO (tord=2) 
C 
C     del and q loops 
C 

C
C     read in of twist-2
C

      if (spin1.eq.0D0) then


       if (nord.eq.1) then
 
      open(unit =11,file='luamp.dat',status='unknown') 
      do j = 1, nx 
            read(11,*) x 
            
            xx(j) = x 

         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            qq(i) = q 
            reu(j,i) = dum1 
            imu(j,i) = dum2 
            reul(j,i) = dum1 
            imul(j,i) = dum2 
         enddo 

      enddo 
      close(unit=11) 

      open(unit =11,file='ldamp.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            red(j,i) = dum1 
            imd(j,i) = dum2
            redl(j,i) = dum1 
            imdl(j,i) = dum2
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamp.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            res(j,i)= dum1 
            ims(j,i) = dum2
            resl(j,i)= dum1 
            imsl(j,i) = dum2
         enddo 
      enddo 
      close(unit=11)       

      else  

      open(unit =11,file='luamp.dat',status='unknown') 
      do j = 1, nx 
            read(11,*) x             
            xx(j) = x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            qq(i) = q 
            reul(j,i) = dum1 
            imul(j,i) = dum2 
            
         enddo 

      enddo 
      close(unit=11) 

      open(unit =11,file='ldamp.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            redl(j,i) = dum1 
            imdl(j,i) = dum2 
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamp.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            resl(j,i)= dum1 
            imsl(j,i) = dum2 
         enddo 
      enddo 
      close(unit=11) 
 
      open(unit =11,file='nlouamp.dat',status='unknown') 
      do j = 1, nx 
            read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            reu(j,i) = dum1 
            imu(j,i) = dum2 
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='nlodamp.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            red(j,i) = dum1 
            imd(j,i) = dum2 
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='nlosamp.dat',status='unknown') 
      do j = 1, nx 
            read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            res(j,i)= dum1 
            ims(j,i) = dum2 
         enddo 
      enddo 
      close(unit=11) 

      endif

C
C     readin of twist-3
C
      open(unit =11,file='luamptw3.dat',status='unknown') 
      do j = 1, nx 
            read(11,*) x 
            
            xx(j) = x 

         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            qq(i) = q 
            reut3(j,i) = dum1 
           imut3(j,i) = dum2 
           
        enddo 

      enddo 
      close(unit=11) 

      open(unit =11,file='ldamptw3.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           redt3(j,i) = dum1 
           imdt3(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamptw3.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           rest3(j,i)= dum1 
           imst3(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='luamptw3d.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
           
           xx(j) = x 

        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           qq(i) = q 
           reut3d(j,i) = dum1 
           imut3d(j,i) = dum2 
            
        enddo 

      enddo 
      close(unit=11) 

      open(unit =11,file='ldamptw3d.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           redt3d(j,i) = dum1 
           imdt3d(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamptw3d.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           rest3d(j,i)= dum1 
           imst3d(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      elseif(spin1.eq.1D0) then

         

         if (nord.eq.1) then
            
      open(unit =11,file='luamp.dat',status='unknown') 
      
      do j = 1, nx 
            read(11,*) x             
            xx(j) = x
           
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            qq(i) = q 
            reu(j,i) = dum1 
            imu(j,i) = dum2 
            reul(j,i) = dum1 
            imul(j,i) = dum2 
         enddo 

      enddo 
      close(unit=11) 

      open(unit =11,file='ldamp.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            red(j,i) = dum1 
            imd(j,i) = dum2
            redl(j,i) = dum1 
            imdl(j,i) = dum2
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamp.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            res(j,i)= dum1 
            ims(j,i) = dum2
            resl(j,i)= dum1 
            imsl(j,i) = dum2
         enddo 
      enddo 
      close(unit=11)       

      open(unit =11,file='luamppol.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            reup(j,i) = dum1 
            imup(j,i) = dum2
            reupl(j,i)= dum1 
            imupl(j,i) = dum2
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='ldamppol.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            redp(j,i) = dum1 
            imdp(j,i) = dum2 
            redpl(j,i) = dum1 
            imdpl(j,i) = dum2
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamppol.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           resp(j,i)= dum1 
           imsp(j,i) = dum2
           respl(j,i)= dum1 
           imspl(j,i) = dum2
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='luampe.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
           xx(j) = x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           qq(i) = q 
           reue(j,i) = dum1 
           imue(j,i) = dum2 
           reuel(j,i) = dum1 
           imuel(j,i) = dum2
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='ldampe.dat',status='unknown') 
      do j = 1, nx 
       read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           rede(j,i) = dum1 
           imde(j,i) = dum2 
           redel(j,i) = dum1 
           imdel(j,i) = dum2
        enddo 
      enddo 
      close(unit=11) 
      open(unit =11,file='lsampe.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           rese(j,i)= dum1 
           imse(j,i) = dum2
           resel(j,i)= dum1 
           imsel(j,i) = dum2
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='luamppole.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           reuep(j,i) = dum1 
           imuep(j,i) = dum2
           reuepl(j,i) = dum1 
           imuepl(j,i) = dum2
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='ldamppole.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           redep(j,i) = dum1 
           imdep(j,i) = dum2 
           redepl(j,i) = dum1 
           imdepl(j,i) = dum2
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamppole.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           resep(j,i)= dum1 
           imsep(j,i) = dum2 
           resepl(j,i)= dum1 
           imsepl(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      else

      open(unit =11,file='luamp.dat',status='unknown') 
      do j = 1, nx 
            read(11,*) x 
            
            xx(j) = x 

         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            qq(i) = q 
            reul(j,i) = dum1 
            imul(j,i) = dum2 
            
         enddo 

      enddo 
      close(unit=11) 

      open(unit =11,file='ldamp.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            redl(j,i) = dum1 
            imdl(j,i) = dum2 
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamp.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            resl(j,i)= dum1 
            imsl(j,i) = dum2 
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='luamppol.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            reupl(j,i) = dum1 
            imupl(j,i) = dum2 
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='ldamppol.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            redpl(j,i) = dum1 
            imdpl(j,i) = dum2 
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamppol.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           respl(j,i)= dum1 
           imspl(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='luampe.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
           xx(j) = x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           qq(i) = q 
           reuel(j,i) = dum1 
           imuel(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='ldampe.dat',status='unknown') 
      do j = 1, nx 
       read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           redel(j,i) = dum1 
           imdel(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 
      open(unit =11,file='lsampe.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           resel(j,i)= dum1 
           imsel(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='luamppole.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           reuepl(j,i) = dum1 
           imuepl(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='ldamppole.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           redepl(j,i) = dum1 
           imdepl(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamppole.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           resepl(j,i)= dum1 
           imsepl(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 
 
      open(unit =11,file='nlouamp.dat',status='unknown') 
      do j = 1, nx 
            read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            reu(j,i) = dum1 
            imu(j,i) = dum2 
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='nlodamp.dat',status='unknown') 
      do j = 1, nx 
         read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            red(j,i) = dum1 
            imd(j,i) = dum2 
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='nlosamp.dat',status='unknown') 
      do j = 1, nx 
            read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            res(j,i)= dum1 
            ims(j,i) = dum2 
         enddo 
      enddo 
      close(unit=11)

      open(unit =11,file='nlouamppol.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           reup(j,i) = dum1 
           imup(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='nlodamppol.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           redp(j,i) = dum1 
           imdp(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='nlosamppol.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           resp(j,i)= dum1 
           imsp(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 
 
      open(unit =11,file='nlogamp.dat',status='unknown') 
      do j = 1, nx 
            read(11,*) x 
         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            reg(j,i) = dum1 
            img(j,i) = dum2 
         enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='nlogamppol.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           regp(j,i) = dum1 
           imgp(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='nlouampe.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
           xx(j) = x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           qq(i) = q 
           reue(j,i) = dum1 
           imue(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='nlodampe.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           rede(j,i) = dum1 
           imde(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='nlosampe.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           rese(j,i)= dum1 
           imse(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 
 
      open(unit =11,file='nlouamppole.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           reuep(j,i) = dum1 
           imuep(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='nlodamppole.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           redep(j,i) = dum1 
           imdep(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 
      open(unit =11,file='nlosamppole.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           resep(j,i)= dum1 
           imsep(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 
 
      open(unit =11,file='nlogampe.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           rege(j,i) = dum1 
           imge(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='nlogamppole.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           regep(j,i) = dum1 
           imgep(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11)

      endif

C
C     readin of twist-3
C
      open(unit =11,file='luamptw3.dat',status='unknown') 
      do j = 1, nx 
            read(11,*) x 
            
            xx(j) = x 

         do i = 1,nq 
            read(11,*) q,dum1,dum2 
            qq(i) = q 
            reut3(j,i) = dum1 
           imut3(j,i) = dum2 
           
        enddo 

      enddo 
      close(unit=11) 

      open(unit =11,file='ldamptw3.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           redt3(j,i) = dum1 
           imdt3(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamptw3.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           rest3(j,i)= dum1 
           imst3(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='luamppoltw3.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           reupt3(j,i) = dum1 
           imupt3(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='ldamppoltw3.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           redpt3(j,i) = dum1 
           imdpt3(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamppoltw3.dat',status='unknown') 
      do j = 1, nx 
       read(11,*) x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          respt3(j,i)= dum1 
          imspt3(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='luampetw3.dat',status='unknown') 
      do j = 1, nx 
          read(11,*) x 
          xx(j) = x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          qq(i) = q 
          reuet3(j,i) = dum1 
          imuet3(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='ldampetw3.dat',status='unknown') 
      do j = 1, nx 
      read(11,*) x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          redet3(j,i) = dum1 
          imdet3(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 
      open(unit =11,file='lsampetw3.dat',status='unknown') 
      do j = 1, nx 
       read(11,*) x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          reset3(j,i)= dum1 
          imset3(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='luamppoletw3.dat',status='unknown') 
      do j = 1, nx 
       read(11,*) x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          reuept3(j,i) = dum1 
          imuept3(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='ldamppoletw3.dat',status='unknown') 
      do j = 1, nx 
       read(11,*) x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          redept3(j,i) = dum1 
          imdept3(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamppoletw3.dat',status='unknown') 
      do j = 1, nx 
       read(11,*) x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          resept3(j,i)= dum1 
          imsept3(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 
      
      open(unit =11,file='luamptw3d.dat',status='unknown') 
      do j = 1, nx 
           read(11,*) x 
           
           xx(j) = x 

        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           qq(i) = q 
           reut3d(j,i) = dum1 
           imut3d(j,i) = dum2 
            
        enddo 

      enddo 
      close(unit=11) 

      open(unit =11,file='ldamptw3d.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           redt3d(j,i) = dum1 
           imdt3d(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamptw3d.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           rest3d(j,i)= dum1 
           imst3d(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='luamppoltw3d.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           reupt3d(j,i) = dum1 
           imupt3d(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='ldamppoltw3d.dat',status='unknown') 
      do j = 1, nx 
        read(11,*) x 
        do i = 1,nq 
           read(11,*) q,dum1,dum2 
           redpt3d(j,i) = dum1 
           imdpt3d(j,i) = dum2 
        enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamppoltw3d.dat',status='unknown') 
      do j = 1, nx 
       read(11,*) x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          respt3d(j,i)= dum1 
          imspt3d(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='luampetw3d.dat',status='unknown') 
      do j = 1, nx 
          read(11,*) x 
          xx(j) = x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          qq(i) = q 
          reuet3d(j,i) = dum1 
          imuet3d(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='ldampetw3d.dat',status='unknown') 
      do j = 1, nx 
      read(11,*) x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          redet3d(j,i) = dum1 
          imdet3d(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 
      open(unit =11,file='lsampetw3d.dat',status='unknown') 
      do j = 1, nx 
       read(11,*) x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          reset3d(j,i)= dum1 
          imset3d(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='luamppoletw3d.dat',status='unknown') 
      do j = 1, nx 
       read(11,*) x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          reuept3d(j,i) = dum1 
          imuept3d(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='ldamppoletw3d.dat',status='unknown') 
      do j = 1, nx 
       read(11,*) x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          redept3d(j,i) = dum1 
          imdept3d(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 

      open(unit =11,file='lsamppoletw3d.dat',status='unknown') 
      do j = 1, nx 
       read(11,*) x 
       do i = 1,nq 
          read(11,*) q,dum1,dum2 
          resept3d(j,i)= dum1 
          imsept3d(j,i) = dum2 
       enddo 
      enddo 
      close(unit=11) 

      endif

      return 
      end 

C
C     subroutine to readin the t-dependence for various nuclei
C

      subroutine readtdep(A1)

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
 
      subroutine ampli(t,i,j) 
 
      integer i,j,n1,n2 
      real*8 t
      integer iord
      parameter(mx=100,mq=100)
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
      integer nx,nq
      common /counter/ nx,nq
      integer icountxq,icountt
      common /counter1/ icountxq,icountt
      real*8 xx(mx),qq(mq)
      common /kins/ xx,qq
      real*8 reu(mx,mq),imu(mx,mq),red(mx,mq),imd(mx,mq) 
      real*8 res(mx,mq),ims(mx,mq),reg(mx,mq),img(mx,mq) 
      real*8 reup(mx,mq),imup(mx,mq),redp(mx,mq),imdp(mx,mq) 
      real*8 resp(mx,mq),imsp(mx,mq),regp(mx,mq),imgp(mx,mq)
      real*8 reue(mx,mq),imue(mx,mq),rede(mx,mq),imde(mx,mq) 
      real*8 rese(mx,mq),imse(mx,mq),rege(mx,mq),imge(mx,mq) 
      real*8 reuep(mx,mq),imuep(mx,mq),redep(mx,mq),imdep(mx,mq) 
      real*8 resep(mx,mq),imsep(mx,mq),regep(mx,mq),imgep(mx,mq)
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

 
      external forms, amptw3

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
            
       call polin2(xx,qq,reu,nx,nq,delt,qsq,resu,erru)
       call polin2(xx,qq,red,nx,nq,delt,qsq,resd,erru)
       call polin2(xx,qq,res,nx,nq,delt,qsq,ress,erru)
       call polin2(xx,qq,imu,nx,nq,delt,qsq,resu1,erru)
       call polin2(xx,qq,imd,nx,nq,delt,qsq,resd1,erru)
       call polin2(xx,qq,ims,nx,nq,delt,qsq,ress1,erru)

            if (iord.eq.2) then
              call polin2(xx,qq,reg,nx,nq,delt,qsq,resg,erru)
              call polin2(xx,qq,img,nx,nq,delt,qsq,resg1,erru)
            endif 

         endif
       
       h1re=A**2D0*(resu+resd+ress+sw*resg)*f1u  
       h1im=A**2D0*(resu1+resd1+ress1+sw*resg1)*f1u
       e1re = 0D0
       e1im = 0D0
       h1tre = 0D0
       h1tim= 0D0
       e1tre= 0D0
       e1tim= 0D0

         else

        if (icountxq.eq.1) then

       call polin2(xx,qq,reu,nx,nq,delt,qsq,resu,erru)
c       print *,"1",icountxq,testt,delt,qsq,resu,erru
       call polin2(xx,qq,red,nx,nq,delt,qsq,resd,erru)
c       print *,"1",icountxq,testt,delt,qsq,resd,erru
       call polin2(xx,qq,res,nx,nq,delt,qsq,ress,erru)
c       print *,"1",icountxq,testt,delt,qsq,ress,erru
       call polin2(xx,qq,imu,nx,nq,delt,qsq,resu1,erru)
c       print *,"1",icountxq,testt,delt,qsq,resu1,erru
       call polin2(xx,qq,imd,nx,nq,delt,qsq,resd1,erru)
c       print *,"1",icountxq,testt,delt,qsq,resd1,erru
       call polin2(xx,qq,ims,nx,nq,delt,qsq,ress1,erru)
c       print *,"1",icountxq,testt,delt,qsq,ress1,erru

            if (iord.eq.2) then
              call polin2(xx,qq,reg,nx,nq,delt,qsq,resg,erru)
              call polin2(xx,qq,img,nx,nq,delt,qsq,resg1,erru)
            endif           

         endif    

         

      h1re=(resu*f1u+resd*f1d+ress*f1s
     > +sw*resg*exg)  
      h1im=(resu1*f1u+resd1*f1d+ress1*f1s
     > +sw*resg1*exg)  

      

      if (icountxq.eq.1) then

       call polin2(xx,qq,reue,nx,nq,delt,qsq,resue,erru)
c       print *,"2",icountxq,testt,delt,qsq,resue,erru
       call polin2(xx,qq,rede,nx,nq,delt,qsq,resde,erru)
c       print *,"2",icountxq,testt,delt,qsq,resde,erru
       call polin2(xx,qq,rese,nx,nq,delt,qsq,resse,erru)
c       print *,"2",icountxq,testt,delt,qsq,resse,erru
       call polin2(xx,qq,imue,nx,nq,delt,qsq,resu1e,erru)
c       print *,"2",icountxq,testt,delt,qsq,resu1e,erru
       call polin2(xx,qq,imde,nx,nq,delt,qsq,resd1e,erru)
c       print *,"2",icountxq,testt,delt,qsq,resd1e,erru
       call polin2(xx,qq,imse,nx,nq,delt,qsq,ress1e,erru)
c       print *,"2",icountxq,testt,delt,qsq,ress1e,erru

            if (iord.eq.2) then
              call polin2(xx,qq,rege,nx,nq,delt,qsq,resge,erru)
              call polin2(xx,qq,imge,nx,nq,delt,qsq,resg1e,erru)
            endif 

         endif 

      e1re=(resue*f2u+resde*f2d+resse*f2s
     > +sw*resge*exg)  
      e1im=(resu1e*f2u+resd1e*f2d+ress1e*f2s
     > +sw*resg1e*exg)   

      if (icountxq.eq.1) then

       call polin2(xx,qq,reup,nx,nq,delt,qsq,resup,erru)
c       print *,"3",icountxq,testt,delt,qsq,resup,erru
       call polin2(xx,qq,redp,nx,nq,delt,qsq,resdp,erru)
c       print *,"3",icountxq,testt,delt,qsq,resdp,erru
       call polin2(xx,qq,resp,nx,nq,delt,qsq,ressp,erru)
c       print *,"3",icountxq,testt,delt,qsq,ressp,erru
       call polin2(xx,qq,imup,nx,nq,delt,qsq,resu1p,erru)
c       print *,"3",icountxq,testt,delt,qsq,resu1p,erru
       call polin2(xx,qq,imdp,nx,nq,delt,qsq,resd1p,erru)
c       print *,"3",icountxq,testt,delt,qsq,resd1p,erru
       call polin2(xx,qq,imsp,nx,nq,delt,qsq,ress1p,erru)
c       print *,"3",icountxq,testt,delt,qsq,ress1p,erru
            if (iord.eq.2) then
              call polin2(xx,qq,regp,nx,nq,delt,qsq,resgp,erru)
              call polin2(xx,qq,imgp,nx,nq,delt,qsq,resg1p,erru)
            endif 

         endif    

      h1tre=(resup*g1+resdp*g1+ressp*g1sea
     > +sw*resgp*exg)  
      h1tim=(resu1p*g1+resd1p*g1+ress1p*g1sea
     > +sw*resg1p*exg) 
     
      if (icountxq.eq.1) then

       call polin2(xx,qq,reuep,nx,nq,delt,qsq,resupe,erru)
c        print *,"4",icountxq,testt,delt,qsq,resupe,erru
       call polin2(xx,qq,redep,nx,nq,delt,qsq,resdpe,erru)
c       print *,"4",icountxq,testt,delt,qsq,resdpe,erru
       call polin2(xx,qq,resep,nx,nq,delt,qsq,resspe,erru)
c       print *,"4",icountxq,testt,delt,qsq,resspe,erru
       call polin2(xx,qq,imuep,nx,nq,delt,qsq,resu1pe,erru)
c       print *,"4",icountxq,testt,delt,qsq,resu1pe,erru
       call polin2(xx,qq,imdep,nx,nq,delt,qsq,resd1pe,erru)
c       print *,"4",icountxq,testt,delt,qsq,resd1pe,erru
       call polin2(xx,qq,imsep,nx,nq,delt,qsq,ress1pe,erru)
c       print *,"4",icountxq,testt,delt,qsq,ress1pe,erru

            if (iord.eq.2) then
              call polin2(xx,qq,regep,nx,nq,delt,qsq,resgpe,erru)
              call polin2(xx,qq,imgep,nx,nq,delt,qsq,resg1pe,erru)
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
c      print *,reuetw3,f2u,redetw3,f2d,resetw3,f2s

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
      subroutine amptw3(t,i,j)


      integer i,j,mx,mq
      parameter(mx=100,mq=100)
      real*8 mp,delt
      real*8 y,tmin,testt,del,qsq,ep1,kfac,t,A,Z,spin1
      common /req/ y,tmin,testt,del,qsq,ep1,kfac,A,Z,mp,spin1
      integer nx,nq
      common /counter/nx,nq
      integer icountxq,icountt
      common /counter1/ icountxq,icountt
      real*8 xx(mx),qq(mq)
      common /kins/ xx,qq 
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
      real*8 reut3(mx,mq),imut3(mx,mq),redt3(mx,mq) 
      real*8 rest3(mx,mq),imst3(mx,mq),imdt3(mx,mq) 
      real*8 reupt3(mx,mq),imupt3(mx,mq),redpt3(mx,mq) 
      real*8 respt3(mx,mq),imspt3(mx,mq),imdpt3(mx,mq)
      real*8 reuet3(mx,mq),imuet3(mx,mq),redet3(mx,mq) 
      real*8 reset3(mx,mq),imset3(mx,mq),imdet3(mx,mq) 
      real*8 reuept3(mx,mq),imuept3(mx,mq),redept3(mx,mq) 
      real*8 resept3(mx,mq),imsept3(mx,mq),imdept3(mx,mq) 

      real*8 reut3d(mx,mq),imut3d(mx,mq),redt3d(mx,mq) 
      real*8 rest3d(mx,mq),imst3d(mx,mq),imdt3d(mx,mq) 
      real*8 reupt3d(mx,mq),imupt3d(mx,mq),redpt3d(mx,mq)
      real*8 respt3d(mx,mq),imspt3d(mx,mq),imdpt3d(mx,mq)
      real*8 reuet3d(mx,mq),imuet3d(mx,mq),redet3d(mx,mq) 
      real*8 reset3d(mx,mq),imset3d(mx,mq),imdet3d(mx,mq) 
      real*8 reuept3d(mx,mq),imuept3d(mx,mq),redept3d(mx,mq) 
      real*8 resept3d(mx,mq),imsept3d(mx,mq),imdept3d(mx,mq)

      real*8 reutw3,imutw3,redtw3,imdtw3,restw3,imstw3 
     >  ,reuptw3,imuptw3,redptw3,imdptw3,resptw3,imsptw3
     >  ,reuetw3,imuetw3,redetw3,imdetw3,resetw3,imsetw3 
     >  ,reueptw3,imueptw3,redeptw3,imdeptw3,reseptw3,imseptw3

      real*8 reul(mx,mq),imul(mx,mq),redl(mx,mq),imdl(mx,mq) 
      real*8 resl(mx,mq),imsl(mx,mq),regl(mx,mq),imgl(mx,mq) 
      real*8 reupl(mx,mq),imupl(mx,mq),redpl(mx,mq),imdpl(mx,mq) 
      real*8 respl(mx,mq),imspl(mx,mq),regpl(mx,mq),imgpl(mx,mq)
      real*8 reuel(mx,mq),imuel(mx,mq),redel(mx,mq),imdel(mx,mq) 
      real*8 resel(mx,mq),imsel(mx,mq),regel(mx,mq),imgel(mx,mq) 
      real*8 reuepl(mx,mq),imuepl(mx,mq),redepl(mx,mq),imdepl(mx,mq) 
      real*8 resepl(mx,mq),imsepl(mx,mq),regepl(mx,mq),imgepl(mx,mq)

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

       call polin2(xx,qq,reul,nx,nq,delt,qsq,resuh,erru)
       call polin2(xx,qq,redl,nx,nq,delt,qsq,resdh,erru)
       call polin2(xx,qq,resl,nx,nq,delt,qsq,ressh,erru)
       call polin2(xx,qq,imul,nx,nq,delt,qsq,resu1h,erru)
       call polin2(xx,qq,imdl,nx,nq,delt,qsq,resd1h,erru)
       call polin2(xx,qq,imsl,nx,nq,delt,qsq,ress1h,erru)

       call polin2(xx,qq,reut3d,nx,nq,delt,qsq,resu2h,erru)
       call polin2(xx,qq,redt3d,nx,nq,delt,qsq,resd2h,erru)
       call polin2(xx,qq,rest3d,nx,nq,delt,qsq,ress2h,erru)
       call polin2(xx,qq,imut3d,nx,nq,delt,qsq,resu3h,erru)
       call polin2(xx,qq,imdt3d,nx,nq,delt,qsq,resd3h,erru)
       call polin2(xx,qq,imst3d,nx,nq,delt,qsq,ress3h,erru)    

         endif
      
      reutw3 = resuh + resu2h
      imutw3 = resu1h + resu3h 
      redtw3 = resdh + resd2h 
      imdtw3 = resd1h + resd3h 
      restw3 = ressh + ress2h 
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

       call polin2(xx,qq,reul,nx,nq,delt,qsq,resuh,erru)
c      print *,"test",icountxq,testt,delt,qsq,resuh,erru
c      read(*,*)
       call polin2(xx,qq,redl,nx,nq,delt,qsq,resdh,erru)
c      print *,"test",icountxq,testt,delt,qsq,resdh,erru
       call polin2(xx,qq,resl,nx,nq,delt,qsq,ressh,erru)
c      print *,"test",icountxq,testt,delt,qsq,ressh,erru
       call polin2(xx,qq,imul,nx,nq,delt,qsq,resu1h,erru)
c      print *,"test",icountxq,testt,delt,qsq,resu1h,erru
       call polin2(xx,qq,imdl,nx,nq,delt,qsq,resd1h,erru)
c      print *,"test",icountxq,testt,delt,qsq,resd1h,erru
       call polin2(xx,qq,imsl,nx,nq,delt,qsq,ress1h,erru)
c      print *,"test",icountxq,testt,delt,qsq,ress1h,erru

       call polin2(xx,qq,reut3d,nx,nq,delt,qsq,resu2h,erru)
c      print *,"test1",icountxq,testt,delt,qsq,resu2h,erru
       call polin2(xx,qq,redt3d,nx,nq,delt,qsq,resd2h,erru)
c      print *,"test1",icountxq,testt,delt,qsq,resd2h,erru
       call polin2(xx,qq,rest3d,nx,nq,delt,qsq,ress2h,erru)
c      print *,"test1",icountxq,testt,delt,qsq,ress2h,erru
       call polin2(xx,qq,imut3d,nx,nq,delt,qsq,resu3h,erru)
c      print *,"test1",icountxq,testt,delt,qsq,resu3h,erru
       call polin2(xx,qq,imdt3d,nx,nq,delt,qsq,resd3h,erru)
c      print *,"test1",icountxq,testt,delt,qsq,resd3h,erru
       call polin2(xx,qq,imst3d,nx,nq,delt,qsq,ress3h,erru)
c      print *,"test1",icountxq,testt,delt,qsq,ress3h,erru

       call polin2(xx,qq,reut3,nx,nq,delt,qsq,resu4,erru)
c    print *,"test2",icountxq,testt,delt,qsq,resu4,erru
       call polin2(xx,qq,redt3,nx,nq,delt,qsq,resd4,erru)
c      print *,"test2",icountxq,testt,delt,qsq,resd4,erru
       call polin2(xx,qq,rest3,nx,nq,delt,qsq,ress4,erru)
c      print *,"test2",icountxq,testt,delt,qsq,ress4,erru
       call polin2(xx,qq,imut3,nx,nq,delt,qsq,resu5,erru)
c      print *,"test2",icountxq,testt,delt,qsq,resu5,erru
       call polin2(xx,qq,imdt3,nx,nq,delt,qsq,resd5,erru)
c      print *,"test2",icountxq,testt,delt,qsq,resd5,erru
       call polin2(xx,qq,imst3,nx,nq,delt,qsq,ress5,erru)
c      print *,"test2",icountxq,testt,delt,qsq,ress5,erru

       call polin2(xx,qq,reuet3,nx,nq,delt,qsq,resu6,erru)
c      print *,"test3",icountxq,testt,delt,qsq,resu6,erru
       call polin2(xx,qq,redet3,nx,nq,delt,qsq,resd6,erru)
c      print *,"test3",icountxq,testt,delt,qsq,resd6,erru
       call polin2(xx,qq,reset3,nx,nq,delt,qsq,ress6,erru)
c      print *,"test3",icountxq,testt,delt,qsq,ress6,erru
       call polin2(xx,qq,imuet3,nx,nq,delt,qsq,resu7,erru)
c      print *,"test3",icountxq,testt,delt,qsq,resu7,erru
       call polin2(xx,qq,imdet3,nx,nq,delt,qsq,resd7,erru)
c      print *,"test3",icountxq,testt,delt,qsq,resd7,erru
       call polin2(xx,qq,imset3,nx,nq,delt,qsq,ress7,erru)
c      print *,"test3",icountxq,testt,delt,qsq,ress7,erru

       call polin2(xx,qq,reupt3,nx,nq,delt,qsq,resu8,erru)
c      print *,"test4",icountxq,testt,delt,qsq,resu8,erru
       call polin2(xx,qq,redpt3,nx,nq,delt,qsq,resd8,erru)
c      print *,"test4",icountxq,testt,delt,qsq,resd8,erru
       call polin2(xx,qq,respt3,nx,nq,delt,qsq,ress8,erru)
c      print *,"test4",icountxq,testt,delt,qsq,ress8,erru
       call polin2(xx,qq,imupt3,nx,nq,delt,qsq,resu9,erru)
c      print *,"test4",icountxq,testt,delt,qsq,resu9,erru
       call polin2(xx,qq,imdpt3,nx,nq,delt,qsq,resd9,erru)
c      print *,"test4",icountxq,testt,delt,qsq,resd9,erru
       call polin2(xx,qq,imspt3,nx,nq,delt,qsq,ress9,erru)
c      print *,"test4",icountxq,testt,delt,qsq,ress9,erru

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

       call polin2(xx,qq,reuel,nx,nq,delt,qsq,resue1,erru)
c      print *,"test9",icountxq,testt,delt,qsq,resue1,erru
       call polin2(xx,qq,redel,nx,nq,delt,qsq,resde1,erru)
c      print *,"test9",icountxq,testt,delt,qsq,resde1,erru
       call polin2(xx,qq,resel,nx,nq,delt,qsq,resse1,erru)
c      print *,"test9",icountxq,testt,delt,qsq,resse1,erru
       call polin2(xx,qq,imuel,nx,nq,delt,qsq,resu1e1,erru)
c      print *,"test9",icountxq,testt,delt,qsq,resu1e1,erru
       call polin2(xx,qq,imdel,nx,nq,delt,qsq,resd1e1,erru)
c      print *,"test9",icountxq,testt,delt,qsq,resd1e1,erru
       call polin2(xx,qq,imsel,nx,nq,delt,qsq,ress1e1,erru)
c      print *,"test9",icountxq,testt,delt,qsq,ress1e1,erru

       call polin2(xx,qq,reuet3d,nx,nq,delt,qsq,resu2e1,erru)
c      print *,"test10",icountxq,testt,delt,qsq,resu2e1,erru
       call polin2(xx,qq,redet3d,nx,nq,delt,qsq,resd2e1,erru)
c      print *,"test10",icountxq,testt,delt,qsq,resd2e1,erru
       call polin2(xx,qq,reset3d,nx,nq,delt,qsq,ress2e1,erru)
c      print *,"test10",icountxq,testt,delt,qsq,ress2e1,erru
       call polin2(xx,qq,imuet3d,nx,nq,delt,qsq,resu3e1,erru)
c      print *,"test10",icountxq,testt,delt,qsq,resu3e1,erru
       call polin2(xx,qq,imdet3d,nx,nq,delt,qsq,resd3e1,erru)
c      print *,"test10",icountxq,testt,delt,qsq,resd3e1,erru
       call polin2(xx,qq,imset3d,nx,nq,delt,qsq,ress3e1,erru)
c      print *,"test10",icountxq,testt,delt,qsq,ress3e1,erru

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

       call polin2(xx,qq,reupl,nx,nq,delt,qsq,resuht,erru)
c      print *,"test5",icountxq,testt,delt,qsq,resuht,erru
c      print *,reupl
       call polin2(xx,qq,redpl,nx,nq,delt,qsq,resdht,erru)
c      print *,"test5",icountxq,testt,delt,qsq,resdht,erru
       call polin2(xx,qq,respl,nx,nq,delt,qsq,ressht,erru)
c      print *,"test5",icountxq,testt,delt,qsq,ressht,erru
       call polin2(xx,qq,imupl,nx,nq,delt,qsq,resu1ht,erru)
c      print *,"test5",icountxq,testt,delt,qsq,resu1ht,erru
       call polin2(xx,qq,imdpl,nx,nq,delt,qsq,resd1ht,erru)
c      print *,"test5",icountxq,testt,delt,qsq,resd1ht,erru
       call polin2(xx,qq,imspl,nx,nq,delt,qsq,ress1ht,erru)
c      print *,"test5",icountxq,testt,delt,qsq,ress1ht,erru

       call polin2(xx,qq,reupt3d,nx,nq,delt,qsq,resu2ht,erru)
c      print *,"test6",icountxq,testt,delt,qsq,resu2ht,erru
       call polin2(xx,qq,redpt3d,nx,nq,delt,qsq,resd2ht,erru)
c      print *,"test6",icountxq,testt,delt,qsq,resd2ht,erru
       call polin2(xx,qq,respt3d,nx,nq,delt,qsq,ress2ht,erru)
c      print *,"test6",icountxq,testt,delt,qsq,ress2ht,erru
       call polin2(xx,qq,imupt3d,nx,nq,delt,qsq,resu3ht,erru)
c      print *,"test6",icountxq,testt,delt,qsq,resu3ht,erru
       call polin2(xx,qq,imdpt3d,nx,nq,delt,qsq,resd3ht,erru)
c      print *,"test6",icountxq,testt,delt,qsq,resd3ht,erru
       call polin2(xx,qq,imspt3d,nx,nq,delt,qsq,ress3ht,erru)
c      print *,"test6",icountxq,testt,delt,qsq,ress3ht,erru

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

       call polin2(xx,qq,reuepl,nx,nq,delt,qsq,resuet,erru)
c      print *,"test7",icountxq,testt,delt,qsq,resuet,erru
       call polin2(xx,qq,redepl,nx,nq,delt,qsq,resdet,erru)
c      print *,"test7",icountxq,testt,delt,qsq,resdet,erru
       call polin2(xx,qq,resepl,nx,nq,delt,qsq,resset,erru)
c      print *,"test7",icountxq,testt,delt,qsq,resset,erru
       call polin2(xx,qq,imuepl,nx,nq,delt,qsq,resu1et,erru)
c      print *,"test7",icountxq,testt,delt,qsq,resu1et,erru
       call polin2(xx,qq,imdepl,nx,nq,delt,qsq,resd1et,erru)
c      print *,"test7",icountxq,testt,delt,qsq,resd1et,erru
       call polin2(xx,qq,imsepl,nx,nq,delt,qsq,ress1et,erru)
c      print *,"test7",icountxq,testt,delt,qsq,ress1et,erru

       call polin2(xx,qq,reuept3d,nx,nq,delt,qsq,resu2et,erru)
c      print *,"test8",icountxq,testt,delt,qsq,resu2et,erru
       call polin2(xx,qq,redept3d,nx,nq,delt,qsq,resd2et,erru)
c      print *,"test8",icountxq,testt,delt,qsq,resd2et,erru
       call polin2(xx,qq,resept3d,nx,nq,delt,qsq,ress2et,erru)
c      print *,"test8",icountxq,testt,delt,qsq,ress2et,erru
       call polin2(xx,qq,imuept3d,nx,nq,delt,qsq,resu3et,erru)
c      print *,"test8",icountxq,testt,delt,qsq,resu3et,erru
       call polin2(xx,qq,imdept3d,nx,nq,delt,qsq,resd3et,erru)
c      print *,"test8",icountxq,testt,delt,qsq,resd3et,erru
       call polin2(xx,qq,imsept3d,nx,nq,delt,qsq,ress3et,erru)
c      print *,"test8",icountxq,testt,delt,qsq,ress3et,erru
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


c--------LOLO
c       reutw3   = 0D0
c       imutw3   = 0D0
c       redtw3   = 0D0
c       imdtw3   = 0D0
c       restw3   = 0D0
c       imstw3   = 0D0
c       reuetw3  = 0D0
c       redetw3  = 0D0
c       resetw3  = 0D0
c       imuetw3  = 0D0
c       imdetw3  = 0D0
c       imsetw3  = 0D0
c       reuptw3  = 0D0
c       redptw3  = 0D0
c       resptw3  = 0D0
c       imuptw3  = 0D0
c       imdptw3  = 0D0
c       imsptw3  = 0D0
c       reueptw3 = 0D0
c       redeptw3 = 0D0
c       reseptw3 = 0D0
c       imueptw3 = 0D0
c       imdeptw3 = 0D0
c       imseptw3 = 0D0

      return

      end

C 
C     functions to compute the various contributions of DVCS, interference and BH 
C     to the UNP, LP and TP Xsections, separated by angular dependence. Also used to form  
C     asymmetries 
C 
      real*8 function tdvcsup(laml,lamp,jj,ii) 
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
c      print *,tiup,h1re,fact,f2p,(2.-2.*y+y**2)
     
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

 
C 
C   function that computes minus 
C   the value of the product of the two BH propagators! 
C   i.e. cor = - P_1 * P_2
C 
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
 
 
C 
C 
C     function that sets the t dependence of the various amplitudes 
 
      subroutine forms (t) 
      implicit none 
      integer marray
      parameter(marray = 40)
      real*8 t,spin1,dy 
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
      real*8 bnew,bnew_g
      
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
      else
       f1p = (gep - t*gmp/(4.*mp1**2))/(1.-t/(4.*mp1**2))
       f2p = (gmp - gep)/(1.-t/(4.*mp**2)) 
C
C     allow for different t-dependences for collider settings
C     Q^2 dependent slopes used for colliders
C
       if (del.gt.0.01D0) then
          f1u = (2.*f1p + f1n)/2.
       else
          f1u = dexp(t*(9D0*(1D0-0.15D0*dlog(qsq/2D0))))
       endif

        bnew = 7.6D0*(1D0-0.15D0*dlog(qsq/2D0))
        bnew = bnew/2.d0
        f1u  = dexp(t*bnew)
 
      endif
      
      
       if (del.gt.0.01D0) then 
ccc       print *,'here'
      f2u = (2.*f2p + f2n)/2. 
      f1d = f1p + 2.*f1n 
      f2d = f2p + 2.*f2n 
      f1s = (ges - t/(4.*mv**2)*gms)/(1-t/(4.*mv**2)) 
      f2s = (gms - ges)/(1.-t/(4.*mv**2)) 
        else
       f2u = f1u
       f1d = f1u
       f2d = f1u
       f1s = f1u
       f2s = f1u
       endif
      g1   = 1./(1.-t/ma**2)**2         
      g1sea = 1./(1.-t/ma**2)**3         
      gpi = 4 * ga3 * mp1**2/(mpi**2 -t)
      if (del.gt.0.01D0) then
         exg = 1./(1.-t/ma**2)**3
      else
         exg = dexp(t*(4.5D0*(1D0-0.15D0*dlog(qsq/2D0))))
      endif
	 
        f2u = f1u
        f1d = f1u
        f2d = f1u
        f1s = f1u
        f2s = f1u
       
       bnew_g = bnew
       bnew = bnew_g
       bnew = bnew/2.d0
       exg  = dexp(t*bnew)
      
      return 
      end 

C
C     2-dim interpolation routine
C

C****************************************************************
      subroutine polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
C****************************************************************

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

      


 
 

 
 
 
 
 
 
 
 
 
 
 
 
 
