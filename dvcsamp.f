CMM   Convolution program for NLO DVCS amplitude see Belitsky et al, hep-ph/9908337
CMM   This code takes as input the ODPDFs  
CMM   see e.g. Golec-Biernat Martin hep-ph/9807497
CMM
CMM   We always take quark singlet as the input, in both polarised and unpolarised cases
CMM   from the evolution code of A Freund. There is one input file for each skewedness value.
CMM   The code convolutes them with the appropriate coefficient functions to produce contributions to 
CMM   the DVCS amplitude from each flavour. These are stored in the output files 
CMM   uamp.dat, damp.dat, samp.dat, gamp.dat: the output is of the form
CMM
CMM   del(1)
CMM   q^2(1) ReA ImA
CMM   .
CMM   .
CMM   q^2(nq) ReA ImA
CMM   del(2)
CMM   .
CMM   .
CMM   del(ndel)
CMM   q^2(1) ReA ImA
CMM   .
CMM   .
CMM   q^2(1) ReA ImA
CMM
CMM  For use by dvcsob.f which includes Bether-Heitler and implements BMNS hep-ph/0004059
CMM   
CMM  How this code works:
CMM   
CMM   1) Subroutine readin(nx,nq,ndel) is called and returns the number of points in X,Q and delta (skewedness)
CMM      It stores the values of the ODPDFS to three dimensional 
CMM      arrays in del' Q' and X' and the values of delta, Q and X in 1D vectors. 
CMM      This are all passed via common blocks /pdfarray/, /delarray/, /qarray/ and /xarray/ 
CMM      to the main program.
CMM   2) The array in x is extended to have 0 and 1 at either end. Output files are opened.  
CMM   3) The calculations are performed in three nested do loops in delta, q and x' respectively. 
CMM      The delta loop specifies the value of delta from the array and stores it to the output.
CMM      The position of the cross over point between DGLAP and ERBL is calculated,
CMM       and the array index of the last point in the the ERBL region is stored in iv.
CMM      In the Q do loop the value of Q is specified, then the GPDs are extrapolated to 
CMM       the point x' = delta and then used to set the GPDs at the point x' = 0 using symmetries. 
CMM       This is achieved using a seven point array around iv and subroutine ratint
CMM      Subroutine wcread calculates the coeff functions for the particular delta and Q  
CMM       (using function tqi tgi and subroutine alpdr to calculate alphas). 
CMM        Argument n sets the order (1 or 2). 
CMM        It passes the coeff functions via common block /wcoeff/ tq, tg, tqm, tgm, tqi, tgi 
CMM        in arrays in x'.
CMM      A patten of minus signs is set (via sq,sg) (according to whether we want axial (iknl =2) 
CMM      or vector(iknl=1) case)
CMM      The 'subtracted' integrands are then defined (taking care of the necessary subtractions 
CMM          (at x' = 0)).
CMM      Using pu,pd,ps,pg for ERBL and pu1,..,pg1 for (just ?) DGLAP in integral 0..1 
CMM      The imaginary piece of the subtracted integrands of the integral del...1 are 
CMM      stored in pui,...,pgi. 
CMM      The special point at x' = delta is treated separately using the symmetries of the 
CMM       coefficient fns.
CMM      The special bin between del and x_iv+1 in the del..1 integral is then calculated 
CMM       and stored in temu1,...temg1 for the real part and temui,...,temgi for the imaginary 
CMM       part.
CMM     Do loop in x' then does the integrations (using gauss) in all the other bins and 
CMM      adds them cumilatively into temu1...
CMM      The arguments of integration routines "gauss" are not the above arrays but 
CMM      interpolations of them using fintrp. 
CMM     The counterterms are then calculated using subroutine counterV or counterA as 
CMM      appropriate to axial (iknl =2) or vector (iknl =1)
CMM     Finally, the real and imaginary peices are put together in reau,...,reag, ximau,....,ximag
CMM     and sent to the output files together with the value of Q^2. 
CMM

      subroutine convolution(tw,gpdtype,unpolopol,n,lambda11,lambda21
     >,nfl,nx1,ndel1,nq1)

c      program convol6
      implicit none
      real*8 del,qone,mu,reqint,imqint,regint,imgint
      integer j,k,ni,li,irt,l
      integer mxx,mq,mdel,m,m1,kt,tw
      integer k1,ktemp,ktemp1,old1
      integer nfl,n,iknl,gpdtype, unpolopol
      integer nx,nq,ndel,nx1,nq1,ndel1
      real*8 nf,lambda1,lambda2,lambda11,lambda21
      real*8 qu,qd,qs,g
      real*8 qar,delar,xar,xarb
      real*8 pu,pd,ps,pg,pu1,ps1,pd1,pg1
      real*8 pui,pdi,psi,pgi
      real*8 tq,tg,tqm,tgm,tqi,tgi
      integer iv,nc1,nc2,ncut,ntemp,iset,iset1,iset2,iset3
      real*8 udel,ddel,sdel,gdel,q2,e,aer,rer,err
      real*8 temu,temd,tems,temg,temu1,temd1,tems1,temg1
      real*8 temui,temdi,temsi,temgi,tempui,tempdi,tempsi
      real*8 gauss,tempu,tempd,temps,tempg,tempgi
      real*8 reau,read,reas,reag,dxdzn,dxdz
      real*8 ximau,ximad,ximas,ximag
      real*8 x0,tem,dz,diffdel,res,comp,comp1
      real*8 sq,sg,qm,smpnol,smpnor,smpsna
      real*8 temp1qu,temp1qd,temp1qs,temp1g,xtemp
      parameter (mxx = 1050, mq = 100, mdel = 100)
CMM   arrays of dimension mxx + 1
C maximum number of points in x, 1050, in x_bj, 20 and in Q^2, 20.
C can be adjusted of course
      common / pdfarray / qu(mdel,mq,0:mxx),qd(mdel,mq,0:mxx),
     >                    qs(mdel,mq,0:mxx),g(mdel,mq,0:mxx)
      common / qarray / qar(mq)
      common / delarray / delar(mdel)
      common / xarray1 / xarb(2,0:mxx,mdel)
      common / xarray / xar(0:mxx)
      common / regfunc / pu(mxx), pd(mxx), ps(mxx), pg(mxx),
     >               pu1(mxx), pd1(mxx), ps1(mxx), pg1(mxx),
     >               pui(mxx), pdi(mxx), psi(mxx), pgi(mxx)
      common / wcoeff / tq(0:mxx), tg(0:mxx), tqm(0:mxx), tgm(0:mxx),
     >                  tqi(0:mxx),tgi(0:mxx)
      common / all / nx,nq,ndel
      common / clambda / lambda1,lambda2
      common / tempar / temp1qu(mdel),temp1qd(mdel),temp1qs(mdel)
     >                  ,temp1g(mdel),xtemp(mdel)
      dimension dxdz(mxx)


C
C     external functions
C

      external alpdr
      external spence
      external counterV,counterA
      external smpnol,smpnor,smpsna,dxdzn

      nx = nx1
      ndel = ndel1
      nq = nq1
      lambda1 = lambda11
      lambda2 = lambda21
      iknl = unpolopol
      k1 = gpdtype


C     load pdf arrays, x array, del and Q value, number of bins, del's and Q's

CMM   Returns values for nx,nq,ndel and
CMM   Fills up 3D-arrays qu,qd,qs,g dim (ndel,nq,nx)
CMM   and 1D arrays delar(ndel), qar(nq), xar(nx)
CMM     define nf as real*8 otherwize dangerous integer division

c      do 99 k1 = 1,gpdtype

c      do 100 n = 1,2

c         do 101 iknl = 1,unpolopol

C
C set number of flavours
C

      nf = dble(nfl)

C
C     readin GPD data files
C

      call readin(k1,n,iknl)

C     opening file for the values of the amplitude for the various quarks and the gluon
C     output mode: first the q value than in the next line del, Real part of amplitude, 
C     Imaginary part of the amplitude. various files for the quark species and the gluon. 


      if (tw.eq.2) then

      if (k1.eq.1) then

      if (n.eq.1.and.iknl.eq.1) then

      open(1,file='luamp.dat',status='unknown')
      open(2,file='ldamp.dat',status='unknown')
      open(3,file='lsamp.dat',status='unknown')
      open(4,file='lgamp.dat',status='unknown')

      elseif(n.eq.1.and.iknl.eq.2) then

      open(1,file='luamppol.dat',status='unknown')
      open(2,file='ldamppol.dat',status='unknown')
      open(3,file='lsamppol.dat',status='unknown')
      open(4,file='lgamppol.dat',status='unknown')
      
      elseif(n.eq.2.and.iknl.eq.1) then

      open(1,file='nlouamp.dat',status='unknown')
      open(2,file='nlodamp.dat',status='unknown')
      open(3,file='nlosamp.dat',status='unknown')
      open(4,file='nlogamp.dat',status='unknown')

      elseif(n.eq.2.and.iknl.eq.2) then

      open(1,file='nlouamppol.dat',status='unknown')
      open(2,file='nlodamppol.dat',status='unknown')
      open(3,file='nlosamppol.dat',status='unknown')
      open(4,file='nlogamppol.dat',status='unknown')

      endif

      else

      if (n.eq.1.and.iknl.eq.1) then

      open(1,file='luampe.dat',status='unknown')
      open(2,file='ldampe.dat',status='unknown')
      open(3,file='lsampe.dat',status='unknown')
      open(4,file='lgampe.dat',status='unknown')

      elseif(n.eq.1.and.iknl.eq.2) then

      open(1,file='luamppole.dat',status='unknown')
      open(2,file='ldamppole.dat',status='unknown')
      open(3,file='lsamppole.dat',status='unknown')
      open(4,file='lgamppole.dat',status='unknown')
      
      elseif(n.eq.2.and.iknl.eq.1) then

      open(1,file='nlouampe.dat',status='unknown')
      open(2,file='nlodampe.dat',status='unknown')
      open(3,file='nlosampe.dat',status='unknown')
      open(4,file='nlogampe.dat',status='unknown')

      elseif(n.eq.2.and.iknl.eq.2) then

      open(1,file='nlouamppole.dat',status='unknown')
      open(2,file='nlodamppole.dat',status='unknown')
      open(3,file='nlosamppole.dat',status='unknown')
      open(4,file='nlogamppole.dat',status='unknown')

      endif

      endif

      else
      
      if (k1.eq.1) then

      if (n.eq.1.and.iknl.eq.1) then

      open(1,file='luamptw3.dat',status='unknown')
      open(2,file='ldamptw3.dat',status='unknown')
      open(3,file='lsamptw3.dat',status='unknown')
      open(4,file='lgamptw3.dat',status='unknown')

      elseif(n.eq.1.and.iknl.eq.2) then

      open(1,file='luamppoltw3.dat',status='unknown')
      open(2,file='ldamppoltw3.dat',status='unknown')
      open(3,file='lsamppoltw3.dat',status='unknown')
      open(4,file='lgamppoltw3.dat',status='unknown')
      
      endif

      else

      if (n.eq.1.and.iknl.eq.1) then

      open(1,file='luampetw3.dat',status='unknown')
      open(2,file='ldampetw3.dat',status='unknown')
      open(3,file='lsampetw3.dat',status='unknown')
      open(4,file='lgampetw3.dat',status='unknown')

      elseif(n.eq.1.and.iknl.eq.2) then

      open(1,file='luamppoletw3.dat',status='unknown')
      open(2,file='ldamppoletw3.dat',status='unknown')
      open(3,file='lsamppoletw3.dat',status='unknown')
      open(4,file='lgamppoletw3.dat',status='unknown')

      endif

      endif
         

      endif

C
C     loop that runs through the various values of del
C

      do 22 k = 1, ndel

         del = delar(k)

C
C     write the del values in output
C

         write(1,*) del
         write(2,*) del
         write(3,*) del
         write(4,*) del

CMM   Now determine iv for this delta (corssover point from ERBL to DGLAP)

         iv = 0

      do 121 j = 1, nx

         xar(j) = xarb(k1,j,k)

         if (xar(j).le.del) then

            iv = iv + 1

            else

            iv = iv

         endif

 121        continue

       xar(nx+1) = 1.0

C
C     initialize jacobain array used in simpson integration!
C

       do 144 j = 1,nx+1

       dz = 1./nx
       x0 = dz*(j-1)

       if (j.le.iv) then

          dxdz(j) = dxdzn(1,j,nx+1,x0,del)

          else

          dxdz(j) = dxdzn(2,j,nx+1,x0,del)
       
       endif

       IF (J.EQ.IV) DIFFDEL = DXDZN(2,J,NX+1,X0,DEL)

 144   continue

      

C
C     outer loop that runs through the various values of q
C

         do 11 ni = 1,nq

         q2 = qar(ni)

C
C Fill the array with the values of the LO and NLO coefficient functions!
C
C
C     Set ratio of factorization to renormalization scale in variable qone!
C
C     Default value Q^2/mu^2 =1 
C
           mu = q2
           qone = q2

           qm = qone**2/mu**2
           
         call wcread(tw,del,qm,q2,n,iknl,nf,nx,iv)

C
C     inner loop that runs through x for the integration
C

         e = del/(2.-del)

C
C     Set switch for real part of T(2*y/del-1) for y>del in the integrand of 
C     the integral from del..1
C     depending whether vector or axial-vector is used (IKNL=1,2)
C     unpol vs pol
C
         if (iknl.eq.1) then

            sq = -1.
            sg = +1.

         else if (iknl.eq.2) then

            sq = +1.
            sg = -1.
            
         else
            stop

         endif

         if (k1.eq.2) then

            if (iknl.eq.1) then
              qu(k,ni,iv) =  -qu(k,ni,1)
              qd(k,ni,iv) =  -qd(k,ni,1)
              qs(k,ni,iv) =  -qs(k,ni,1)
              g(k,ni,iv) =  g(k,ni,1)
              else
              qu(k,ni,iv) =  qu(k,ni,1)
              qd(k,ni,iv) =  qd(k,ni,1)
              qs(k,ni,iv) =  qs(k,ni,1)
              g(k,ni,iv) =  -g(k,ni,1)
              endif
          endif

CMM   defining subtracted integrands for both integrals
CMM   tqm is T^q in integrand of  integral from del to 1


         do 33 l = 1, nx+1

            if (l.lt.iv) then
               
           pu(l) = dxdz(l)*tq(l)*4./9.*(qu(k,ni,l)-qu(k,ni,iv))/e
           pd(l) = dxdz(l)*1./9.*tq(l)*(qd(k,ni,l)-qd(k,ni,iv))/e 
           ps(l) = dxdz(l)*1./9.*tq(l)*(qs(k,ni,l)-qs(k,ni,iv))/e 
           pg(l) = dxdz(l)*1./nf*tg(l)*(g(k,ni,l)-g(k,ni,iv))/(e**2)

            elseif(l.gt.iv) then

CMM     Same thing but also define the integrands for the second integral. 
CMM     Combine the original integral 
C       del..1 with the real part of the integrand T(2*y/del-1)*(q(y)-qdel).

               pu1(l-iv+1) = dxdz(l)*(tqm(l)*4./9.*qu(k,ni,l)/e + 
     >                  sq*tq(l)*4./9.*(qu(k,ni,l)-qu(k,ni,iv))/e)
               pd1(l-iv+1) = dxdz(l)*(1./9.*tqm(l)*qd(k,ni,l)/e + 
     >                  sq*1./9.*tq(l)*(qd(k,ni,l)-qd(k,ni,iv))/e) 
               ps1(l-iv+1) = dxdz(l)*(1./9.*tqm(l)*qs(k,ni,l)/e + 
     >                  sq*1./9.*tq(l)*(qs(k,ni,l)-qs(k,ni,iv))/e) 
               pg1(l-iv+1) = dxdz(l)*(1./nf*tgm(l)*g(k,ni,l)/(e**2) + 
     >                  sg*1./nf*tg(l)*(g(k,ni,l)-g(k,ni,iv))/(e**2))

CMM  calculate imaginary piece from del...1 in first integral
CMM  Subtractions needed here too

        pui(l-iv+1) = dxdz(l)*tqi(l)*4./9.*(qu(k,ni,l)- qu(k,ni,iv))/e
        pdi(l-iv+1) = dxdz(l)*1./9.*tqi(l)*(qd(k,ni,l)- qd(k,ni,iv))/e 
        psi(l-iv+1) = dxdz(l)*1./9.*tqi(l)*(qs(k,ni,l)- qs(k,ni,iv))/e 
        pgi(l-iv+1) = dxdz(l)*1./nf*tgi(l)*(g(k,ni,l)-g(k,ni,iv))/e**2
           
            endif

 33         continue  

CMM   value of integrand at y=del for the integral del..1
CMM   Special point...not calculated above
CMM   Algabraic expressions for tq(1) are given specially in wcread
CMM   See notes for the calculation.
CMM   Here we use the fact that tq(1) = tqm(iv), i.e. T(y=0) = Tm (y=del)
CMM   To LO T = 1/(1-x) Tm = 1/(1+x), x = 2y/del -1 to T(y=0) = 1/2 =Tm(y=del)
CMM   The use of pu1(iv) is just for code convenience, these integrands correspond to yy = del > xar(iv)
C
C Since semi-open simpson is used these points are not directly used in the integration but are still put here
C for bookkeeping purposes!
C

            pu1(1) = diffdel*tq(1)*4./9.*qu(k,ni,iv)/e
            pd1(1) = diffdel*1./9.*tq(1)*qd(k,ni,iv)/e
            ps1(1) = diffdel*1./9.*tq(1)*qs(k,ni,iv)/e 
            pg1(1) = diffdel*1./nf*tg(1)*g(k,ni,iv)/(e**2)

C
C     Values of the integrand for 0..del integral at y=del
C
            pu(iv) = 0.0
            pd(iv) = 0.0
            ps(iv) = 0.0
            pg(iv) = 0.0


C
C    value of the integrands of the imaginary part in the integral from del..1.
C    at y=del.
C

            pui(1) = 0.0
            pdi(1) = 0.0
            psi(1) = 0.0 
            pgi(1) = 0.0

C
C     \int^1_0 dy T(2*y/del-1)*(q(y)-q(del) -/+ \int^1_del dy T(1-2*y/del)*q(y)
C     and reverse the sign for the gluon.
C

C
C     Set # of points in ERBL and DGLAP region,first ERBL then DGLAP
C

            nc1 = iv
            nc2 = nx-iv+2
C
C     spacing of points
C
            dz = 1./nx

C
C     do the integration from 0..del
C

          temu = smpnor(nc1,dz,pu,err)

          temd = smpnor(nc1,dz,pd,err)
          
          tems = smpnor(nc1,dz,ps,err)
          
          if (n.eq.2) then

          temg = smpnor(nc1,dz,pg,err)

          endif

C
C     DO the integration from del..1
C

        temu1 = smpnol(nc2,dz,pu1,err)

        temd1 = smpnol(nc2,dz,pd1,err)

        tems1 = smpnol(nc2,dz,ps1,err)

          if (n.eq.2) then

        temg1 = smpnol(nc2,dz,pg1,err)

          endif

CMM    Separate calc for for V and A not needed. 
CMM    Correct choice made within wcread

C
C Integration for imaginary part. del..1 only!
C

          if (n.eq.2) then

        temui = smpnol(nc2,dz,pui,err)

        temdi = smpnol(nc2,dz,pdi,err)

        temsi = smpnol(nc2,dz,psi,err)     

        temgi = smpnol(nc2,dz,pgi,err)

        elseif (n.eq.1.and.tw.eq.3) then

        temui = smpnol(nc2,dz,pui,err)

        temdi = smpnol(nc2,dz,pdi,err)

        temsi = smpnol(nc2,dz,psi,err)     

        temgi = 0D0

          endif
              
CMM      Unpolarized case
C
C     Initialize variables
C
           reqint = 0.

           imqint = 0.

           regint = 0.

           imgint = 0.

C
C     compute value of the counterterms for unpolarized case
C

           if (iknl.eq.1) then
              
      call counterV(tw,del,qone,mu,nf,reqint,imqint,regint,imgint,n)            

C
C     Assemble final value of the real parts
C

           reau = temu - temu1 + 4./9.*qu(k,ni,iv)*reqint/e

           read = temd - temd1 + 1./9.*qd(k,ni,iv)*reqint/e

           reas = tems - tems1 + 1./9.*qs(k,ni,iv)*reqint/e

           if (n.eq.2) then

           reag = temg + temg1 + 1./(nf*e**2)*g(k,ni,iv)*regint

           else

           reag = 0D0

           endif

C
C     Assembly of imaginary part!
C

C
C     There are no integrals for the imaginary part in LO (n=1) twist-2
C

           If (n.eq.1.and.tw.eq.2) then

              temui = 0.0

              temdi = 0.0

              temsi = 0.0

              temgi = 0.0

           endif

C
C     put together counterterms and result of integral for imaginary part
C

           ximau = 4./9.*qu(k,ni,iv)*imqint/e + temui

           ximad = 1./9*qd(k,ni,iv)*imqint/e  + temdi

           ximas = 1./9.*qs(k,ni,iv)*imqint/e + temsi

           if (n.eq.2) then

           ximag = 1./(nf*e**2)*g(k,ni,iv)*imgint + temgi

           else

           ximag = 0D0

           endif

C
C     Same as above but now for the polarized amplitudes
C

           else if(iknl.eq.2) then

      call counterA(tw,del,qone,mu,nf,reqint,imqint,regint,imgint,n)


           reau = temu + temu1 + 4./9.*qu(k,ni,iv)*reqint/e

           read = temd + temd1 + 1./9.*qd(k,ni,iv)*reqint/e

           reas = tems + tems1 + 1./9.*qs(k,ni,iv)*reqint/e

           if (n.eq.2) then 

           reag = temg - temg1 + 1./(nf*e**2)*g(k,ni,iv)*regint

           else

           reag = 0D0

           endif

           if (n.eq.1) then

              temui = 0.0

              temdi = 0.0

              temsi = 0.0

              temgi = 0.0

              endif

           ximau = 4./9.*qu(k,ni,iv)*imqint/e + temui

           ximad = 1./9*qd(k,ni,iv)*imqint/e  + temdi

           ximas = 1./9.*qs(k,ni,iv)*imqint/e + temsi

           if (n.eq.2) then

           ximag = 1./(nf*e**2)*g(k,ni,iv)*imgint + temgi

           else

           ximag = 0D0

           endif

           else
              write(6,*) 'iknl not rec. stop'
              stop
           endif
           
C
C     Write output of final results in files
C

           write(1,102) q2**2, reau, ximau
           write(2,102) q2**2, read, ximad
           write(3,102) q2**2, reas, ximas
           write(4,102) q2**2, reag, ximag


 11     continue
 22     continue
 

        close(28)
        close(4)
        close(3)
        close(2)
        close(1)

c 101    continue
c 100    continue
c 99    continue
        return

 102    FORMAT(3(E15.8,1X))
        end


C
C     Dick Robert's code for alpha_s
C
      real*8 function alpdr(n,nf,q)
      implicit none
      external alp
      integer n,iord
      real*8 nf
      real*8 qsct,qsdt,alambda,xflav,lambda1,lambda2
      common / clambda / lambda1, lambda2 
      common/alpin/alambda,xflav,qsct,qsdt,iord
      real*8 al2,q2,t,alp,q
      iord =n-1
      xflav = dble(nf)
      qsdt = 6.9696D0  !! mc = 1.32 values in PDG 2000
      qsct = 73.96D0   !! mb = 4.3

      alambda = lambda2
      if (iord.eq.0) then
          alambda = lambda1
      endif    
      al2=alambda*alambda
      q2 = q*q
      t=dlog(q2/al2)
      alpdr=alp(t)
      return
      end 

      real*8  FUNCTION	alp(T)
      implicit none
      REAL*8 alambda,xflav,qsct,qsdt,pi
      REAL*8 tol,tt,t,qsdtt,qsctt,al,al2,qs
      integer iord,ith
      real*8 b0,b1,f,fp,del,alfqc4,alfqc5,alfqs5,alfinv
      real*8 x1,x2,as2,as,alfqs3,alfqc3
      COMMON/alpin/alambda,xflav,qsct,qsdt,iord
      DATA PI/3.141592653589793/
      DATA TOL/.00005/
      ITH=0
      TT=T
      qsdtt=qsdt/4.
      qsctt=qsct/4.    
      AL=alambda
      AL2=AL*AL
      XFLAV=4.
      QS=AL2*EXP(T)

      if(qs.lt.0.5d0) then   !!  running stops below 0.5
          qs=0.5d0
          t=dlog(qs/al2)
          tt=t
      endif

      IF(QS.gt.QSCTT) GO	TO 12  
      IF(QS.lt.QSDTT) GO	TO 312  
   11 CONTINUE
      B0=11-2.*XFLAV/3. 
      X1=4.*PI/B0
      IF(IORD.eq.0) then
      ALP=X1/T
      ELSE
      B1=102.-38.*XFLAV/3.
      X2=B1/B0**2
      AS2=X1/T*(1.-X2*dlog(T)/T)
    5 AS=AS2
      F=-T+X1/AS-X2*dlog(X1/AS+X2)
      FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
      AS2=AS-F/FP
      DEL=ABS(F/FP/AS)
      IF((DEL-TOL).GT.0.) go to 5
      ALP=AS2
      ENDIF
      IF(ITH.EQ.0) RETURN
      GO TO (13,14,15) ITH
      GO TO 5
   12 ITH=1
      T=dlog(QSCTT/AL2)
      GO TO 11
   13 ALFQC4=ALP
      XFLAV=5.   
      ITH=2
      GO TO 11
   14 ALFQC5=ALP
      ITH=3
      T=TT
      GO TO 11
   15 ALFQS5=ALP
      ALFINV=1./ALFQS5+1./ALFQC4-1./ALFQC5
      ALP=1./ALFINV
      RETURN

  311 CONTINUE
      B0=11-2.*XFLAV/3. 
      X1=4.*PI/B0
      IF(IORD.eq.0) then
      ALP=X1/T
      ELSE
      B1=102.-38.*XFLAV/3.
      X2=B1/B0**2
      AS2=X1/T*(1.-X2*dlog(T)/T)
   35 AS=AS2
      F=-T+X1/AS-X2*dlog(X1/AS+X2)
      FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
      AS2=AS-F/FP
      DEL=ABS(F/FP/AS)
      IF((DEL-TOL).GT.0.) go to 35
      ALP=AS2
      endif
      IF(ITH.EQ.0) RETURN
      GO TO (313,314,315) ITH
  312 ITH=1
      T=dlog(QSDTT/AL2)
      GO TO 311
  313 ALFQC4=ALP
      XFLAV=3.   
      ITH=2
      GO TO 311
  314 ALFQC3=ALP
      ITH=3
      T=TT
      GO TO 311
  315 ALFQS3=ALP
      ALFINV=1./ALFQS3+1./ALFQC4-1./ALFQC3
      ALP=1./ALFINV
      RETURN
      END



C
C     Coefficient functions for LO tq1 and NLO tq2,tq2a,tg2,tg2a
C
CMM   See eq.(14-17) in BMNS hep-ph/9908337

      real*8 function tq1(xx,qm)
      implicit none
      real*8 tq2,tq2a,tq3r,tq3i,tq3ar,tq3ai
      real*8 tg2,tg2a,tg3r,tg3i,tg3ar,tg3ai
      real*8 tq1t3,tq1r,tq1i
      real*8 xx, qm
      real*8 pi
      PARAMETER (pi = 3.141592653589793)

      tq1 = 1./(1.-xx)

      return

      entry tq1t3(xx,qm)

      tq1t3 = (log(2.)-log(1.-xx))/(1.+xx)

      return

      entry tq2(xx,qm)

      tq2 = ((2.*log((1.-xx)/2.) + 3.)*(log(qm) +
     > log((1.-xx)/2.)/2. - 3./4.) - 27./4.)/(2.*(1.-xx)) -
     > 3.*log((1.-xx)/2.)/(2.*(1.+xx))

      return

      entry tq2a(xx,qm)

      tq2a = ((2.*log((1.-xx)/2.) + 3.)*(log(qm) +
     > log((1.-xx)/2.)/2. - 3./4.) - 27./4.)/(2.*(1.-xx)) -
     > log((1.-xx)/2.)/(2.*(1.+xx))

      return

CMM   Needed for del....1, i.e. x >1
CMM   since xx > 1, t > xi now 
CMM   take minus sign out of log according to the sign of the +ie prescription
C
C     **3*r = real part of coefficient function
C     **3*i = imaginary part of coefficient function
C

      entry tq1r(xx,qm)

      tq1r = (log(2.)-log(xx-1.))/(1.+xx)

      return

      entry tq1i(xx,qm)

      tq1i = pi/(1.+xx)

      return

      entry tq3r(xx,qm)

      tq3r = ((2.*log((xx-1.)/2.) + 3.)*(log(qm) +
     > log((xx-1.)/2.)/2. - 3./4.) - 27./4. - pi**2)/(2.*(1.-xx))  -
     > 3.*log((xx-1.)/2.)/(2.*(1.+xx))

      return

      entry tq3i(xx,qm)

      tq3i = -pi/2.*(2.*(log(qm)+log((xx-1.)/2.))/(1.-xx) -
     > 3./(1.+xx))
      
      return

      entry tq3ar(xx,qm)

      tq3ar = ((2.*log((xx-1.)/2.) + 3.)*(log(qm) +
     > log((xx-1.)/2.)/2. - 3./4.) - 27./4. - pi**2)/(2.*(1.-xx)) -
     > log((xx-1.)/2.)/(2.*(1.+xx))

      return

      entry tq3ai(xx,qm)

      tq3ai = -pi/2.*(2.*(log(qm)+log((xx-1.)/2.))/(1.-xx)  -
     > 1./(1.+xx))

      return

      entry tg2(xx,qm)

      tg2 = -((1./(1.-xx**2) + log((1.-xx)/2.)/(1.+xx)**2)*
     > (log(qm) + log((1.-xx)/2.) - 2.) -
     > (log((1.-xx)/2.)/(1.+xx))**2/2.                    )/2. + 
     > ((log(qm) + log((1.-xx)/2.) - 2.)/(1.-xx) +  
     > log((1.-xx)/2.)/(1.+xx))/2.

      return

      entry tg2a(xx,qm)

      tg2a = ((1./(1.-xx**2) + log((1.-xx)/2.)/(1.+xx)**2)*
     > (log(qm) + log((1.-xx)/2.) - 2.) - 
     > (log((1.-xx)/2.)/(1.+xx))**2/2.)/2. 

      return

      entry tg3r(xx,qm)

      tg3r = -((1./(1.-xx**2) + log((xx-1.)/2.)/(1.+xx)**2)*
     > (log(qm) + log((xx-1.)/2.) - 2.) -
     > (log((xx-1.)/2.)/(1.+xx))**2/2. - (pi/(1.+xx))**2/2.)/2.  +
     > ((log(qm) + log((xx-1.)/2.) - 2.)/(1.-xx)  +
     > log((xx-1.)/2.)/(1.+xx))/2.

      return

      entry tg3i(xx,qm)

      tg3i = -pi*(1./(1.-xx**2) -  
     > (log(qm) +log((xx-1.)/2.) -2.)/(1.+xx)**2 )/2.

      return

      entry tg3ar(xx,qm)

      tg3ar =((1./(1.-xx**2) + log((xx-1.)/2.)/(1.+xx)**2)*
     > (log(qm) + log((xx-1.)/2.) - 2.)  -
     > (log((xx-1.)/2.)/(1.+xx))**2/2. - 
     > (pi/(1+xx))**2/2.)/2.

      return

      entry tg3ai(xx,qm)

      tg3ai = -pi*((log(qm) + log((xx-1.)/2.) -2.)/(1.+xx)**2 
     > + 1./(1.-xx**2))/2.

      return
      end

C
C     subroutine to read in of GPDs, q's and del's
C
      subroutine readin(k1,n,iknl)
      implicit none
      integer ndel,nq,nx,n,iknl
      integer k1
      real*8 pi
      PARAMETER (pi = 3.141592653589793)
      integer mxx,mq,mdel,int,num
      integer k,ik,i,i1,ntemp,ntemp1,dig
      real*8 qu,qd,qs,g,qar,delar,xarb
      character*2 temp
      character*78 temp1
      parameter (mxx = 1050, mq = 100, mdel = 100)
      common / pdfarray / qu(mdel,mq,0:mxx),qd(mdel,mq,0:mxx),
     >                    qs(mdel,mq,0:mxx),g(mdel,mq,0:mxx)
      common / qarray / qar(mq)
      common / delarray / delar(mdel)
      common / xarray1 / xarb(2,0:mxx,mdel)
      common / all / nx,nq,ndel
       
C
C     opening of GPD data files
C

      if (k1.eq.1) then

      if (n.eq.1.and.iknl.eq.1) then

C
C     assembly of filename. up to 100 (0-99) files possible
C
         do 101 i = 1,ndel

         write(temp,333) i

         temp1 = 'ldel'//temp//'.dat'

         i1 = i+10

      open(i1,file=temp1,status='old')

 101  continue
      
      elseif(n.eq.1.and.iknl.eq.2) then

      do 102 i = 1,ndel

         write(temp,333) i

         temp1 = 'ldel'//temp//'pol.dat'

         i1 = i+10

      open(i1,file=temp1,status='old')

 102  continue


       elseif(n.eq.2.and.iknl.eq.1) then

          do 103 i = 1,ndel

         write(temp,333) i

         temp1 = 'ndel'//temp//'.dat'

         i1 = i+10

      open(i1,file=temp1,status='old')

 103  continue
      

      elseif(n.eq.2.and.iknl.eq.2) then

      do 104 i = 1,ndel

         write(temp,333) i

         temp1 = 'ndel'//temp//'pol.dat'

         i1 = i+10

      open(i1,file=temp1,status='old')

 104  continue

      endif

      else

      if (n.eq.1.and.iknl.eq.1) then

      do 105 i = 1,ndel

         write(temp,333) i

         temp1 = 'ldel'//temp//'e.dat'

         i1 = i+10

      open(i1,file=temp1,status='old')

 105  continue   
      
      elseif(n.eq.1.and.iknl.eq.2) then

         do 106 i = 1,ndel

         write(temp,333) i

         temp1 = 'ldel'//temp//'pole.dat'

         i1 = i+10

      open(i1,file=temp1,status='old')

 106  continue  

       elseif(n.eq.2.and.iknl.eq.1) then

          do 107 i = 1,ndel

         write(temp,333) i

         temp1 = 'ndel'//temp//'e.dat'

         i1 = i+10

      open(i1,file=temp1,status='old')

 107  continue  

      elseif(n.eq.2.and.iknl.eq.2) then

         do 108 i = 1,ndel

         write(temp,333) i

         temp1 = 'ndel'//temp//'pole.dat'

         i1 = i+10

      open(i1,file=temp1,status='old')

 108  continue  

      endif

      endif

C
C     run through the data files and write quarksinglet, gluon etc. into
C     respective arrays
C

C
C     each file corresponds to one value of del = x_bj
C      

        do 1 k = 1,ndel

           int = k + 10
C
C     readin x_bj
C

           read(int,*) delar(k)

           do 21 ik = 1, nq
C
C     readin q value
C
              read(int,*) qar(ik)

C
C     readin values of quarksinglets and gluons
C
           do 2 i = 1,nx

            read (int,*) xarb(k1,i,k),g(k,ik,i),qu(k,ik,i),
     >      qd(k,ik,i),qs(k,ik,i)

 2         continue

C
C     fix value of pdf at x=1, i.e. at 0!
C

            g(k,ik,nx+1) = 0.0
           qu(k,ik,nx+1) = 0.0
           qd(k,ik,nx+1) = 0.0
           qs(k,ik,nx+1) = 0.0

 21        continue

           close(int)
           
 1       continue

         return

 333     FORMAT(I2.2)

         end

C
C     read in of coefficint functions of either 1st or 2nd. order at a 
C     certain del=x_bj and q.
C
CMM   rewritten will only run if iknl =1,2 and n = 1,2 otherwize crashes
CMM   designed to put the right bit of the coeff fn into the integral concerned
CMM   use arrays tq,tg for the int_0^1 and tqm,tgm and tqi,tgi for re and im 
CMM   parts of the second integral int_del^1

      subroutine wcread(tw,del,qm,q2,n,iknl,nf,nx,iv)

      implicit none
      real*8 nf,del,q2
      integer n,iknl,nx,iv,tw
      real*8 tq,tg,tqm,tgm,xar
      real*8 tqi,tgi
      real*8 cf
      real*8 pi,qm,e,t
      integer i
      PARAMETER (pi = 3.141592653589793)
      integer mxx,mq,mdel
      parameter (mxx = 1050, mq = 100, mdel = 100)
      common / wcoeff / tq(0:mxx), tg(0:mxx), tqm(0:mxx), tgm(0:mxx),
     >                  tqi(0:mxx),tgi(0:mxx)
      common / xarray / xar(0:mxx)
      real*8   tq1,tq2,tg2,tq2a,tg2a,alpdr,alp,tq1t3,tq1r,tq1i
      external tq1,tq2,tg2,tq2a,tg2a,alpdr,tq1t3,tq1r,tq1i
      real*8   tq3r,tq3i,tg3r,tg3i,tq3ar,tq3ai,tg3ar,tg3ai
      external tq3r,tq3i,tg3r,tg3i,tq3ar,tq3ai,tg3ar,tg3ai

CMM   parameter qm = Q^2/mu^2 determines the factorization scale dependence.

      cf = 4./3.

CMM   This is eta NOT x_bj !

      e = del/(2.-del)

C
C     compute alphas_s/2pi
C

      alp = alpdr(n,nf,q2)/(2.*pi)

CMM   First loop is for 0...del of 0..1 integral
C     leaving out the points y=0 and y=del They are dealt with 
C     in the generation of the integrands in the main program!
C

      do 2 i = 2,iv-1

         t = (2.*xar(i) - del)/(2.-del)
C
C     n = 1 LO, n = 2 NLO
C
         if (n.eq.1) then

CMM   same for unpolarized and polarized cases.

            if (tw.eq.2) then 

             tq(i) =  tq1(t/e,qm)
             tg(i) = 0.0

             elseif (tw.eq.3) then
                
                tq(i) =  tq1t3(t/e,qm)
                tg(i) = 0.0

            endif    

         elseif(n.eq.2) then
C
C     difference between polarized and unpolarized coefficient functions in NLO
C
            if (iknl.eq.1) then

               tq(i) =  tq1(t/e,qm) + alp*cf*tq2(t/e,qm)
               tg(i) =                alp*nf*tg2(t/e,qm)

            elseif(iknl.eq.2) then

               tq(i)  = tq1(t/e,qm) + alp*cf*tq2a(t/e,qm)
               tg(i)  =               alp*nf*tg2a(t/e,qm)

            else

               write(6,*) 'iknl not rec. in wcread... stop'
               stop

            endif
         else

            write(6,*) 'n - order not rec. in wcread... stop'
            stop

         endif

 2       continue

C
C     fill arrays for coefficient functions for del..1 integral
C     f.ex. tq refers to the del..1 integral originating from the 0..1 
C     integral and tqm corresponds to the coefficient function in the del..1
C     integral where the argument of the fucntion is not 2y/x_bj - 1 but
C     rather 1-2y/x_bj. del = x_bj.
C

         do 1 i = iv+1,nx

            t = (2.*xar(i) - del)/(2.-del)

C
C     First LO (n=1). pol. and unpol. case have the same coefficient function
C     in this case! Coefficient function is onyl real!
C

            if (n.eq.1) then

               if (tw.eq.2) then
C
C     twist-2 coefficient functions
C

               if(iknl.eq.1) then

                  tq(i)  =  tq1(t/e,qm)
                  tg(i) = 0.
                  tqm(i) =  tq1(-t/e,qm)
                  tgm(i) = 0.

               elseif(iknl.eq.2) then

                  tq(i)  =  tq1(t/e,qm)
                  tqm(i) =  tq1(-t/e,qm)
                  tg(i) =  0.
                  tgm(i) = 0.
               else

                  write(6,*) 'iknl not rec. in wcread... stop'
                  stop

               endif

               elseif (tw.eq.3) then
C
C     twist-3 coefficient functions 
C

                  if(iknl.eq.1) then

C
C     real parts in LO
C

                  tq(i)  =  tq1r(t/e,qm)
                  tg(i) = 0.
                  tqm(i) =  tq1(-t/e,qm)
                  tgm(i) = 0.

C
C     imaginary part in LO
C

                  tqi(i) = tq1i(t/e,qm)
                  tgi(i) = 0.


               elseif(iknl.eq.2) then

C
C     real parts in LO
C

                  tq(i)  =  tq1r(t/e,qm)
                  tqm(i) =  tq1(-t/e,qm)
                  tg(i) =  0.
                  tgm(i) = 0.

C
C     imaginary parts in LO
C

                  tqi(i) = tq1i(t/e,qm)
                  tgi(i) = 0.

               else

                  write(6,*) 'iknl not rec. in wcread... stop'
                  stop

               endif


                  endif

C
C      Now NLO (n=2). There exist now also imaginary parts for the 
C      coefficient functions f.ex. tqi. Correct sig is determined
C      from the sign of the +ie prescription.
C

            elseif(n.eq.2) then

               if(iknl.eq.1) then

               tq(i) =  tq1(t/e,qm)  + alp*cf*tq3r(t/e,qm)
               tg(i) =                 alp*nf*tg3r(t/e,qm)
               tqm(i) = tq1(-t/e,qm) + alp*cf*tq2(-t/e,qm)
               tgm(i) =                alp*nf*tg2(-t/e,qm)

CMM    Imaginary parts               

               tqi(i) =                alp*cf*tq3i(t/e,qm)
               tgi(i) =                alp*nf*tg3i(t/e,qm)

               elseif(iknl.eq.2) then

               tq(i)  =   tq1(t/e,qm)  + alp*cf*tq3ar(t/e,qm)
               tg(i)  =                  alp*nf*tg3ar(t/e,qm)
               tqm(i) =   tq1(-t/e,qm) + alp*cf*tq2a(-t/e,qm)
               tgm(i) =                  alp*nf*tg2a(-t/e,qm)

CMM    Imaginary parts               

               tqi(i) =                  alp*cf*tq3ai(t/e,qm)
               tgi(i) =                  alp*nf*tg3ai(t/e,qm)

               else

                  write(6,*) 'iknl not rec. in wcread... stop'
                  stop

               endif
            else

               write(6,*) 'n - order not rec. in wcread... stop'
               stop
            endif

 1       continue

CMM    value of the coefficient function at y=0. otherwise 
CMM    same principle as above. 
CMM    tqm is supposed to be T(-t) but tq(y=del) = tqm(y=0) so calculate the 
CMM    latter instead.

         if (n.eq.1) then

            if (iknl.eq.1) then

               tq(1) =  1./2.
               tg(1) = 0.0
               tgm(1) = 0.0

            else if(iknl.eq.2) then

               tq(1)  = 1./2.
               tg(1) = 0.0
               tgm(1) = 0.0

            else
               write(6,*) 'iknl not rec. in wcread... stop'
               stop
            endif

         else if(n.eq.2) then

            if (iknl.eq.1) then

               tq(1)  =  1./2. + alp*cf*(3./4.*log(qm) - 3./2.)
               tg(1)  =          alp*nf*(3.*log(qm) - 9.)/16.

            elseif(iknl.eq.2) then

               tq(1)  =  1./2. + alp*cf*(3./4.*log(qm) - 2.)
               tg(1)  =          alp*nf*(log(qm) - 3.)/16.

            else
               write(6,*) 'iknl not rec. in wcread... stop'
               stop
            endif
         else
            write(6,*) 'n - order not rec. in wcread... stop'
            stop
         endif

         return
         end


C
C     interpolation and extrapolation routine based on rational functions
C

      SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=10,TINY=1.E-25,MXX=1050)
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

C
C     subroutine to compute counterterm of the coefficient functions
C     first unpolarized case
C

      subroutine counterV(tw,del, qone, mu, nf, reqint, imqint, regint, 
     &                    imgint,n)
      implicit none
      integer tw
      real*8 del, qone, mu, nf
      real*8 reqint, imqint, regint, imgint
      real*8 reqint1, imqint1, regint1, imgint1
      real*8 reqint2, imqint2, regint2, imgint2
      real*8 arg1,arg2,alp
      real*8 cf 
      real*8 Spence,alpdr
      external Spence,alpdr
      real*8 pi
      PARAMETER (pi = 3.14159265358979324D0)
      integer n,iswit

      cf = 4./3.
C
C     some arguments of logs
C
      arg1 = (1.-del)/del
      arg2 = qone**2/mu**2

C
C     first LO
C

      if (tw.eq.2) then

      reqint1 = -log(arg1)/2.

      imqint1 = pi/2. 

      regint1 = 0.0

      imgint1 = 0.0

      elseif (tw.eq.3) then

      reqint1 = pi**2/6. - Spence(del)/2. - log(del)**2/4

      imqint1 = -pi*log(del)/2.

      regint1 = 0.0

      imgint1 = 0.0

      endif

C
C     Now NLO!
C
C     Multiply the imaginary parts with the proper 
C     colorfactors cf,nf!
C

      reqint2 = cf/4. * ( pi**2/2. + log(arg2) * (pi**2 
     &- 3.*log(arg1) - (log(arg1))**2) - 3.*Spence(-arg1) 
     &+ log(arg1)*(3.*log(del) + pi**2 + 9.) -log(arg1)**3/3.) 

      imqint2 = -cf*pi/4.*(pi**2/3. + 9. + 3.*log(del) - log(arg1)**2 
     & - log(arg2)*(3. + 2.*log(arg1)))

      regint2 = nf/4.*(-1. + pi**2/4.*(4./3. - del) + (2.-del)
     & *log(arg1)*(1. - log(arg1)/4.) 
     & + log(arg2)/2.*(1. - (2.-del)*log(arg1)) 
     & + Spence(-arg1) - log(del)*log(arg1))

      imgint2 = -nf*pi/8.*((2.-del) * (2. - log(arg1) - log(arg2))
     &  - 2.*log(del))

C
C     value of alpha_s at q
C
      alp = alpdr(n,nf,mu)/2./pi

C
C     markers for LO and NLO
C

      if (n.eq.1) then

         iswit = 0

      elseif (n.eq.2) then

         iswit =1

      else
         write(6,*) 'n not one or two in counterV stop'
         stop
      endif   

C
C     assembly of counter terms, bothe for real and imaginary parts
C

      reqint = del*(reqint1 + dble(iswit) * alp * reqint2)

      imqint = del*(imqint1 + dble(iswit) * alp * imqint2)

      regint = del*(regint1 + dble(iswit) * alp * regint2)

      imgint = del*(imgint1 + dble(iswit) * alp * imgint2)

      return

      end

C
C     same subroutine as above but now for the polarized case
C

      subroutine counterA(tw,del, qone, mu, nf, reqinta, imqinta, 
     &                    reginta, imginta,n)
      implicit none
      integer tw
      real*8 del, qone, mu, nf
      real*8 reqinta, imqinta, reginta, imginta
      real*8 reqint1, imqint1, regint1, imgint1
      real*8 reqint2, imqint2, regint2, imgint2
      real*8 arg1,arg2,alp
      real*8 cf 
      real*8 Spence,alpdr
      external Spence,alpdr
      real*8 pi
      PARAMETER (pi = 3.141592653589793D0)
      integer n,iswit

      cf = 4./3.

      arg1 = (1.-del)/del
      arg2 = qone**2/mu**2


      if (tw.eq.2) then

      reqint1 = -log(arg1)/2.

      imqint1 = pi/2. 

      regint1 = 0.

      imgint1 = 0. 

      elseif(tw.eq.3) then

      reqint1 = pi**2/6. - Spence(del)/2. - log(del)**2/4

      imqint1 = -pi*log(del)/2.

      regint1 = 0.0

      imgint1 = 0.0

         endif

      reqint2 = cf/4.*( pi**2/6. + log(arg2)*
     & (pi**2 - 3.*log(arg1) - (log(arg1))**2) - Spence(-arg1)
     & + log(arg1)*(log(del) + pi**2 + 9.) -log(arg1)**3/3.) 

      imqint2 = -cf*pi/4.*(pi**2/3. + 9. + log(del) - log(arg1)**2 
     & - log(arg2)*(3. + 2.*log(arg1)))

      regint2 = nf/4.*(1. + pi**2*del/4. + del*log(arg1)*
     & (1.-log(arg1)/4.) - log(arg2)/2.*(1. + del*log(arg1)))

      imgint2 = -nf*pi*del/8.*(2. - log(arg1) - log(arg2))

      alp = alpdr(n,nf,mu)/2./pi

      if (n.eq.1) then

         iswit = 0

      elseif (n.eq.2) then

         iswit =1

      else
         write(6,*) 'n not one or two in counterV stop'
         stop
      endif   

      reqinta = del*(reqint1 + dble(iswit) * alp * reqint2)

      imqinta = del*(imqint1 + dble(iswit) * alp * imqint2)

      reginta = del*(regint1 + dble(iswit) * alp * regint2)

      imginta = del*(imgint1 + dble(iswit) * alp * imgint2)

      return

      end


C
C
C     fucntion to compute the spence function appearing in the counterterm 
C     functions!
C

      REAL*8 FUNCTION SPENCE(X)
      IMPLICIT NONE
      REAL*8 X,C,Z1,HF,PI,PI3,PI6,PI12
      REAL*8 H,T,S,A,Y,ALFA,B0,B1,B2
      INTEGER I
      DIMENSION C(0:19)

      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159265358979324D0)
      PARAMETER (PI3 = PI**2/3., PI6 = PI**2/6., PI12 = PI**2/12.)

      DATA C( 0) / 0.42996 69356 08136 97D0/
      DATA C( 1) / 0.40975 98753 30771 05D0/
      DATA C( 2) /-0.01858 84366 50145 92D0/
      DATA C( 3) / 0.00145 75108 40622 68D0/
      DATA C( 4) /-0.00014 30418 44423 40D0/
      DATA C( 5) / 0.00001 58841 55418 80D0/
      DATA C( 6) /-0.00000 19078 49593 87D0/
      DATA C( 7) / 0.00000 02419 51808 54D0/
      DATA C( 8) /-0.00000 00319 33412 74D0/
      DATA C( 9) / 0.00000 00043 45450 63D0/
      DATA C(10) /-0.00000 00006 05784 80D0/
      DATA C(11) / 0.00000 00000 86120 98D0/
      DATA C(12) /-0.00000 00000 12443 32D0/
      DATA C(13) / 0.00000 00000 01822 56D0/
      DATA C(14) /-0.00000 00000 00270 07D0/
      DATA C(15) / 0.00000 00000 00040 42D0/
      DATA C(16) /-0.00000 00000 00006 10D0/
      DATA C(17) / 0.00000 00000 00000 93D0/
      DATA C(18) /-0.00000 00000 00000 14D0/
      DATA C(19) /+0.00000 00000 00000 02D0/

      IF(X .EQ. 1) THEN
       H=PI6
      ELSEIF(X .EQ. -1) THEN
       H=-PI12
      ELSE
       T=-X
       IF(T .LE. -2) THEN
        Y=-1/(1+T)
        S=1
        A=-PI3+HF*(LOG(-T)**2-LOG(1+1/T)**2)
       ELSEIF(T .LT. -1) THEN
        Y=-1-T
        S=-1
        A=LOG(-T)
        A=-PI6+A*(A+LOG(1+1/T))
       ELSE IF(T .LE. -HF) THEN
        Y=-(1+T)/T
        S=1
        A=LOG(-T)
        A=-PI6+A*(-HF*A+LOG(1+T))
       ELSE IF(T .LT. 0) THEN
        Y=-T/(1+T)
        S=-1
        A=HF*LOG(1+T)**2
       ELSE IF(T .LE. 1) THEN
        Y=T
        S=1
        A=0
       ELSE
        Y=1/T
        S=-1
        A=PI6+HF*LOG(T)**2
       ENDIF
       H=Y+Y-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 19,0,-1

       B0=C(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=-(S*(B0-H*B2)+A)
      ENDIF
      SPENCE=H
      RETURN
      END

      FUNCTION SMPNOL (NX, DX, FN, ERR)
C                                            DP Left-Open Simpson Integration
C     Inputs:  Nx, Dx, Fn(1:Nx)                             [F(1) is not used]
C     Output:  Err                                          (Error estimate)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FN(NX)

      MS = MOD(NX, 2)
      IF (NX .LE. 1 .OR. NX .GT. 1000) THEN
         PRINT *, 'NX =', NX, ' OUT OF RANGE IN SMPNOL!'
         STOP
      ELSEIF (NX .EQ. 2) THEN
         TEM = DX * FN(2)
      ELSEIF (NX .EQ. 3) THEN
         TEM = DX * FN(2) * 2.
      ELSE
         IF (MS .EQ. 0) THEN
            TEM = DX * (23.* FN(2) - 16.* FN(3) + 5.* FN(4)) / 12.
            TMP = DX * (3.* FN(2) - FN(3)) / 2.
            ERR = ABS(TEM - TMP)
            TEM = TEM + SMPSNA (NX-1, DX, FN(2), ER1)
            ERR = ABS(ER1) + ERR
         ELSE
            TEM = DX * (8.* FN(2) - 4.* FN(3) + 8.* FN(4)) / 3.
            TMP = DX * (3.* FN(2) + 2.* FN(3) + 3.* FN(4)) / 2.
            ERR = ABS(TEM - TMP)
            TEM = TEM + SMPSNA (NX-4, DX, FN(5), ER1)
            ERR = ABS(ER1) + ERR
         ENDIF
      ENDIF

      SMPNOL = TEM
      RETURN
C                        ****************************
      END
C
      FUNCTION SMPNOR (NX, DX, FN, ERR)
C                                              DP Right-Open Simpson Integration
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FN(NX)

      MS = MOD(NX, 2)
      IF (NX .LE. 1 .OR. NX .GT. 1000) THEN
         PRINT *, 'NX =', NX, ' OUT OF RANGE IN SMPNOR!'
         STOP
      ELSEIF (NX .EQ. 2) THEN
         TEM = DX * FN(Nx-1)
      ELSEIF (NX .EQ. 3) THEN
         TEM = DX * FN(NX-1) * 2.
      ELSE
        IF (MS .EQ. 0) THEN
         TEM = DX * (23.* FN(NX-1) - 16.* FN(NX-2) + 5.* FN(NX-3)) / 12.
         TMP = DX * (3.* FN(NX-1) - FN(NX-2)) / 2.
         ERR = ABS(TEM - TMP)
         TEM = TEM + SMPSNA (NX-1, DX, FN(1), ER1)
         ERR = ER1 + ERR
        ELSE
         TEM = DX * (8.* FN(NX-1) - 4.* FN(NX-2) + 8.* FN(NX-3)) / 3.
         TMP = DX * (3.* FN(NX-1) + 2.* FN(NX-2) + 3.* FN(NX-3)) / 2.
         ERR = ABS(TEM - TMP)
         TEM = TEM + SMPSNA (NX-4, DX, FN(1), ER1)
         ERR = ER1 + ERR
        ENDIF
      ENDIF

      SMPNOR = TEM
      RETURN
C                        ****************************
      END

      FUNCTION SMPSNA (NX, DX, F, ERR)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MAXX = 1000)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DIMENSION F(NX)
      DATA IW1, IW2, TINY / 2*0, 1.E-35 /
C
      IF (DX .LE. 0.) THEN
c        CALL WARNR(IW2,NWRT,'DX cannot be < 0. in SMPSNA', 'DX', DX,
c     >               D0, D1, 0)
        SMPSNA = 0.
        RETURN
      ENDIF

      IF (NX .LE. 0 .OR. NX .GT. MAXX) THEN
c        CALL WARNI(IW1, NWRT, 'NX out of range in SMPSNA', 'NX', NX,
c     >               1, MAXX, 1)
        SIMP = 0.
      ELSEIF (NX .EQ. 1) THEN
        SIMP = 0.
      ELSEIF (NX .EQ. 2) THEN
        SIMP = (F(1) + F(2)) / 2.
        ERRD = (F(1) - F(2)) / 2.
      ELSE
        MS = MOD(NX, 2)

C For odd # of intervels, compute the diff between the Simpson 3/8-rule result
C for the last three bins and the regular Simpson rule result for the next to
C last two bins.  The problem is thereby reduced to a even-bin one in all cases.

        IF (MS .EQ. 0) THEN
          ADD = (9.*F(NX) + 19.*F(NX-1) - 5.*F(NX-2) + F(NX-3)) / 24.
          NZ = NX - 1
        ELSE
          ADD = 0.
          NZ = NX
        ENDIF

        IF (NZ .EQ. 3) THEN
          SIMP = (F(1) + 4.* F(2) + F(3)) / 3.
          TRPZ = (F(1) + 2.* F(2) + F(3)) / 2.
        ELSE
          SE = F(2)
          SO = 0
          NM1 = NZ - 1
          DO 60 I = 4, NM1, 2
            IM1 = I - 1
            SE = SE + F(I)
            SO = SO + F(IM1)
c	  print *, F(I), F(IM1)
   60     CONTINUE
          SIMP = (F(1) + 4.* SE + 2.* SO + F(NZ)) / 3.
          TRPZ = (F(1) + 2.* (SE + SO) + F(NZ)) / 2.
        ENDIF

        ERRD = TRPZ - SIMP
        SIMP = SIMP + ADD

      ENDIF
C
      SMPSNA = SIMP * DX

      IF (ABS(SIMP) .GT. TINY) THEN
        ERR = ERRD / SIMP
      ELSE
        ERR = 0.
      ENDIF
C
      RETURN
C                        ****************************
      END

C
C     function that computes the jacobian in going from a non-equidistant 
C     grid spacing to an equidistant one, necessary for simpson integration!
C     functional from of transform is different depending on whether one is
C     in the ERBL or DGLAP region!
C

      function dxdzn(i,n,nx,x,del)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      dim = 1./DEL

      if (i.eq.1) then

      IF (DEL.LT.0.1) THEN
      MS = MOD(NX,4)

      IF (MS.EQ.0) THEN

      NC = NX/8

      ELSE

      NC = (NX/4+1)/2

      ENDIF
      a = 10.
      s = 4.*((1.*(NX-1))/(1.*(NX-4)))
      ELSE
      MS = MOD(NX,2)

      IF (MS.EQ.0) THEN

      NC = NX/4

      ELSE

      NC = (NX/2+1)/2

      ENDIF
      a = 7.
      s = 2.*((1.*(NX-1))/(1.*(NX-2)))
      ENDIF

      xp = x*s

      a1 = a*(xp-1./2.)**2
      b = dexp(-a/4.)
      b1 = dexp(-a1)

      if (n.LE.NC) then

      tem3 = DEL*(1.-2.*xp)*s*(a*b1 + 4.*b)/2.

      ELSE

      tem3 = - DEL*(1.-2.*xp)*s*(a*b1 + 4.*b)/2.

      ENDIF

      if (abs(tem3).lt.1E-10) tem3 = 0.0

      dxdzn = tem3

      else

      b = - LOG(dim)

      if (del.lt.0.1) then

      a = 1./4.*(NX-4)/(1.*(NX-1))
      a1 = LOG(a)

      else

      a = 1./2.*(NX-2)/(1.*(NX-1))   
      a1 = LOG(a)

      endif

      c = b/a1
      
      dxdzn = c*x**(c-1.)

      endif

      return     
c
      end
C



