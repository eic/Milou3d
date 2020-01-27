      subroutine diffcrossdvcs (ordi,t,x,q,result1)
C 
C     Definition block 
C 
 
      implicit none 
 
C 
C     arrays and common blocks for del,q and x 
C 
       
      real*8 q,t,x,result1 
      integer ordi


C 
C     Kinematical constants and parameters used throughout the program 
C      

      external forms

      real*8 S,pit
      real*8 dM,tmax,tmin,eps 
      real*8 phimin,phimax 
      real*8 ymax,ysim 
      real*8 y,testt,del,qsq 
      common /req/ y,tmin,testt,del,qsq 
      parameter(pit = 3.141592653589793)
 
C 
C     parameters for parametrization of t-dependence (form factors) 
C 
      real*8 al,kp,kn,ksea,mp,mn,mv,ma,mpi,ga3 
      data al,kp,kn,ksea,mp,mn,mv,ma,mpi,ga3/7.29735D-3, 
     &        1.79,-1.91,0.0,0.938272,0.939565, 
     &        0.840,0.900,0.134977,1.267/ 
       
C 
C     initiation of phi integration for diff. cross section 
C 
      real*8 corint1,corcosi,cor1,corcos,result
      external cor1,corcos
C 
C     flag to indicate LO or NLO 
C 
       
      integer iyflag 
c      integer ordn
c      common /order/ ordi

C 
C     Definition of the various contribution to the diff. Xsection from DVCS,  
C     interference and BH 
C       
      real*8 tbhup,tdvcsup,tiup 
      external tbhup,tdvcsup,tiup
C 
C     set up of integration routine 
C 
      external D01AKF 
      INTEGER LW,IW,LIW,IFAIL 
      PARAMETER(LW=4000,LIW=LW/4) 
      REAL*8 epsabs,epsrel 
      real*8 abserr,w 
      DIMENSION w(lw),iw(liw) 
C 
C     set up for expressions for diff. Xsection and asymmetries 
C 
      
      real*8 difxuw 

 
C     initiate error boundaries for integral routine 
 
      EPSABS = 0.0 
      EPSREL =1.0D-5 
      IFAIL =-1 

         qsq = q

         testt = t

         del = x
   
C 
C     max y und t 
C 
      tmax = -1.0D0 
      ymax = 0.95D0 
      eps = 1.0D0 

C 
C     Calculating parameters depending on del and q! 
C 
 
crs      S = 4.*27.5*920. 
      S = 4.*27.5*820. 
 
      iyflag = 0 
 
      tmin = -del**2 * mp**2 / (1. -del - del*mp**2/qsq) 
      y = qsq/(del*S) * (1.-del*(1. -mp**2/qsq)) 

c      print*,del,mp,qsq,y
  
C 
C     y warning if y out of ex. range. Skip the calculation for this Q^2 value. 
C  
            if (y.gt.ymax) then 
               iyflag = iyflag + 1 
c      write(6,*) 'y = ',y, ' out of natural range for x,qsq',del,qsq 
c               y = ymax 
c               write(6,*) 'For convenience set y = ymax = 0.95 but' 
c               write(6,*) 'WARNING: Results unreliable for HERA ENERGY' 
c               write(6,*) ' ' 
               goto 200
            endif 
 
C 
C     compute diff cross section unpol., pol. ,transversely pol. 
c 
 
C 
C     phi integrals 
C 
      phimin = 0.0 
      phimax = 2.0*3.141592653589793 
      call D01AKF(cor1,phimin,phimax,epsabs,epsrel,result,abserr,w, 
     &            lw,iw,liw,ifail) 
      corint1 = result 
      
      call D01AKF(corcos,phimin,phimax,epsabs,epsrel,result,abserr,w, 
     &            lw,iw,liw,ifail) 
      corcosi = result 

c      print*,corint1,corcosi

  
C 
C     compute dM from phas space 
C 

      dM = al**3*del*y**2/(8*pit*qsq**2* 
     >   (1 + 4*del**2*mp**2/qsq)**0.5) 
C 
C     diff. Xsection 
C 

c      corcosi = 0.0
c      corint1 = 0.0


      difxuw =(2*pit*tdvcsup(ordi,testt)+corcosi*tiup(ordi,testt) 
     > + corint1*tbhup(ordi,testt))*dM

c      call forms(t) 
c
c      difxuw =(corint1*tbhup(ordi,testt))*dM

      result1 = difxuw



c      print*,pit,tdvcsup(ordi,testt),corcosi,tiup(ordi,testt),
c     &       corint1,tbhup(ordi,testt),dM,result1

 
 200  return 
      end 
 

  

C 
C     functions to compute the various contributions of DVCS, interference and BH 
C     to the UNP, LP and TP Xsections, separated by angular dependence. Also used to form  
C     asymmetries 
C 
      real*8 function tdvcsup(ordi,t) 
C 
C     definition block 
C 
      implicit none 
       
      real*8 t 
      real*8 h1re,h1im,h1tre,h1tim 
      real*8 e1re,e1im,e1tre,e1tim 
      common /amps1/h1re,h1im,e1re,e1im,h1tre,h1tim,e1tre,e1tim 
      real*8 phicap 
      common /CAP/ phicap 
      real*8 fact 
      real*8 mp 
      data mp/0.938272/ 
      real*8 pi 
      parameter(pi = 3.141592653589793) 
      integer ordi 
 
      real*8 y,tmin,testt,del,qsq 
      common /req/ y,tmin,testt,del,qsq 
      real*8 tbhup,tiup 
 
      real*8 f1p,f2p,f1n,f2n 
      common /form2/ f1p,f2p,f1n,f2n 
 
      external ampli 
C
C  unpol. DVCS^2
C 
      
      call ampli(ordi,t,del,qsq)
  
      tdvcsup = 2.*(2.-2.*y + y**2)/(y**2*(2.-del)**2 * qsq) * 
     &  (4.*(1.-del)*(h1re**2 + h1im**2 + h1tre**2 + h1tim**2)  
     &  - del**2*(2.*(h1re*e1re +h1im*e1im +h1tre*e1tre +h1tim*e1tim))  
     &  - (del**2 + (2.-del)**2 * t/4./mp**2) * (e1re**2 + e1im**2) 
     &  - del**2 * t/4./mp**2 *(e1tre**2 + e1tim**2)) 
 
      return 
 
C
C     unpol. BH^2
C

      entry tbhup(ordi,t) 
 
      tbhup = -2.*(2.-2.*y + y**2)/(y**2*t) *( 
     &   4.*(1.-del)/del**2 *(1.-tmin/t) * f1p**2 +2.*(f1p+f2p)**2 
     &   + (t/tmin -1.) *f2p**2) 
 
      return 
 
      entry tiup(ordi,t) 
 
C 
C     unpol. interference term
C 
  
      call ampli(ordi,t,del,qsq) 
 
      fact = 8.*((1.-y)*(1.-del)*(1.-tmin/t))**(0.5)/ 
     &       (del*((-t*qsq))**(0.5)) 

      tiup = -(2.-2.*y+y**2) * (fact/y**3) * 
     & (f1p*h1re + del/(2.-del)*  
     & (f1p + f2p) * h1tre - t/4./mp**2 * f2p * e1re) 

      return 
  
      end 
 
C 
C     subroutine which computes the amplitudes h1 etc. 
C 

C
C     you need to write the functions with the LO and NLO 
C     parametrizations of the amplitudes, which are then read in below.
C

C
C     you can have an alternative t-depndence than the one I chose here
C     like a a global exponential like exg in forms!
C
 
      subroutine ampli(ord,t,x,q) 
 
      real*8 t,x,q,i,j 
      integer ord 
      external forms 
      real*8 f1u,f2u,f1d,f2d 
      real*8 f1s,f2s 
      real*8 g1,g1sea,gpi,exg 

      common /form1/ f1u,f1d,f1s,exg,f2u,f2d,f2s,g1,g1sea,gpi 
      real*8 reu,imu,red,imd,res,ims,reg,img,reup,imup, 
     > redp,imdp,resp,imsp,regp,imgp,reue,imue,rede,imde,rese,imse
     > ,rege,imge,reuep,imuep,redep,imdep,resep,imsep,regep,imgep     
      
      real*8 h1re,h1im,e1re,e1im,h1tre,h1tim,e1tre,e1tim 
      common /amps1/h1re,h1im,e1re,e1im,h1tre,h1tim,e1tre,e1tim 
       
      external reu,imu,red,imd,res,ims,reg,img,reup,imup, 
     > redp,imdp,resp,imsp,regp,imgp,reue,imue,rede,imde,rese,imse
     > ,rege,imge,reuep,imuep,redep,imdep,resep,imsep,regep,imgep

 
      double precision resu1,resu2
      
      if (ord.eq.1) then 
 
         sw = 0.0 
 
         else 
 
            sw = 1D0 
 
       endif 

       i = x

       j = q

       call forms(t) 

       ord = 1

c       resu1 = regep(ord,i,j)
c       resu2 = imgep(ord,i,j)
c       print*,ord,i,j,resu1,resu2
c       print*,f1u,f1d,f1s,exg,f2u,f2d,f2s,g1,g1sea,gpi,sw 
c       return

crs      h1re= reu(i,j)*f1u+red(ord,i,j)*f1d+res(ord,i,j)*f1s
      h1re= reu(ord,i,j)*f1u+red(ord,i,j)*f1d+res(ord,i,j)*f1s
     > +sw*reg(ord,i,j)*exg  
      h1im=imu(ord,i,j)*f1u+imd(ord,i,j)*f1d+ims(ord,i,j)*f1s
     > +sw*img(ord,i,j)*exg  
      e1re=(reue(ord,i,j)*f2u+rede(ord,i,j)*f2d+rese(ord,i,j)*f2s
     > +sw*rege(ord,i,j)*exg)  
      e1im=(imue(ord,i,j)*f2u+imde(ord,i,j)*f2d+imse(ord,i,j)
     > *f2s+sw*imge(ord,i,j)*exg)  
      h1tre=reup(ord,i,j)*g1+redp(ord,i,j)*g1+resp(ord,i,j)*g1sea
     > +sw*regp(ord,i,j)*exg  
      h1tim=imup(ord,i,j)*g1+imdp(ord,i,j)*g1+imsp(ord,i,j)*g1sea
     > +sw*imgp(ord,i,j)*exg  
      e1tre=reuep(ord,i,j)*gpi+redep(ord,i,j)*gpi
     > +resep(ord,i,j)*gpi+sw*regep(ord,i,j)*exg  
      e1tim=imuep(ord,i,j)*gpi+imdep(ord,i,j)*gpi
     > +imsep(ord,i,j)*gpi+sw*imgep(ord,i,j)*exg 
 

      return 
 
      end 
 
C 
C     functions that store the different angular depnedences of the various DVCS,  
C     interference and BH terms 
C 

 
 
      real*8 function corcos(phi) 
      real*8 phi,cor1,arg 
      external cor1 
      real*8 y,tmin,testt,del,qsq 
      common /req/ y,tmin,testt,del,qsq 
      arg = phi 
      corcos = DCOS(arg)*cor1(phi) 
      return 
      end 
 
 
C 
C     function that computes the value of the product of the two BH propagators! 
C     now for the unweigted cross section! 
C 
      real*8 function cor1(phi) 
      real*8 phi,brac,twophi,t 
      real*8 y,tmin,testt,del,qsq 
      real*8 c1,c2,c3 
      common /req/ y,tmin,testt,del,qsq 
      t = testt 
      twophi = 2.*phi 
      c0 = -(t-tmin)*(1.-del)*(1.-y)/qsq 
      if(c0.lt.0.0) then 
         write(6,*) 'c0 negative...for' 
         write(6,*) 'phi,twophi,y,tmin,t,del,qsq' 
         write(6,*)  phi,twophi,y,tmin,t,del,qsq 
         write(6,*) 'STOP' 
         stop 
      else    
      c1 = 2.*(c0)**0.5 
      c2 = -t/qsq*(1.-del*(2.-y))  
      c3 = -t/qsq*(1.-y-del*(2.-y)) 
      brac = (1.-y+c1*dcos(phi)+c2)*(1.+c1*dcos(phi)+c3)  
c      print *,"cor",1.-y,c1,-t*c2/qsq,-t*c3*qsq 
c      read(*,*) 
      cor1 = y**2/brac 
      endif 
      return 
      end 
 
 
C 
C 
C     function that sets the t dependence of the various amplitudes 
 
      subroutine forms (t) 
      implicit none 
      real*8 t 
      real*8 kp,kn,ksea,mp,mn,mv,ma,mpi,ga3 
      data kp,kn,ksea,mp,mn,mv,ma,mpi,ga3/ 
     &       1.79,-1.91,-2.0,0.938272,0.939565, 
     &       0.840,0.900,0.134977,1.267/ 
      real*8 gep,gmp,gmn,gen 
      real*8 f1u,f2u,f1d,f2d 
      real*8 ges,gms,f1s,f2s 
      real*8 g1,g1sea,gpi,exg 
      common /form1/ f1u,f1d,f1s,exg,f2u,f2d,f2s,g1,g1sea,gpi 
      real*8 f1p,f2p,f1n,f2n 
      common /form2/ f1p,f2p,f1n,f2n 
      gep = 1./(1.-t/mv**2)**2    
      gmp = gep * (1. + kp)    
      gmn = gep * kn          
      gen = 0.0               
      ges = 1./(1.-t/mv**2)**3 
      gms = (1. + ksea) * ges  
      f1p = (gep + t*gmp/(4.*mp**2))/(1.+t/(4.*mp**2)) 
      f2p = (gmp - gep)/(1.+t/(4.*mp**2)) 
      f1n = t*gmn/(4.*mn**2)/(1.+t/(4.*mn**2)) 
      f2n = gmn/(1.+t/(4.*mn**2)) 
      f1u = (2.*f1p + f1n)/2. 
      f2u = (2.*f2p + f2n)/2. 
      f1d = f1p + 2.*f1n 
      f2d = f2p + 2.*f2n 
      f1s = (ges + t/(4.*mv**2)*gms)/(1.+t/(4.*mv**2)) 
      f2s = (gms - ges)/(1.+t/(4.*mv**2)) 
      g1   = 1./(1.-t/ma**2)**2         
      g1sea = 1./(1.-t/ma**2)**3         
      gpi = 4 * ga3 * mp**2/(mpi**2 -t)  
      exg = 1./(1.-t/ma**2)**3 
c      exg = 4.5*Exp(4.5*t) 
      return 
      end 
 
 

 
 
 
 
 
 
 
 
 
 
 
 
 
