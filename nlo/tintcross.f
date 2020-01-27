      subroutine tintcross
      implicit none

      integer ordint
      real    xint, q2int
      common / varint / ordint, xint, q2int
      
      real tmin, tmax
      common / tvar / tmin, tmax


      character name*19
      real     res

      double precision funct, funct1, D01AHF
      external         funct, funct1, D01AHF

      double precision funct2, phi1, phi2, D01DAF 
      external         funct2, phi1, phi2, D01DAF 

      double precision YA,YB, ANS, ABSACC

      double precision  A, B, EPSR, RELERR
      integer NPTS, NLIMIT, IFAIL

      real s, fact, factgp, pi, y, alpha, w
      integer nbins

      real jacob

      double precision w_min, w_max, x_min, x_max


      double precision result

      real vect_w_lo(100),  vect_w_nlo(100),
     &     vect_q2_lo(100), vect_q2_nlo(100)


      real vect_q2_ep(100)

      real mp

      integer ic

      
      A = -1.0


      EPSR = 0.01
      NLIMIT = 10000

      s = 90200.
      mp = 0.938272

      fact = 0.38937966*1000000.
      pi    = 3.1415
      alpha = 1./137.

      nbins = 100

* calculate single differential cross sections d\sigma/dQ2
* which means to integrate over t and x

      w_min =  30.
      w_max = 120.

      do ic = 1, nbins

         q2int  = 4.5+real(ic)*1.5/real(nbins)
         ordint = 1

         YA = dble(q2int)/(w_max**2)
         YB = dble(q2int)/(w_min**2)

c         YA = 3.0
c         YB = 7.0


         tmin = -1.0D0
c         tmax = dble(-1.* xint**2*mp**2/(1.-xint-xint*mp**2/q2int))

         ABSACC = 0.01

         result = D01DAF(YA,YB,PHI1,PHI2,funct2,ABSACC,ANS,NPTS,IFAIL)

         vect_q2_ep(ic) = real(ans)*fact

c         print*,ic,vect_q2_ep(ic),npts,ifail

         print*,ya,yb,result

      enddo

      goto 20

*     cross section as a function of W

      A = -1.0

      nbins = 100

      do ic = 1, nbins

*
*** vect_w_lo
*
         w = 30.+ (real(ic)*90./real(nbins))
         q2int  = 4.5
         xint   = q2int/(w**2)
         y =q2int/(s*xint)
crs         B = dble(-1.* xint**2*mp**2/(1.-xint-xint*mp**2/q2int))
         B = -0.01
         factgp = (alpha/pi) * ((1.+(1.-y)**2)/(2.*y*q2int))
         jacob = q2int/((xint**2)*s)

         ordint = 1
         ifail = 0

         result = D01AHF(A,B,EPSR,NPTS,RELERR,funct,NLIMIT,IFAIL)
         res = real(result)*fact/(factgp*jacob)
c         res = real(result)*fact

         vect_w_lo(ic) = res

*
*** vect_w_nlo
*
         ordint = 2
         ifail = 0

         result = D01AHF(A,B,EPSR,NPTS,RELERR,funct,NLIMIT,IFAIL)

         res = real(result)*fact/(factgp*jacob)
c         res = real(result)*fact

         vect_w_nlo(ic) = res

*
*** vect_q2_lo
*
         w = 75.
         q2int  = 2.0+real(ic)*17.9/real(nbins)
         xint   = q2int/(w**2)
         y =q2int/(s*xint)
         factgp = (alpha/pi) * ((1.+(1.-y)**2)/(2.*y*q2int))
         jacob = q2int/((xint**2)*s)

crs         B = dble(-1.* xint**2*mp**2/(1.-xint-xint*mp**2/q2int))
         B = -0.01

         ordint = 1
         ifail = 0

         result = D01AHF(A,B,EPSR,NPTS,RELERR,funct,NLIMIT,IFAIL)
         res = real(result)*fact/(factgp*jacob)
c         res = real(result)*fact

         vect_q2_lo(ic) = res

*
*** vect_q2_nlo
*

         ordint = 2
         ifail = 0

         result = D01AHF(A,B,EPSR,NPTS,RELERR,funct,NLIMIT,IFAIL)
         res = real(result)*fact/(factgp*jacob)
c         res = real(result)*fact

         vect_q2_nlo(ic) = res
       
      enddo

*
** for tests
*


c         q2int  = 4.0
c         xint   = 0.0003
c         y =q2int/(s*xint)
c
c         ordint = 1
c         ifail = 0
c
c         result = D01AHF(A,B,EPSR,NPTS,RELERR,funct,NLIMIT,IFAIL)
c         res = real(result)*fact
c
c         print*,res


*
*
*


 20   continue


 101  FORMAT(F16.8)

      name  = 'txt/vect_q2_ep.txt'
      print*,name
      
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic=1,nbins
         write(95,101) vect_q2_ep(ic)
      enddo
      CLOSE(95)



      name  = 'txt/vect_w_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic=1,nbins
         write(95,101) vect_w_lo(ic)
      enddo
      CLOSE(95)

      name  = 'txt/vect_w_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic=1,nbins
         write(95,101) vect_w_nlo(ic)
      enddo
      CLOSE(95)

      name  = 'txt/vect_q2_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic=1,nbins
         write(95,101) vect_q2_lo(ic)
      enddo
      CLOSE(95)

      name  = 'txt/vect_q2_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic=1,nbins
         write(95,101) vect_q2_nlo(ic)
      enddo
      CLOSE(95)


      return
      end


