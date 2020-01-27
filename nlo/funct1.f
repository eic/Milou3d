      double precision function funct1(x)
      implicit none
*
* routine calculates the cross section integrated over t and 
* serves itself as function to be integrated over x
*

      double precision x
      double precision result

* variables needed to integrate over t

      double precision A, B, EPSR, RELERR
      integer NPTS,NLIMIT,IFAIL

      double precision funct, D01AHF 
      external         funct, D01AHF

* variables which are needed inside the function funct 
* when integrating over t

      integer ordint
      real    xint, q2int
      common / varint / ordint, xint, q2int

* variables set for the integration over t

      A = -1.0
      B = -0.1
      EPSR = 0.01
      NLIMIT = 10000


      xint = real(x)

      print*,x

c      q2int  !  these variables are set in the routine before
c      ordint !  tintcross

      result = D01AHF(A,B,EPSR,NPTS,RELERR,funct,NLIMIT,IFAIL)

      if (ifail.ne.0) then
         print*,'funct1: something wrong ',ifail
      endif


      funct1 = result

      print*,'result of integration over t ',result


      return
      end
