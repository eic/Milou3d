      double precision function phi2(Y)
      implicit none
      double precision Y

      real tmin, tmax
      common / tvar / tmin, tmax

      integer ordint
      real    xint, q2int
      common / varint / ordint, xint, q2int
 
      double precision mp, dq2

      mp = 0.938D0
      dq2 = dble(q2int)


c      phi2 = dble(tmax)

c       phi2 = -1.* Y**2*mp**2/(1.-Y-Y*mp**2/dq2)

c      phi2 = -1.*10**(-10.)

c       print*,Y,dq2,phi2


      phi2 = 0.

      return
      end


