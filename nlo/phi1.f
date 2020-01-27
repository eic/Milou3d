      double precision function phi1(Y)
      implicit none
      double precision Y

      real tmin, tmax
      common / tvar / tmin, tmax


      phi1 = dble(tmin)


      return
      end
