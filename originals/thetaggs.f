      real function thetaggs
      include ?
c      thetaggs = 
c     +  abs(atan2(sqrt(prgbel(1)**2+prgbel(2)**2),prgbel(3))-
c     +      atan2(sqrt(pvgbel(1)**2+pvgbel(2)**2),pvgbel(3)))

      thetaggs =
     +  abs(atan2(sqrt(prgprf(1)**2+prgprf(2)**2),prgprf(3))-
     +      atan2(sqrt(pvgprf(1)**2+pvgprf(2)**2),pvgprf(3)))

      return
      end
