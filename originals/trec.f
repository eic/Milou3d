C------------------------------------------------
      REAl FUNCTION trec
C------------------------------------------------

      implicit real(a-h,o-z)
     
      integer num

      include ?
      
      real q22t, w2t, x2t,  y2t, t2t
      real S
      real epin, elepin
      real ergam,esclep,thrgam,thsclep,phrgam,phsclep,pzsclep,pzrgam

      call smear(sm1,0.03)
      call smear(sm2,0.03)

      call smear(sm3,1.e-3)
      call smear(sm4,1.e-3)

      call smear(sm5,3.e-3)
      call smear(sm6,3.e-3)

      xee_temp = plolab(4)
      xeg_temp = prglab(4)
      xthe_temp = atan2(sqrt(plolab(1)**2+plolab(2)**2),plolab(3))
      xphe_temp = atan2(plolab(2),plolab(1))
      xthg_temp = atan2(sqrt(prglab(1)**2+prglab(2)**2),prglab(3))
      xphg_temp = atan2(prglab(2),prglab(1))

      xee_temp = xee_temp*(1+sm1)
      xeg_temp = xeg_temp*(1+sm2)
      xthe_temp = xthe_temp*(1+sm3/xthe_temp)
      xthg_temp = xthg_temp*(1+sm4/xthg_temp)
      xphe_temp = xphe_temp*(1+sm5/xphe_temp)
      xphg_temp = xphg_temp*(1+sm6/xphg_temp)




      px_spa = xee_temp*sin(xthe_temp)*cos(xphe_temp)
      py_spa = xee_temp*sin(xthe_temp)*sin(xphe_temp)
      pz_spa = xee_temp*cos(xthe_temp)
      e_spa =  xee_temp

      px_lar = xeg_temp*sin(xthg_temp)*cos(xphg_temp)
      py_lar = xeg_temp*sin(xthg_temp)*sin(xphg_temp)
      pz_lar = xeg_temp*cos(xthg_temp)
      e_lar =  xeg_temp
 

      epin   = 920.
      elepin = 27.55

      S = 4*epin*elepin

      ergam  =e_lar
      esclep =e_spa
      thrgam  =atan2(sqrt(px_lar**2+py_lar**2),pz_lar)
      thsclep =atan2(sqrt(px_spa**2+py_spa**2),pz_spa)
      phrgam  =atan2(py_lar,px_lar)
      phsclep =atan2(py_spa,px_spa)
      pzsclep = esclep*cos(thsclep)
      pzrgam = ergam*cos(thrgam)
cc      print *,pzsclep,pzrgam

      q22t =4*(elepin**2)*sin(thrgam)*(1.+cos(thsclep))/
     &        (sin(thrgam)+sin(thsclep)-sin(thrgam+thsclep))

      x2t  = (elepin/epin)*((sin(thrgam)+sin(thsclep)+
     &        sin(thrgam+thsclep))/
     &        (sin(thrgam)+sin(thsclep)-sin(thrgam+thsclep)))

      y2t  = Q22t/(x2t*S)
      if (y2t.gt.0) then
        W2t    = sqrt(y2t*S)
      else
         W2t    = 290.
      endif

      t2t =   (ergam*sin(thrgam)*cos(phrgam)
     &      +  esclep*sin(thsclep)*cos(phsclep))**2
     &      + (ergam*sin(thrgam)*sin(phrgam)
     &      +  esclep*sin(thsclep)*sin(phsclep))**2

      trec=t2t 

      return
      end
CC-----------------------------
      SUBROUTINE SMEAR(SM1,SIG)
CC-----------------------------

      IMPLICIT NONE

      REAL SM1, SIG, RNDM
      INTEGER I

      DO I = 1, 12
         SM1 = SM1 + RNDM(1.)
      ENDDO

      SM1 = SIG*(SM1-6.0)

      RETURN
      END
