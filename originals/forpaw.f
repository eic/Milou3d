*CMZ :          27/05/2004  16.27.44  by  H1 Saclay
*-- Author :    Unknown   17/12/2003



*--------------------------------
      subroutine ntinit
*--------------------------------

      include 'forpaw.common'

* Initialization of ntuple and histos

      REAL hmemor
      PARAMETER( NWPAWC=500000 )
      COMMON/pawc/hmemor(NWPAWC)

* ----------------------------------------------------------


*---Create the Ntuple:
      IDNT = 1
      CALL HBNT(IDNT,'essai',' ')

*---Define the ntuple :
      NTID = 1

      call hbname(NTID,'dvcsvar1',xntp,'xntp,'//
     + 'qntp,yntp,phintp,tntp,'//
     + 'plilab(4),plolab(4),ppilab(4),ppolab(4),'//
     + 'pvglab(4),prglab(4),pvisr(4),'//
     + 'pliprf(4),ploprf(4),ppiprf(4),ppoprf(4),'//
     + 'pvgprf(4),prgprf(4),'//
     + 'plibel(4),plobel(4),ppibel(4),ppobel(4),'//
     + 'pvgbel(4),prgbel(4),'//
     + 'wbhntp,wdvcsntp,wintntp,wtotntp,'//
     + 'ym,phibel,resol,phibelrec,phibelgen')

      return
      end

*---------------------------------------------
      subroutine ntend
*---------------------------------------------

* close the file and write out ntuple and histo :

*      call hrout(1,icycle,' ')
*      call hrout(500,icycle,' ')

      call hrout(0,icycle,' ')

      call hrend('toto')

      LUNNT = 31

      close(LUNNT)

      return
      end

