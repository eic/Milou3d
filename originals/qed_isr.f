*CMZ :          25/01/2004  10.41.30  by  H1 Saclay
*-- Author :    Unknown   17/12/2003


***********************************************************************
      subroutine qed_isr(SEP,Ee,Egam,weight)
C     ==================================

C 23.01.04
C changed such that Egamma is an input variable

***********************************************************************
*                                                                     *
*  qed_isr : this routine makes a photon radiation on an electron     *
*            (or positron) line in the collinear approximation        *
*                                                                     *
*  Author: Laurent Favart                                  03/11/99   *
*                                                                     *
*  Last updated:                                           03/11/99   *
*                                                                     *
***********************************************************************
*                                                                     *
* INPUT  Ee:   incoming electron energy                               *
*                                                                     *
* OUTPUT Ee:     electron energy after possible photon radiation      *
*        Egam:   radiated photon energy (=0 if no photon radiated)    *
*        weight: weight of the event. The cross section has to be     *
*                multiplied by this weight.                           *
*                                                                     *
* Remarks: This routine can only be used once per electron line. Two  *
*   photons radiation from the same line would imply interferences    *
*   which are not included.                                           *
*                                                                     *
*          The collinear approximation means that the photon is always*
*   emmitted along the electron direction and that the transverse     *
*   momentum of the electron does not change. This approximation only *
*   make sense if the ISR and FSR can be separeted i.e. if the angle  *
*   between the incoming ands the outgoing electron is large. It is   *
*   the case in DIS but not in photoproduction at HERA.               *
*                                                                     *
*          To simulate the photon, you have to decide above which     *
*   energy (Egam) you want to simulate it and include it in the final *
*   state. In most of the cases it is a waste of time to simulate     *
*   very small energy photon.                                         *
*                                                                     *
***********************************************************************

      implicit real(a-h,o-z)

c Pi (API) and the electron mass (AME)
      real API/3.14159/,AME/0.0005/

      real zmin

      zmin = 1e-4

      radcor=1.
      z=1.
      weight=1.

      zk = Egam/Ee

      if (zk.gt.zmin) then
       weight = (1./137.) / (2.*API)
       weight = weight * ( 1. + (1-zk)**2) / zk
       weight = weight * log(sep/(4.*AME**2))
      else
       aa = 3. +4.*log(zmin)-4.*zmin +zmin**2
       bb = ( (1./137.) / (2.*API) )* log(sep/(4.*AME**2))
       cc = 2. + bb*aa
       weight = cc / (2. *zmin)
      endif

      weight = weight / Ee

      z=1.-zk

      Ee=Ee*z

  999 return
      end
