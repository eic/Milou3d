*CMZ :          30/05/2004  08.59.08  by  H1 Saclay
*-- Author :    Unknown   17/12/2003


****************************************************************
*
      subroutine calc_kinem2(ibad,ioutkin,isca)
*
****************************************************************

      IMPLICIT REAL*8(A-H,O-Z)

      include 'dvcs.common'
      include 'forpaw.common'

      integer ibad,ioutkin

      integer debug
      common/debug/idebug

      integer iprodsca
      common /iprodsca/ iprodsca

      real elepi,EGAMR,sweight,srad,elepin,ehadi
      common /RADGEN/ elepi,EGAMR,sweight,srad,elepin,ehadi

      common/dvcs_VAR/x_main,q_main,phi_main,t_main,ym_main


      real rndm,pi_real,th_lo,en_lo,th_rg,en_rg,en_su,w_eg
      real*8 mp

* 4-vectors : (px, py,pz; E)

* -> Lab frame :
      dimension pli_l(4)   ! incoming lepton in lab frame
      dimension plo_l(4)   ! outcoming lepton in lab frame
      dimension pvgam_l(4) ! virtual photon in lab frame
      dimension ppi_l(4)   ! incoming proton in lab frame
      dimension ppo_l(4)   ! outcoming proton
      dimension prgam_l(4) ! outcoming photon

* -> Proton rest frame :
      dimension pli_tmp(4)   ! incoming lepton in p-rest frame
      dimension plo_tmp(4)   ! outcoming lepton in p-rest frame
      dimension pvgam_tmp(4) ! virtual photon in p-rest frame
      dimension ppi_tmp(4)   ! incoming proton in p-rest frame
      dimension ppo_tmp(4)   ! outcoming proton
      dimension prgam_tmp(4) ! outcoming photon
      dimension toto_li(4)
      dimension toto_lo(4)
      dimension toto_vgam(4)

* -> In p-rest frame, with z = - q_1^3
*    (cf fig. 1 of Belitsky et al.)
      dimension pli_b(4)   ! incoming lepton
      dimension plo_b(4)   ! outcoming lepton
      dimension pvgam_b(4) ! virtual photon
      dimension ppi_b(4)   ! incoming proton
      dimension prgam_b(4) ! outcoming (real) photon
      dimension ppo_b(4)   ! outcoming proton


      dimension beta(3),betar(3)
      dimension ptmp(4)
      dimension rot(3)

      logical do_xaxis
      common/xaxis/do_xaxis,i_xaxis

      ibad     = 0
      ioutkin  = 0
      iprodsca = 0

      mp = 0.93827d0
      pi = datan(1.d0) * 4.0d0

      elep    = dble(elepi)
      eproton = dble(stETARG)
      shat    = 4.0d0 * elep * eproton

cls      print *,elep,eproton,sweight

      q2  = qntp
      xb  = xntp
      t   = tntp
      phi = phintp

cep      dym = ym
      dym = ym_main


      y = q2/(xb*(shat - mp**2))

      yntp = y

cc      if (yntp.gt.stYMAX) then
      if (yntp.gt.stYMAX.or.yntp.lt.stYMIN) then ! s.fazio October 2010
      ioutkin = 1
      goto 999
      endif

* -> incoming lepton in lab frame
      pli_l(4) =  elep
      pli_l(3) = -elep
      pli_l(1) = 0.d0
      pli_l(2) = 0.d0
      call v_copy(pli_l,plilab)

* ->  radiated photon in lab frame
      pvisr(1) = 0.
      pvisr(2) = 0.
      pvisr(3) = -real(egamr)
      pvisr(4) =  real(egamr)

* -> incoming proton in lab frame
      ppi_l(4) = dsqrt(eproton**2 + mp**2)
      ppi_l(3) = eproton
      ppi_l(1) = 0.0d0
      ppi_l(2) = 0.0d0
      call v_copy(ppi_l,ppilab)

* -> outcoming lepton in lab frame
      plo_l(4) = elep*(1.-y) + q2 / (4.*elep)
      thetal   = dacos(1. - 2.*elep*(1.-y)/plo_l(4) )
      phil     = dble(rndm(0.)) *2.d0 * pi
      plo_l(1) = plo_l(4) * dsin(thetal) * dcos(phil)
      plo_l(2) = plo_l(4) * dsin(thetal) * dsin(phil)
      plo_l(3) = plo_l(4) * dcos(thetal)
      call v_copy(plo_l,plolab)

* --------------------------------------
ccc      if (plo_l(4).lt.stELMIN) then
ccc       ioutkin = 1
ccc       goto 999
ccc      endif
* --------------------------------------

* -> virtual photon in lab frame
      do i=1,4
       pvgam_l(i) = pli_l(i) - plo_l(i)
      enddo
      call v_copy(pvgam_l,pvglab)

      if (idebug.eq.1) then
       write(6,*) '     '
       write(6,*) 'START'
       write(6,*) '     '
       write(6,*) ' ------   In lab frame :'
       write(6,*) 'Incoming lepton :'
       call v_print(pli_l)
       print *,elepi,abs(stELEP)
       write(6,*) 'Incoming proton : '
       call v_print(ppi_l)
       print *,eproton,stETARG
       write(6,*) 'Outcoming lepton : '
       call v_print(plo_l)
       print *,y,q2,elep*(1.-y) + q2 / (4.*elep),
     +         1. - 2.*elep*(1.-y)/plo_l(4)
       write(6,*) 'Virtual photon : '
       call v_print(pvgam_l)
      endif

* --- Boost to p-rest frame, same system axis ---

      boostp = eproton / mp
      beta(1) =  -ppi_l(1) / ppi_l(4)
      beta(2) =  -ppi_l(2) / ppi_l(4)
      beta(3) =  -ppi_l(3) / ppi_l(4)
cls>>
cls      beta(1) =  0.0d0
cls      beta(2) =  0.0d0
cls      beta(3) = -eproton / dsqrt(eproton**2 + 0.938**2)
cls>>
      do i=1,3
       betar(i) = -beta(i)
      enddo

* -> Boost incoming lepton
      do i=1,4
       ptmp(i)    = pli_l(i)
      enddo
      call boost(ptmp,beta)
      do i=1,4
       pli_tmp(i) = ptmp(i)
      enddo
      call v_copy(pli_tmp,pliprf)

* -> Boost outcoming lepton
      do i=1,4
       ptmp(i)    = plo_l(i)
      enddo
      call boost(ptmp,beta)
      do i=1,4
       plo_tmp(i) = ptmp(i)
      enddo
      call v_copy(plo_tmp,ploprf)

* -> Boost incoming proton
      do i=1,4
       ptmp(i)    = ppi_l(i)
      enddo
      call boost(ptmp,beta)
      do i=1,4
       ppi_tmp(i) = ptmp(i)
      enddo
      call v_copy(ppi_tmp,ppiprf)

* -> Boost virtual photon
      do i=1,4
       ptmp(i)      = pvgam_l(i)
      enddo
      call boost(ptmp,beta)
      do i=1,4
       pvgam_tmp(i) = ptmp(i)
      enddo
      call v_copy(pvgam_tmp,pvgprf)

      if (idebug.eq.1) then
       write(6,*) ' '
       write(6,*) '----- In p-rest frame : '
       write(6,*) 'incoming lepton :'
       call v_print(pli_tmp)
       write(6,*) 'incoming proton :'
       call v_print(ppi_tmp)
       write(6,*) 'outcoming lepton :'
       call v_print(plo_tmp)
       write(6,*) 'virtual photon :'
       call v_print(pvgam_tmp)
      endif

* --- Rotate the system axis, z = -q_1^3 ---

      do i=1,3
       rot(i) = -pvgam_tmp(i)
      enddo

      do_xaxis = .true.
      i_xaxis  = 1
      call rotate(pli_tmp,rot,1,pli_b)
      do_xaxis = .false.
      call v_copy(pli_b,plibel)

      call rotate(plo_tmp,rot,1,plo_b)
      call rotate(pvgam_tmp,rot,1,pvgam_b)
      call rotate(ppi_tmp,rot,1,ppi_b)

      call v_copy(plo_b,plobel)
      call v_copy(pvgam_b,pvgbel)
      call v_copy(ppi_b,ppibel)

      if (idebug.eq.1) then
       write(6,*) ' '
       write(6,*) '-----  in Belitsky frame : '
       write(6,*) 'incoming lepton :'
       call v_print(pli_b)
       write(6,*) 'incoming proton :'
       call v_print(ppi_b)
       write(6,*) 'outcoming lepton :'
       call v_print(plo_b)
       write(6,*) 'virtual photon :'
       call v_print(pvgam_b)
      endif

* -> Outcoming proton, in Belitsky frame :

      if (stIELAS.eq.1) then

        ppo_b(4) = mp - t / (2.0d0 * mp)
        qpo_b    = dsqrt( t*(t-4*mp**2) ) / (2.d0*mp)

        p2   = qpo_b
        q13  = - (q2 / (2.*mp*xb))*dsqrt(1.+4.*mp**2*xb**2/q2)
        q20  = (t + q2/xb ) / (2.*mp)

      else

        ppo_b(4) = (dym**2+mp**2) / (2.*mp)
        ppo_b(4) = ppo_b(4) - t / (2.0d0 * mp)

        qpo_b = dsqrt(ppo_b(4)**2 -dym**2)
        p2   = qpo_b
        q13  = - (q2 / (2.*mp*xb))*dsqrt(1.+4.*mp**2*xb**2/q2)
        q20  = (t + q2/xb -(dym**2-mp**2) ) / (2.*mp)

      endif

      ccc  = q13**2-q20**2+p2**2
      ccc  = ccc / (2. * p2 * q13)
      cost = ccc

      if ( (cost.gt.1).or.(cost.lt.-1.) ) then
       write(6,*) 'pb cost ',cost,dym,xb,q2,t,phi
       ibad    = 1
       ioutkin = 1
       goto 999
      endif

      thetan = dacos(cost)

CLS>> NEW
      phi_p = dble(atan2(plobel(2),plobel(1))) -
     +        dble(rndm(0.)) *2.d0 * pi
CLS   =========================================

      ppo_b(1) = qpo_b * dcos(phi_p) * dsin(thetan)
      ppo_b(2) = qpo_b * dsin(phi_p) * dsin(thetan)
      ppo_b(3) = qpo_b * dcos(thetan)

      call v_copy(ppo_b,ppobel)
      if (idebug.eq.1) then
       write(6,*) 'Outcoming proton :'
       call v_print(ppo_b)
       call v_print4(ppobel)
      endif

* -> Outcoming (real) photon, in Belitsky frame :
      do i=1,4
       prgam_b(i) = pli_b(i) + ppi_b(i) - plo_b(i) - ppo_b(i)
      enddo

      call v_copy(prgam_b,prgbel)

* --- Rotate back the system axis
      call rotate(ppo_b,rot,-1,ppo_tmp)
      call v_copy(ppo_tmp,ppoprf)
      call rotate(prgam_b,rot,-1,prgam_tmp)
      call v_copy(prgam_tmp,prgprf)

      if (idebug.eq.1) then
       write(6,*) ' '
       write(6,*) '----- Back to standard axis system (checks) '
       call rotate(pvgam_b,rot,-1,toto_vgam)
       write(6,*) 'virt gamma back to standard axis (check) '
       call v_print(toto_vgam)
       call rotate(pli_b,rot,-1,toto_li)
       write(6,*) 'incoming lepton back to st axis (check) '
       call v_print(toto_li)
       call rotate(plo_b,rot,-1,toto_lo)
       write(6,*) 'outcoming lepton back to st axis (check) '
       call v_print(toto_lo)
       write(6,*) 'Outcoming proton :'
       call v_print4(ppoprf)
       write(6,*) 'Real photon :'
       call v_print4(prgprf)
      endif

* --- Boost back from p-rest frame to lab frame

      call boost(ppo_tmp,betar)
      do i=1,4
       ppo_l(i)   = ppo_tmp(i)
      enddo

      call boost(prgam_tmp,betar)
      do i=1,4
       prgam_l(i) = prgam_tmp(i)
      enddo

      call v_copy(ppo_l,ppolab)
      call v_copy(prgam_l,prglab)

      if (idebug.eq.1) then
       write(6,*) ' '
       write(6,*) '----- Back to lab frame (checks) '
       call boost(toto_vgam,betar)
       write(6,*) 'Virtual photon : '
       call v_print(toto_vgam)
       call boost(toto_li,betar)
       write(6,*) 'Incoming lepton :'
       call v_print(toto_li)
       call boost(toto_lo,betar)
       write(6,*) 'Outcoming lepton :'
       call v_print(toto_lo)
       write(6,*) 'Outcoming proton :'
       call v_print4(ppolab)
       write(6,*) 'Real photon : '
       call v_print4(prglab)
       write(6,*) '   '
      endif

      if (idebug.eq.1) then
       call check(plilab,ppilab,plolab,ppolab,prglab)
       call check(pliprf,ppiprf,ploprf,ppoprf,prgprf)
       call check(plibel,ppibel,plobel,ppobel,prgbel)
       write(6,*) '   '
       write(6,*) 'END'
      endif

      pi_real = atan(1.0) * 4.0

      th_lo = atan2(sqrt(plolab(1)**2+plolab(2)**2),plolab(3))
      th_lo = th_lo*180./pi_real
      en_lo = plolab(4)
      th_rg = atan2(sqrt(prglab(1)**2+prglab(2)**2),prglab(3))
      th_rg = th_rg*180./pi_real
      en_rg = prglab(4)

      if (idebug.eq.1) print *,th_lo,en_lo,th_rg,en_rg

      en_su = en_lo+en_rg

      w_eg =      (plolab(4)+prglab(4))**2
      w_eg = w_eg-(plolab(1)+prglab(1))**2
      w_eg = w_eg-(plolab(2)+prglab(2))**2
      w_eg = w_eg-(plolab(3)+prglab(3))**2
      w_eg = sqrt(w_eg)

      if (idebug.eq.1) then
      print *,'M(eg) invariant mass = ',w_eg
      endif

      if (th_lo.lt.stTHLMIN) ioutkin = 1
      if (th_lo.gt.stTHLMAX) ioutkin = 1
cls      if (en_lo.lt.stELMIN)  ioutkin = 1
      if (th_rg.lt.stTHGMIN) ioutkin = 1
      if (th_rg.gt.stTHGMAX) ioutkin = 1
cls      if (en_rg.lt.stEGMIN)  ioutkin = 1

CLS>> NEW
      if (en_su.lt.stELMIN)  ioutkin = 1 ! v3 => see README file...

      if (w_eg.lt.stM12MIN)  ioutkin = 1
      if (en_lo.lt.stELMB)   ioutkin = 1
      if (en_rg.lt.stEGMB)   ioutkin = 1

      if (stIRAD.eq.1) then
      if (EGAMR.lt.stEIMIN)  ioutkin = 1
      endif

      isca = iprodsca

ccc      print *,'th_lo=',th_lo,'ioutkin=',ioutkin

999   continue

      return
      end

