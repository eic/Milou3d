*CMZ :          13/07/2004  17.06.26  by  H1 Saclay
*-- Author :    Unknown   17/12/2003


* ---------------------------------
      real function fdvcs*8(x)
* ---------------------------------

      IMPLICIT REAL*8(A-H,O-Z)

      include 'dvcs.common'
      include 'forpaw.common'

      double precision mp

      PARAMETER (MXDIM = 50 )
      DIMENSION X(MXDIM)

      dimension res(100)

      common /setting/ nset
      common /tintin/ SIGN,bti,rti

      common /counter1/ icountxq,icountt

C ------------------------------------------------------
      common/dvcs_BOUNDS/xmin,xmax,qmin,qmax,tmin,tmax,
     >       phimin,phimax,pphimin,pphimax,thmin,thmax

      logical Angle_Integ
      common/dvcs_IN/spin,Z,A,s,tchar,xlaml,xlamp,
     >               Angle_Integ,iord,icount,nx,nq,IPRO

      common/dvcs_OUT/res

      common/dvcs_VAR/x_main,q_main,phi_main,t_main,ym_main
C ------------------------------------------------------

      real elepi,EGAMR,sweight,srad,elepin,ehadi
      common /RADGEN/ elepi,EGAMR,sweight,srad,elepin,ehadi

      logical stFIXED_save
      common /FIXED/ stFIXED_save

      integer i_outrange
      common/out_of_range/i_outrange

      integer igen_exp
      common /genexp/ igen_exp

      logical first
      data first /.true./


      real shat,ytest,mel
      real wmin_test
      real smphibel
      REAL sig_smphibel

      fdvcs = 0.d0

      icountxq = 1
      icountt  = 1

*  X(1) -> x
*  X(2) -> q
*  X(3) -> t
*  X(4) -> Egamma (ISR)

*  X(5) -> phi
*  X(6) -> pphi
*  X(7) -> theta

cls>>
cls>> generate in linx,linq,lint
cls>>

      if (igen_exp.eq.1) then
      xx  = dexp(X(1))
      else
      xx  = X(1)
      endif
      if (igen_exp.eq.1) then
      q   = dexp(X(2))
      else
      q   = X(2)
      endif
      t   = X(3)

      x_main = xx
      q_main = q
      t_main = t

      if (Angle_Integ) then
       phi   =0.
       pphi  =0.
       theta =0.
       nset = 2
      else
       phi   = X(5)
       pphi  = X(6)
       theta = X(7)
      endif

      if (srad.eq.0.) then
       egamr = 0.
      else
       egamr = sngl(X(4))
      endif

      phi_main = phi

      xpi = datan(1.d0) * 4.0d0

C
C --- IN CASE OF RADGEN
C
C ---------------------------------------------
      elepi   = elepin
      sweight = 1.
      snew    = s

      if (srad.eq.1.) then

 13   continue
      elepi = elepin
      call QED_ISR(sngl(s),elepi,EGAMR,sweight)
      if (elepi.lt.5.) goto 13

cc      write(6,*) 'x4 eg we ',X(4),egamr,sweight

       mp     = 0.93827d0

      elep    = dble(elepi)
      eproton = dble(ehadi)
      if (stFIXED_save) then
      snew    = 2.0d0 * elep * mp + mp**2
      else
      snew    = 4.0d0 * elep * eproton
      endif

      endif

      weight = dble(sweight)
C ---------------------------------------------

      i_outrange = 0

CLS>> NEW


        ym_main = 0.93827d0

C ----- For dissociated proton :
      IF (stIELAS.eq.0) THEN
       call PXMASS (ym_main,ENHA)
      ENDIF

       mp     = 0.93827d0

      if (stIELAS.eq.1) then
       t_min  = -mp**2*x_main**2/(1.d0-x_main)
      else
cep       t_min = x_main*(ym_main**2-mp**2)/(1.d0-x_main)
       t_min = -ym_main**2 * x_main/(1.d0-x_main) + mp**2 * x_main
      endif


      if (t_main.gt.t_min) i_outrange = 1

      if (t_main.gt.t_min) then
ccc       write(6,*) 'out of tmin',t_min,t_main
      endif

      if (t_main.gt.t_min) goto 999


      xntp      = sngl(x_main)
      qntp      = sngl(q_main)
      phintp    = sngl(phi_main)
      tntp      = sngl(t_main)

      shat    = sngl(snew)
      ytest   = qntp/(xntp*(shat - sngl(mp)**2))


      wmin_test = 0.2*sqrt(qntp/xntp)
ccc      print *,ym_main,wmin_test
cls      if (sngl(ym_main).gt.wmin_test) then
cls      i_outrange = 1
cls      goto 999
cls      endif


      if ((ytest.gt.1.).or.(ytest.lt.0.)) then
        i_outrange = 1
        goto 999
      endif

CLS>>
      if (stFIXED_save) then
      call  calc_kinem2_fixed(ibad,ioutkin,isca)
      else
      call  calc_kinem2(ibad,ioutkin,isca)
      endif
cc      print *,ioutkin  
      if (ioutkin.eq.1)    goto 999
CLS>>

      mel = 0.000511
      if (qntp.lt.((mel**2*ytest**2)/(1.-ytest))) then
        i_outrange = 1
        goto 999
      endif

ccc      i_outrange = 0


      phi_newcalc_0 = dble(calcphirb(1.,0))
      phi_newcalc = dble(calcphirb(1.,1))
      
      phibel_0=sngl(phi_newcalc_0)
      phibel=sngl(phi_newcalc)

      resoltest = 100*(phibel_0-phibel)/phibel_0
      
      resol = real(resoltest)
      phibelrec=real(phibel)
      phibelgen=real(phibel_0)

ccc      phi_newcalc_0 = 0.
   
ccc      print *,phi_newcalc_0
 
c-------------------------------------
      call dvcsob(spin,Z,A,nx,nq,icount,xx,q,snew,t,phi,pphi,
     >     theta,tchar,xlaml,xlamp,iord,res,phi_newcalc_0)
c-------------------------------------
      icount = 2

      if (i_outrange.eq.1) goto 999


      if (ipro.eq.1) resu = res(4)   !  BH
      if (ipro.eq.2) resu = res(2)   !  DVCS
      if (ipro.eq.3) resu = res(3)   !  INT
      xint = 1.*res(3)
      if ((xint+res(4)+res(2)).lt.0.) xint=0.
      if (ipro.eq.4) resu = res(4)+res(2)+xint
      if (ipro.eq.4)
     +  resu = (res(4)+res(2))*
     +   (1+0.2*(-1.D0*SIGN)*DCOS(phi_newcalc_0))
ccc      print *,SIGN
ccc      if (ipro.eq.4) resu = res(7)   !  CA without TW3
      if (ipro.eq.5) resu = res(11)  ! SSA without TW3

ccc      print *,resu

      if (resu.lt.0.d0) then
        write(6,*) 'warning FDVCS<0'
        i_outrange = 1
      endif

      if (i_outrange.eq.1) goto 999

      if (igen_exp.eq.1) then
      resu = resu * xx
      resu = resu * q
      endif

      fdvcs = resu
      fdvcs = fdvcs*weight

ccc      write(6,*) X(1),X(2),X(3),fdvcs

c --- Fill histograms :

      CALL XHFILL(1,XX,fdvcs)
      CALL XHFILL(2,q,fdvcs)
      CALL XHFILL(3,t,fdvcs)
      CALL XHFILL(4,dble(egamr),fdvcs)

999   continue

cc      write(6,*) x_main,q_main,t_main,fdvcs

      return
      end

