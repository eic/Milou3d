*CMZ :          02/06/2004  07.53.41  by  H1 Saclay
*-- Author :    Unknown   17/12/2003



      subroutine gener_evt


      implicit double precision (a-h,o-z)

      include 'dvcs.common'
      include 'forpaw.common'

      external fdvcs

      common /counter1/ icountxq,icountt

      dimension res(100)

      common/dvcs_OUT/res
      common/dvcs_VAR/x_main,q_main,phi_main,t_main,ym_main

      real*8 MPY

      integer ibad,ioutkin,isca

      logical FIRST
      data FIRST/.TRUE./

* --- COMMON FOR THE NTUPLE :  -> DVCSPaw.
* ----------------------------------------------------------

      MXTRY  = 50

*---Define the ntuple :
      NTID = 1

      if (FIRST) then
       if (stIELAS.eq.0) call DEFINE_RESO
       FIRST=.FALSE.
      endif

100   continue

      icountxq = 1
      icountt  = 1


      CALL SPRING( FDVCS, MXTRY )

      if (i_outrange.eq.1) goto 100

      xntp      = sngl(x_main)
      qntp      = sngl(q_main)
      phintp    = sngl(phi_main)
      tntp      = sngl(t_main)

      bhsq       = res(4)
      dvcssq     = res(2)
      xintdvcsbh = res(3)
      tot        = bhsq+dvcssq+xintdvcsbh

      if (stFIXED) then
      call calc_kinem2_fixed(ibad,ioutkin,isca)
      else
      call calc_kinem2(ibad,ioutkin,isca)
      endif

      wbhntp   = sngl(bhsq)
      wdvcsntp = sngl(dvcssq)
      wintntp  = sngl(xintdvcsbh)
      wtotntp  = sngl(tot)

      ym       = sngl(ym_main)

cls      print *,'gen : prglab(4) , xntp  ',prglab(4), xntp

      if (ibad.ne.0)     goto 100
      if (ioutkin.ne.0)  goto 100
      if (isca.ne.0)     goto 100


      call filllujet(ifail)

      call hfnt(NTID)


      return
      end
