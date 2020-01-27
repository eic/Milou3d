      subroutine routines

      write(*,*) 'loading routines'
      return
      end
***********************************************************************

      subroutine bhold(t,x,q2,res)
      implicit none

      double precision t, x, y, q2, res
      double precision jacob, ge, gm, tau, m_p
      double precision pi, alpha, s
      double precision fact

      pi    = 3.1415D0
      alpha = dble(1./137.)
      s     = 90200.0D0
      m_p   = dble(0.938)


      y = q2/(s*x)


      ge  = (1.+(dabs(t)/(0.71D0)))**(-2.0D0)
      gm  = 2.7D0*ge
      tau = dabs(t)/(4.0D0*(m_p**2))

      jacob = (1.D0/(x*s))

      res = (((alpha**3)*s*(y**2)*(1.D0+(1.D0-y)**2))/
     &      (pi*(q2**2)*dabs(t)*(1.D0-y)))
     &     *((ge**2+tau*(gm**2))/(1.D0+tau))*jacob*2.D0*pi


c      print*,ge,gm,tau


      return
      end

******************************************************************
c     function fbh(x)
      subroutine fbh

      integer ordint
      real    xint, q2int
      double precision dt, dx, dq2, dres

      real    mp

      q2int=4.0 
c     xint=x
      wmin=30.
      wmax=120.
      dq2 = dble(q2int)
      write(*,*) 'W range',wmin,wmax

      mp=0.938

      fact = 0.38937966*1000000.

      iter=500
      aiter=real(npt)

      sigt=0.

c W integration
      xmin=q2int/(wmax**2-q2int)
      xmax=q2int/(wmin**2-q2int)
      write(*,*) 'x range',xmin,xmax
      alxmi=alog(xmin)
      alxma=alog(xmax)

      do ix=1,iter
       aix=1.*ix
       xb1=exp((aix*(alxma-alxmi)/aix)+alxmi)
       xb2=exp(((aix-1.)*(alxma-alxmi)/aix)+alxmi)
       x=(xb1-xb2)/2. + xb2

c      x=0.0003 
       dx=dble(x)

c t integration
       tmin=mp**2*x**2/(1.-x**2)
       tmax=1.0

       altmi=alog(tmin)
       altma=alog(tmax)

       do it=1,iter
        ait=1.*it
        tb1=exp((ait*(altma-altmi)/ait)+altmi)
        tb2=exp(((ait-1.)*(altma-altmi)/ait)+altmi)
        t=(tb1-tb2)/2. + tb2

        dt=dble(t) 


        ordint=1
c       call diffcrossdvcs(ordint,dt,dx,dq2,dres)
        call bhold(dt,dx,dq2,dres)

c       print*, dx,dq2,dt,dres

        air=(tb1-tb2)
        air=air*(xb1-xb2)
c       sigt=sigt+real(dres)*air*t*x
        sigt=sigt+res*air*t*x
c       sigt=sigt+real(dres)*air*t

       enddo
c     print*, dx,dq2,dt
      enddo

c      print*, ordint,dx,dq2,dt,dres

      sig = sigt*fact
      write(*,*) sig,' nb'

      return
      end




