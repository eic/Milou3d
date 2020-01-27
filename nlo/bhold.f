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
