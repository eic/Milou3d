      function funct(x)

      real    xint, q2int
      double precision dt, dx, dq2, dres

      q2int=4.5 
      dx  = dble(x)
      dq2 = dble(q2int)

      tmin=0.01
      tmax=1.0

c t integration
      altmi=alog(tmin)
      altma=alog(tmax)
      iter=100
      aiter=real(npt)

      sigt=0.

      do it=1,iter
       ait=1.*it
       tb1=exp((ait*(altma-altmi)/ait)+altmi)
       tb2=exp(((ait-1.)*(altma-altmi)/ait)+altmi)
       t=(tb1-tb2)/2. + tb2

       dt=dble(t) 

       print*, ordint,dx,dq2,dt

       ordint=1
c      call diffcrossdvcs(ordint,dt,dx,dq2,dres)
       call bhold(ordint,dt,dx,dq2,dres)

       air=(tb1-tb2)/ait
       sigt=sigt+real(dres)*air*t

      enddo

c      print*, ordint,dx,dq2,dt,dres

      funct = sigt

      return
      end




