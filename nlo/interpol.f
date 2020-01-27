      real function interpol(ord,vect1,x,q2,xi,q2i)
      
      integer ord
      double precision xi, q2i
      integer ic,jc, kc, xc, q2c
      integer nc
      double precision value
      double precision x(10), q2(10)
      double precision vect1(2,10,10)

* check values
c      print*,'station 1',xi,q2i

      if ((xi.lt.x(1)).or.(xi.gt.x(10))) then
         print*,'x not fitting ',xi,x(1),x(10)
      endif
      if ((q2i.lt.q2(1)).or.(q2i.gt.q2(10))) then
         print*,'q2 not fitting ',q2i,q2(1),q2(10)
      endif
      
c      print*,'station 2'

* determine bin

      do ic = 1,9
c         print*,'station x',ic,xi,x(ic),x(ic+1)
         if (xi.gt.x(ic).and.xi.le.x(ic+1)) then
            xc = ic
            goto 10
         endif
      enddo



 10   continue
c      print*,'station 3',q2i,q2      

      do jc = 1,9
         if (q2i.gt.q2(jc).and.q2i.le.q2(jc+1)) then
            q2c = jc
            goto 20
         endif
      enddo



 20   continue

c      print*,'station 4a'
c      print*,'station 4',xc,q2c,xi,q2i


* do interpolation

      if (q2i.lt.q2(q2c)
     &  + (((q2(q2c+1)-q2(q2c))/(x(xc+1)-x(xc)))*(xi-x(xc)))) then
         kc = 1
      else
         kc = 2
      endif

c      print*,'station 5',kc
c      print*,q2(q2c),q2(q2c+1),x(xc),x(xc+1)
c      print*,vect1(ord,xc,q2c),vect1(ord,xc,q2c+1)
c      print*,vect1(ord,xc+1,q2c),vect1(ord,xc+1,q2c+1)

      if (kc.eq.1) then


         value = real(vect1(ord,xc+1,q2c))
     &          + (q2i-q2(q2c)) *
     &            ((vect1(ord,xc+1,q2c+1) 
     &             -vect1(ord,xc+1,q2c))/
     &               (q2(q2c+1)-q2(q2c)))
     &             + (xi-x(xc+1))  *
     &            ((vect1(ord,xc+1,q2c) 
     &             -vect1(ord,xc,q2c))/
     &              (x(xc+1) - x(xc)))  
c         value=0.
      elseif (kc.eq.2) then
         value = vect1(ord,xc,q2c+1)
     &          + (q2i-q2(q2c+1)) *
     &            ((vect1(ord,xc,q2c+1) 
     &             -vect1(ord,xc,q2c))/
     &           (q2(q2c+1)-q2(q2c)))
     &          + (xi-x(xc)) *
     &           ((vect1(ord,xc+1,q2c+1) 
     &            -vect1(ord,xc,q2c+1))/
     &           (x(xc+1)-x(xc)))
c         value=0.0
      endif


      interpol = real(value)



      return
      end
