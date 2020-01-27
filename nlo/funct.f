      double precision function funct(x)
      implicit none
      double precision x

      integer ordint
      real    xint, q2int
      common / varint / ordint, xint, q2int

      double precision dt, dx, dq2, dres

      dx  = dble(xint) 
      dq2 = dble(q2int)
      dt = x

      print*, ordint,dx,dq2,dt

      call diffcrossdvcs(ordint,dt,dx,dq2,dres)

c      print*, ordint,dx,dq2,dt,dres

c      funct = 2.
      


      funct = dres

c      print*,funct
      

      return
      end




