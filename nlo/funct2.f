      double precision function funct2(var1,var2)

      double precision var1, var2

      integer ordint
      real    xint, q2int
      common / varint / ordint, xint, q2int

      double precision dt, dx, dq2, dres


      dt  = var1
      dx  = var2
      dq2 = dble(q2int)

c      call diffcrossdvcs(ordint,dt,dx,dq2,dres)

c      call bhold(var1,var2,dq2,dres)

c      print*,ordint,var1,var2,dq2,dres

      dres = 3.0D0



      funct2 = dres

      return
      end
