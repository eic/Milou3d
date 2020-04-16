      function interpolate3D(nx,nq,nt,xx,qq,tt,x,q2,t,arr) result(amp)

      implicit none

      integer mx,mq,mt
      parameter(mx=100,mq=100,mt=100)
    
      integer nx,nq,nt
      real*8 xx(mx)
      real*8 qq(mq)
      real*8 tt(mt)
      real*8 x,q2,t
      real*8 arr(mx,mq,mt)


      real point(3)
      integer arr_size(3)
      real grid(nx+nq+nt) 
      real arr_new(nx,nq,nt)

      integer i,j,k

      real amp

      real FINT

C POINT

      point(1) = LOG10(x)
      point(2) = LOG10(q2)
      point(3) = t

C SIZE
    
      arr_size(1) = nx
      arr_size(2) = nq
      arr_size(3) = nt

C GRID

      do i = 1, nx
        grid(i) = LOG10(xx(i)) 
      enddo

      do j = 1, nq
        grid(nx+j) = LOG10(qq(j)) 
      enddo

      do k = 1, nt
        grid(nx+nq+k) = tt(k) 
      enddo

      do i = 1, nx
      do j = 1, nq
      do k = 1, nt
        arr_new(i,j,k) = arr(i,j,k) 
      enddo
      enddo
      enddo

C INTERPOLATE

      amp = fint(3,point,arr_size,grid,arr_new)  

      end function
