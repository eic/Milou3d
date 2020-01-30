C     Read single file

      function readinsingle(nx,nq,nt,path,isim) result(arr)

      implicit none

      integer nx,nq,nt
      character(len=*) path
      logical isim
      real*8 arr(nx,nq,nt) 

      integer i,j,k
      real*8 x,q2,t
      real*8 dum1,dum2

      open(unit=11,file=path,status='unknown')

      do j = 1, nx
            read(11,*) x

         do i = 1, nq
            read(11,*) q2
            
            do k = 1, nt
                read(11,*) t,dum1,dum2

                if(isim) then
                    arr(j,i,k) = dum1
                else
                    arr(j,i,k) = dum2
                end if

            enddo
         enddo
      enddo

      close(unit=11)

      end function 


