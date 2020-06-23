C     Read single file

      function readinsingle(nx,nq,nt,path,isim) result(arr)

      implicit none

      include 'dvcs.common'

      integer nx,nq,nt
      character(len=*) path
      logical isim
      real*8 arr(nx,nq,nt) 

      integer i,j,k
      real*8 x,q2,t
      real*8 dum1,dum2

      open(unit=11,file=stPATHGRID(1:stPATHGRIDLenght)//path,
     >  status='unknown')

      do i = 1, nx
            read(11,*) x

         do j = 1, nq
            read(11,*) q2
            
            do k = 1, nt
                read(11,*) t,dum1,dum2

                if(isim) then
                    arr(i,j,k) = dum2
                else
                    arr(i,j,k) = dum1
                end if

            enddo
         enddo
      enddo

      close(unit=11)

      end function 


