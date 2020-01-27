      subroutine readvect(name,vect_1,vect_2)
      implicit none

      character name*19
      real vect_1(10,2), vect_2(10,10,2)

      character*72 LINE
      real rval, rval1, rval2, rval3
      integer ic, jc


c      print*,name

      
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      
c      print*,name

      do ic = 1, 10
c         print*,name
         READ(95,100) LINE
c         print*,line
         READ(LINE,101) rval
         vect_1(ic,1) = rval
         do jc = 1, 10
            READ(95,100) LINE
c            print*,line
            READ(LINE,102) rval1, rval2, rval3  
c            print*,'num ',rval1,rval2,rval3
            vect_1(jc,2) = rval1
            vect_2(ic,jc,1) = rval2
            vect_2(ic,jc,2) = rval3
c            print*,'vec ',vect_1(jc,2),vect_2(ic,jc,1),vect_2(ic,jc,2)
         enddo
c         print*,'end of loop'
      enddo

 100  FORMAT(A72)
 101  FORMAT(F14.8)
 102  FORMAT(3F16.8)

      CLOSE(95)

c      print*,'end of routine ',name

      return
      end
