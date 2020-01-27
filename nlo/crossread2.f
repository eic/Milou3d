      real function crossread2(vect1_bh,vect2_bh,vect3_bh,vect4_bh)

      implicit none

      integer mbins
      parameter(mbins=100)

      integer ic
      real rval
      character name*17, LINE*72

      real vect1_bh(mbins),  vect2_bh(mbins),
     &     vect3_bh(mbins),  vect4_bh(mbins)



c      print*,'station 1'
      name = 'txt/vect1_bh.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect1_bh(ic) = rval
      enddo
      close(95)

 100  FORMAT(A72)
 101  FORMAT(F18.8)



      return
      end
