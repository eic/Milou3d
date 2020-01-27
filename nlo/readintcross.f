      real function readintcross(vect_w_lo,vect_q2_lo,
     &                           vect_w_nlo,vect_q2_nlo,
     &                           vect_q2_ep)
      implicit none


      real vect_w_lo(100),  vect_q2_lo(100),
     &     vect_w_nlo(100), vect_q2_nlo(100)

      real vect_q2_ep(100)


      integer nbins
      parameter(nbins=100)

      integer ic
      real rval
      character name*19, LINE*72

      name = 'txt/vect_w_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, nbins
         READ(95,100) LINE
         READ(LINE,101)  rval
         vect_w_lo(ic) = rval
      enddo
      close(95)

      name = 'txt/vect_q2_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, nbins
         READ(95,100) LINE
         READ(LINE,101)  rval
         vect_q2_lo(ic) = rval
      enddo
      close(95)

      name = 'txt/vect_w_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, nbins
         READ(95,100) LINE
         READ(LINE,101)  rval
         vect_w_nlo(ic) = rval
      enddo
      close(95)

      name = 'txt/vect_q2_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, nbins
         READ(95,100) LINE
         READ(LINE,101)  rval
         vect_q2_nlo(ic) = rval
      enddo
      close(95)

*
* ep cross section
*
      name = 'txt/vect_q2_ep.txt'
      
      print*,name
      
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, 100
         READ(95,100) LINE
         READ(LINE,101)  rval
         vect_q2_ep(ic) = rval

         print*,rval,vect_q2_ep(ic)

      enddo
      close(95)

 100  FORMAT(A72)
 101  FORMAT(F16.8)

      return
      end



