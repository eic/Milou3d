      real function crossread(vect1_lo,vect2_lo,vect3_lo,vect4_lo,
     &                     vect1_nlo,vect2_nlo,vect3_nlo,vect4_nlo)
      implicit none

      integer mbins
      parameter(mbins=100)

      integer ic
      real rval
      character name*17, LINE*72

      real vect1_lo(mbins),  vect2_lo(mbins),
     &     vect3_lo(mbins),  vect4_lo(mbins),
     &     vect1_nlo(mbins), vect2_nlo(mbins),
     &     vect3_nlo(mbins), vect4_nlo(mbins)


c      print*,'station 1'
      name = 'txt/vect1_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect1_lo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect1_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect1_nlo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect2_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect2_lo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect2_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect2_nlo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect3_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect3_lo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect3_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect3_nlo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect4_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect4_lo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect4_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect4_nlo(ic) = rval
      enddo
      close(95)



 100  FORMAT(A72)
 101  FORMAT(F16.8)



      return
      end
