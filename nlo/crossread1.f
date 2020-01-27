      real function crossread1(vect5_lo,vect6_lo,vect7_lo,vect8_lo,
     &                     vect5_nlo,vect6_nlo,vect7_nlo,vect8_nlo)
      implicit none

      integer mbins
      parameter(mbins=100)

      integer ic
      real rval
      character name*17, LINE*72

      real vect5_lo(mbins),  vect6_lo(mbins),
     &     vect7_lo(mbins),  vect8_lo(mbins),
     &     vect5_nlo(mbins), vect6_nlo(mbins),
     &     vect7_nlo(mbins), vect8_nlo(mbins)


c      print*,'station 1'
      name = 'txt/vect5_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect5_lo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect5_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect5_nlo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect6_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect6_lo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect6_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect6_nlo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect7_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect7_lo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect7_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect7_nlo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect8_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect8_lo(ic) = rval
      enddo
      close(95)

c      print*,'station 1'
      name = 'txt/vect8_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do ic = 1, mbins
         READ(95,100) LINE
         READ(LINE,101) rval
         vect8_nlo(ic) = rval
      enddo
      close(95)



 100  FORMAT(A72)
 101  FORMAT(F16.8)



      return
      end
