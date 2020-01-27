      program dvcsnlo
      implicit none
      
      double precision t, x, q, res1, res2 , reu
      integer ord

      real fact, fres1, fres2
      
      external reu



*amplitudes

c      call amplitudes

* threefold differential cross section
      
c      call crosssection

* t integrated cross section

      call tintcross

      print*,'here i am'

      end

