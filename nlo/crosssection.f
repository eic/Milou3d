      subroutine crosssection
      implicit none

      integer mbins
      parameter (mbins=100)

      real vect1_lo(mbins), vect2_lo(mbins), 
     &     vect3_lo(mbins), vect4_lo(mbins)

      real vect1_nlo(mbins), vect2_nlo(mbins), 
     &     vect3_nlo(mbins), vect4_nlo(mbins)

      real vect5_lo(mbins), vect6_lo(mbins), 
     &     vect7_lo(mbins), vect8_lo(mbins)

      real vect5_nlo(mbins), vect6_nlo(mbins), 
     &     vect7_nlo(mbins), vect8_nlo(mbins)


      real vect1_bh(mbins), vect2_bh(mbins),
     &     vect3_bh(mbins), vect4_bh(mbins)

      double precision x, q2, t, res
      
      real q2_1, q2_2
      real x_1,  x_2
      real xmin, xmax
      real t_const
      integer ord
      
      real fact

      character name*17

      integer ic

      fact = 0.38937966*1000000.

      q2_1 = 4.0
      q2_2 = 9.0
      x_1  = 0.0003
      x_2  = 0.007

      xmin = 0.000288
      xmax = 0.00719

      t_const = -0.01


      q2 = dble(q2_1)
      x  = dble(x_1)
      t  = -0.01
      ord = 1
      call diffcrossdvcs(ord,t,x,q2,res)
c      print*,res

c      return

      do ic = 1, mbins
         t = dble(-1.*real(ic)/real(mbins))

         q2 = dble(q2_1)
         x  = dble(x_1)
         ord = 1
         call diffcrossdvcs(ord,t,x,q2,res)
         vect1_lo(ic) = real(res)*fact

         q2 = dble(q2_2)
         x  = dble(x_1)
         ord = 1
         call diffcrossdvcs(ord,t,x,q2,res)
         vect2_lo(ic) = real(res)*fact

         q2 = dble(q2_1)
         x  = dble(x_2)
         ord = 1
         call diffcrossdvcs(ord,t,x,q2,res)
         vect3_lo(ic) = real(res)*fact

         q2 = dble(q2_2)
         x  = dble(x_2)
         ord = 1
         call diffcrossdvcs(ord,t,x,q2,res)
         vect4_lo(ic) = real(res)*fact


         q2 = dble(q2_1)
         x  = dble(x_1)
         ord = 2
         call diffcrossdvcs(ord,t,x,q2,res)
         vect1_nlo(ic) = real(res)*fact

         q2 = dble(q2_2)
         x  = dble(x_1)
         ord = 2
         call diffcrossdvcs(ord,t,x,q2,res)
         vect2_nlo(ic) = real(res)*fact

         q2 = dble(q2_1)
         x  = dble(x_2)
         ord = 2
         call diffcrossdvcs(ord,t,x,q2,res)
         vect3_nlo(ic) = real(res)*fact

         q2 = dble(q2_2)
         x  = dble(x_2)
         ord = 2
         call diffcrossdvcs(ord,t,x,q2,res)
         vect4_nlo(ic) = real(res)*fact

*
* BH cross check
*

         q2 = dble(q2_1)
         x  = dble(x_1)
         call bhold(t,x,q2,res)
         vect1_bh(ic) = real(res)*fact

      enddo

      do ic = 1, mbins

         q2 = dble(2.0+real(ic)*17.9/real(mbins))
         x  = dble(x_1)
         t  = dble(t_const)
         ord = 1
         call diffcrossdvcs(ord,t,x,q2,res)
         vect5_lo(ic) = real(res)*fact

         q2 = dble(2.0+real(ic)*17.9/real(mbins))
         x  = dble(x_1)
         t  = dble(t_const)
         ord = 2
         call diffcrossdvcs(ord,t,x,q2,res)
         vect5_nlo(ic) = real(res)*fact

         q2 = dble(2.0+real(ic)*17.9/real(mbins))
         x  = dble(x_2)
         t  = dble(t_const)
         ord = 1
         call diffcrossdvcs(ord,t,x,q2,res)
         vect6_lo(ic) = real(res)*fact

         q2 = dble(2.0+real(ic)*17.9/real(mbins))
         x  = dble(x_2)
         t  = dble(t_const)
         ord = 2
         call diffcrossdvcs(ord,t,x,q2,res)
         vect6_nlo(ic) = real(res)*fact

         q2 = dble(q2_1)
         x  = dble(10**(log10(xmin)
     &            +real(ic)*(log10(xmax)-log10(xmin))/real(mbins)))
         t  = dble(t_const)
         ord = 1
         call diffcrossdvcs(ord,t,x,q2,res)
         vect7_lo(ic) = real(res)*fact

         q2 = dble(q2_1)
         x  = dble(10**(log10(xmin)
     &            +real(ic)*(log10(xmax)-log10(xmin))/real(mbins)))
         t  = dble(t_const)
         ord = 2
         call diffcrossdvcs(ord,t,x,q2,res)
         vect7_nlo(ic) = real(res)*fact

         q2 = dble(q2_2)
         x  = dble(10**(log10(xmin)
     &            +real(ic)*(log10(xmax)-log10(xmin))/real(mbins)))
         t  = dble(t_const)
         ord = 1
         call diffcrossdvcs(ord,t,x,q2,res)
         vect8_lo(ic) = real(res)*fact

         q2 = dble(q2_2)
         x  = dble(10**(log10(xmin)
     &            +real(ic)*(log10(xmax)-log10(xmin))/real(mbins)))
         t  = dble(t_const)
         ord = 2
         call diffcrossdvcs(ord,t,x,q2,res)
         vect8_nlo(ic) = real(res)*fact

      enddo

*
* write cross section to files
*

 101  FORMAT(F16.8)
 102  FORMAT(F18.8)

      name  = 'txt/vect1_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect1_lo(ic)
      enddo
      CLOSE(95)


      name  = 'txt/vect2_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect2_lo(ic)
      enddo
      CLOSE(95)


      name  = 'txt/vect3_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect3_lo(ic)
      enddo
      CLOSE(95)


      name  = 'txt/vect4_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect4_lo(ic)
      enddo
      CLOSE(95)


      name  = 'txt/vect1_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect1_nlo(ic)
      enddo
      CLOSE(95)


      name  = 'txt/vect2_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect2_nlo(ic)
      enddo
      CLOSE(95)


      name  = 'txt/vect3_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect3_nlo(ic)
      enddo
      CLOSE(95)


      name  = 'txt/vect4_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect4_nlo(ic)
      enddo
      CLOSE(95)

      name  = 'txt/vect5_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect5_lo(ic)
      enddo
      CLOSE(95)

      name  = 'txt/vect5_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect5_nlo(ic)
      enddo
      CLOSE(95)

      name  = 'txt/vect6_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect6_lo(ic)
      enddo
      CLOSE(95)

      name  = 'txt/vect6_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect6_nlo(ic)
      enddo
      CLOSE(95)

      name  = 'txt/vect7_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect7_lo(ic)
      enddo
      CLOSE(95)

      name  = 'txt/vect7_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect7_nlo(ic)
      enddo
      CLOSE(95)

      name  = 'txt/vect8_lo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect8_lo(ic)
      enddo
      CLOSE(95)

      name  = 'txt/vect8_nlo.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,101) vect8_nlo(ic)
      enddo
      CLOSE(95)

      name  = 'txt/vect1_bh.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do ic = 1, mbins
         write(95,102) vect1_bh(ic)
      enddo
      CLOSE(95)

       
      return
      end
