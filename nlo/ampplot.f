      subroutine ampplot(xdummy)
      implicit none
      
      real xdummy
      character name*17
      real xg(10), q2g(10)

      integer nbins
      parameter (nbins=50)

      real vectreu(2,nbins,nbins),   vectimu(2,nbins,nbins), 
     &     vectred(2,nbins,nbins),   vectimd(2,nbins,nbins),
     &     vectres(2,nbins,nbins),   vectims(2,nbins,nbins),
     &     vectreg(2,nbins,nbins),   vectimg(2,nbins,nbins),
     &     vectreup(2,nbins,nbins),  vectimup(2,nbins,nbins), 
     &     vectredp(2,nbins,nbins),  vectimdp(2,nbins,nbins),
     &     vectresp(2,nbins,nbins),  vectimsp(2,nbins,nbins),
     &     vectregp(2,nbins,nbins),  vectimgp(2,nbins,nbins),
     &     vectreue(2,nbins,nbins),  vectimue(2,nbins,nbins),
     &     vectrede(2,nbins,nbins),  vectimde(2,nbins,nbins),
     &     vectrese(2,nbins,nbins),  vectimse(2,nbins,nbins),
     &     vectrege(2,nbins,nbins),  vectimge(2,nbins,nbins),
     &     vectreuep(2,nbins,nbins), vectimuep(2,nbins,nbins),
     &     vectredep(2,nbins,nbins), vectimdep(2,nbins,nbins),
     &     vectresep(2,nbins,nbins), vectimsep(2,nbins,nbins),
     &     vectregep(2,nbins,nbins), vectimgep(2,nbins,nbins)

      real x, q2
      integer i, j, k

      character LINE*72

      real rval

      data xg / 0.000288, 0.000355, 0.0005, 0.00075, 0.001,
     &          0.002,    0.00355,  0.005,  0.006,   0.0072 /

      data q2g /  1.9993960, 3.0765160, 4.3848360, 5.9194890,
     &            7.6895290, 9.6907690, 11.923209, 14.386849,
     &           17.073424, 19.998784/ 


 100  FORMAT(A72)
 101  FORMAT(F16.8)

      print*,'station 1'

      name = 'txt/vectreu.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectreu(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 2'

      name = 'txt/vectimu.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimu(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 3'

      name = 'txt/vectred.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectred(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 4'
 
      name = 'txt/vectimd.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimd(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 5'

      name = 'txt/vectres.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectres(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 6'

      name = 'txt/vectims.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectims(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 7'

      name = 'txt/vectreg.txt' 
c      print*,'yo 1 ',name
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
c      print*,'yo 2 ',name
      do i= 1, 2
c         print*,'yo 3 ',i
         do j = 1, nbins
c            print*,'yo 4 ',i,j
            do k = 1, nbins
c               print*,'yo 5 ',i,j,k
               READ(95,100) LINE
               READ(LINE,101) rval
               vectreg(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 8'

      name = 'txt/vectimg.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimg(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 9'

      name = 'txt/vectreup.txt' 
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectreup(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 10'

      name = 'txt/vectimup.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimup(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 11'

      name = 'txt/vectredp.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectredp(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 12'

      name = 'txt/vectimdp.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimdp(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 13'

      name = 'txt/vectresp.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectresp(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 14'

      name = 'txt/vectimsp.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimsp(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 15'

      name = 'txt/vectregp.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectregp(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 16'

      name = 'txt/vectimgp.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimgp(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 17'

      name = 'txt/vectreue.txt' 
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectreue(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 18'

      name = 'txt/vectimue.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimue(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 19'

      name = 'txt/vectrede.txt' 
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectrede(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 20'

      name = 'txt/vectimde.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimde(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 21'

      name = 'txt/vectrese.txt'  
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectrese(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 22'

      name = 'txt/vectimse.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimse(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 23'

      name = 'txt/vectrege.txt' 
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectrege(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 24'

      name = 'txt/vectimge.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimge(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 25'

      name = 'txt/vectreuep.txt' 
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectreuep(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 26'

      name = 'txt/vectimuep.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimuep(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 27'

      name = 'txt/vectredep.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectredep(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 28'

      name = 'txt/vectimdep.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimdep(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 29'

      name = 'txt/vectresep.txt' 
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectresep(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 30'

      name = 'txt/vectimsep.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimsep(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 31'

      name = 'txt/vectregep.txt' 
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectregep(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station 32'

      name = 'txt/vectimgep.txt'
      OPEN (UNIT=95,FILE=name,STATUS='OLD')
      do i= 1, 2
         do j = 1, nbins
            do k = 1, nbins
               READ(95,100) LINE
               READ(LINE,101) rval
               vectimgep(i,j,k) = rval
            enddo
         enddo
      enddo
      close(95)

      print*,'station end of reading'


*
* filling the histogramms
*
      do j = 1, nbins
         do k = 1, nbins
            x   = 10**(log10(xg(1))
     &           + real(j)*(log10(xg(10))-log10(xg(1)))/real(nbins))
            q2  = q2g(1)+real(k)*(q2g(10)-q2g(1))/real(nbins)
            

            call hfill(101,alog10(x),q2,vectreu(1,j,k))
            call hfill(102,alog10(x),q2,vectimu(1,j,k))
            call hfill(103,alog10(x),q2,vectred(1,j,k))
            call hfill(104,alog10(x),q2,vectimd(1,j,k))
            call hfill(105,alog10(x),q2,vectres(1,j,k))
            call hfill(106,alog10(x),q2,vectims(1,j,k))
            call hfill(107,alog10(x),q2,vectreg(1,j,k))
            call hfill(108,alog10(x),q2,vectimg(1,j,k))
            call hfill(109,alog10(x),q2,vectreup(1,j,k))
            call hfill(110,alog10(x),q2,vectimup(1,j,k))
            call hfill(111,alog10(x),q2,vectredp(1,j,k))
            call hfill(112,alog10(x),q2,vectimdp(1,j,k))
            call hfill(113,alog10(x),q2,vectresp(1,j,k))
            call hfill(114,alog10(x),q2,vectimsp(1,j,k))
            call hfill(115,alog10(x),q2,vectregp(1,j,k))
            call hfill(116,alog10(x),q2,vectimgp(1,j,k))
            call hfill(117,alog10(x),q2,vectreue(1,j,k))
            call hfill(118,alog10(x),q2,vectimue(1,j,k))
            call hfill(119,alog10(x),q2,vectrede(1,j,k))
            call hfill(120,alog10(x),q2,vectimde(1,j,k))
            call hfill(121,alog10(x),q2,vectrese(1,j,k))
            call hfill(122,alog10(x),q2,vectimse(1,j,k))
            call hfill(123,alog10(x),q2,vectrege(1,j,k))
            call hfill(124,alog10(x),q2,vectimge(1,j,k))
            call hfill(125,alog10(x),q2,vectreuep(1,j,k))
            call hfill(126,alog10(x),q2,vectimuep(1,j,k))
            call hfill(127,alog10(x),q2,vectredep(1,j,k))
            call hfill(128,alog10(x),q2,vectimdep(1,j,k))
            call hfill(129,alog10(x),q2,vectresep(1,j,k))
            call hfill(130,alog10(x),q2,vectimsep(1,j,k))
            call hfill(131,alog10(x),q2,vectregep(1,j,k))
            call hfill(132,alog10(x),q2,vectimgep(1,j,k))


            call hfill(201,alog10(x),q2,vectreu(2,j,k))
            call hfill(202,alog10(x),q2,vectimu(2,j,k))
            call hfill(203,alog10(x),q2,vectred(2,j,k))
            call hfill(204,alog10(x),q2,vectimd(2,j,k))
            call hfill(205,alog10(x),q2,vectres(2,j,k))
            call hfill(206,alog10(x),q2,vectims(2,j,k))
            call hfill(207,alog10(x),q2,vectreg(2,j,k))
            call hfill(208,alog10(x),q2,vectimg(2,j,k))
            call hfill(209,alog10(x),q2,vectreup(2,j,k))
            call hfill(210,alog10(x),q2,vectimup(2,j,k))
            call hfill(211,alog10(x),q2,vectredp(2,j,k))
            call hfill(212,alog10(x),q2,vectimdp(2,j,k))
            call hfill(213,alog10(x),q2,vectresp(2,j,k))
            call hfill(214,alog10(x),q2,vectimsp(2,j,k))
            call hfill(215,alog10(x),q2,vectregp(2,j,k))
            call hfill(216,alog10(x),q2,vectimgp(2,j,k))
            call hfill(217,alog10(x),q2,vectreue(2,j,k))
            call hfill(218,alog10(x),q2,vectimue(2,j,k))
            call hfill(219,alog10(x),q2,vectrede(2,j,k))
            call hfill(220,alog10(x),q2,vectimde(2,j,k))
            call hfill(221,alog10(x),q2,vectrese(2,j,k))
            call hfill(222,alog10(x),q2,vectimse(2,j,k))
            call hfill(223,alog10(x),q2,vectrege(2,j,k))
            call hfill(224,alog10(x),q2,vectimge(2,j,k))
            call hfill(225,alog10(x),q2,vectreuep(2,j,k))
            call hfill(226,alog10(x),q2,vectimuep(2,j,k))
            call hfill(227,alog10(x),q2,vectredep(2,j,k))
            call hfill(228,alog10(x),q2,vectimdep(2,j,k))
            call hfill(229,alog10(x),q2,vectresep(2,j,k))
            call hfill(230,alog10(x),q2,vectimsep(2,j,k))
            call hfill(231,alog10(x),q2,vectregep(2,j,k))
            call hfill(232,alog10(x),q2,vectimgep(2,j,k))

         enddo
      enddo


      return
      end
