****************************************************************
*                                                              *      
*                         REU                                  * 
*                                                              *
****************************************************************

      double precision function reu(ord,xi,q2i)

      implicit none

c      include 'func.include'      

      integer ord
      double precision xi, q2i
c      double precision amp
c      double precision x(10), q2(10)
c      double precision interpol

c      logical first

c      double precision vreu(2,10,10),   vimu(2,10,10), 
c     &                 vred(2,10,10),   vimd(2,10,10),
c     &                 vres(2,10,10),   vims(2,10,10),
c     &                 vreg(2,10,10),   vimg(2,10,10),
c     &                 vreup(2,10,10),  vimup(2,10,10), 
c     &                 vredp(2,10,10),  vimdp(2,10,10),
c     &                 vresp(2,10,10),  vimsp(2,10,10),
c     &                 vregp(2,10,10),  vimgp(2,10,10),
c     &                 vreue(2,10,10),  vimue(2,10,10),
c     &                 vrede(2,10,10),  vimde(2,10,10),
c     &                 vrese(2,10,10),  vimse(2,10,10),
c     &                 vrege(2,10,10),  vimge(2,10,10),
c     &                 vreuep(2,10,10), vimuep(2,10,10),
c     &                 vredep(2,10,10), vimdep(2,10,10),
c     &                 vresep(2,10,10), vimsep(2,10,10),
c     &                 vregep(2,10,10), vimgep(2,10,10)
c
c      common /vectors/ vreu,   vimu, 
c     &                 vred,   vimd,
c     &                 vres,   vims,
c     &                 vreg,   vimg,
c     &                 vreup,  vimup, 
c     &                 vredp,  vimdp,
c     &                 vresp,  vimsp,
c     &                 vregp,  vimgp,
c     &                 vreue,  vimue,
c     &                 vrede,  vimde,
c     &                 vrese,  vimse,
c     &                 vrege,  vimge,
c     &                 vreuep, vimuep,
c     &                 vredep, vimdep,
c     &                 vresep, vimsep,
c     &                 vregep, vimgep



c      data x / 0.000288, 0.000355, 0.0005, 0.00075, 0.001,
c     &         0.002,    0.00355,  0.005,  0.006,   0.0072 /
c
c      data q2 /  1.9993960, 3.0765160, 4.3848360, 5.9194890,
c     &           7.6895290, 9.6907690, 11.923209, 14.386849,
c     &          17.073424, 19.998784/ 

c      data first /.true./
c      save first
      
      print*,'this is reu'

c      if (first) then
c         first = .false.
c         call readnlo
c         print*,'read vector'
c      endif
      
c      amp = interpol(ord,vreu,x,q2,xi,q2i)

c      reu = dble(amp)

      reu = 0.0

      return
      end

c      include 'interpol.f'
c      include 'readnlo.f'
