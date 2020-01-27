      subroutine readnlo
      implicit none

      integer ic, jc

      real ldamp_1_dat(10,2),      ldampe_1_dat(10,2),
     &     ldamppol_1_dat(10,2),   ldamppole_1_dat(10,2),
     &     lgamp_1_dat(10,2),      lgampe_1_dat(10,2),
     &     lgamppol_1_dat(10,2),   lgamppole_1_dat(10,2),
     &     lsamp_1_dat(10,2),      lsampe_1_dat(10,2),
     &     lsamppol_1_dat(10,2),   lsamppole_1_dat(10,2),
     &     luamp_1_dat(10,2),      luampe_1_dat(10,2),
     &     luamppol_1_dat(10,2),   luamppole_1_dat(10,2)
     
      real nlodamp_1_dat(10,2),    nlodampe_1_dat(10,2),
     &     nlodamppol_1_dat(10,2), nlodamppole_1_dat(10,2),
     &     nlogamp_1_dat(10,2),    nlogampe_1_dat(10,2),
     &     nlogamppol_1_dat(10,2), nlogamppole_1_dat(10,2),
     &     nlosamp_1_dat(10,2),    nlosampe_1_dat(10,2),
     &     nlosamppol_1_dat(10,2), nlosamppole_1_dat(10,2),
     &     nlouamp_1_dat(10,2),    nlouampe_1_dat(10,2),
     &     nlouamppol_1_dat(10,2), nlouamppole_1_dat(10,2)

      real ldamp_2_dat(10,10,2),      ldampe_2_dat(10,10,2),
     &     ldamppol_2_dat(10,10,2),   ldamppole_2_dat(10,10,2),
     &     lgamp_2_dat(10,10,2),      lgampe_2_dat(10,10,2),
     &     lgamppol_2_dat(10,10,2),   lgamppole_2_dat(10,10,2),
     &     lsamp_2_dat(10,10,2),      lsampe_2_dat(10,10,2),
     &     lsamppol_2_dat(10,10,2),   lsamppole_2_dat(10,10,2),
     &     luamp_2_dat(10,10,2),      luampe_2_dat(10,10,2),
     &     luamppol_2_dat(10,10,2),   luamppole_2_dat(10,10,2)
     
      real nlodamp_2_dat(10,10,2),    nlodampe_2_dat(10,10,2),
     &     nlodamppol_2_dat(10,10,2), nlodamppole_2_dat(10,10,2),
     &     nlogamp_2_dat(10,10,2),    nlogampe_2_dat(10,10,2),
     &     nlogamppol_2_dat(10,10,2), nlogamppole_2_dat(10,10,2),
     &     nlosamp_2_dat(10,10,2),    nlosampe_2_dat(10,10,2),
     &     nlosamppol_2_dat(10,10,2), nlosamppole_2_dat(10,10,2),
     &     nlouamp_2_dat(10,10,2),    nlouampe_2_dat(10,10,2),
     &     nlouamppol_2_dat(10,10,2), nlouamppole_2_dat(10,10,2)

      double precision vreu(2,10,10),   vimu(2,10,10), 
     &                 vred(2,10,10),   vimd(2,10,10),
     &                 vres(2,10,10),   vims(2,10,10),
     &                 vreg(2,10,10),   vimg(2,10,10),
     &                 vreup(2,10,10),  vimup(2,10,10), 
     &                 vredp(2,10,10),  vimdp(2,10,10),
     &                 vresp(2,10,10),  vimsp(2,10,10),
     &                 vregp(2,10,10),  vimgp(2,10,10),
     &                 vreue(2,10,10),  vimue(2,10,10),
     &                 vrede(2,10,10),  vimde(2,10,10),
     &                 vrese(2,10,10),  vimse(2,10,10),
     &                 vrege(2,10,10),  vimge(2,10,10),
     &                 vreuep(2,10,10), vimuep(2,10,10),
     &                 vredep(2,10,10), vimdep(2,10,10),
     &                 vresep(2,10,10), vimsep(2,10,10),
     &                 vregep(2,10,10), vimgep(2,10,10)

      common /vectors/ vreu,   vimu, 
     &                 vred,   vimd,
     &                 vres,   vims,
     &                 vreg,   vimg,
     &                 vreup,  vimup, 
     &                 vredp,  vimdp,
     &                 vresp,  vimsp,
     &                 vregp,  vimgp,
     &                 vreue,  vimue,
     &                 vrede,  vimde,
     &                 vrese,  vimse,
     &                 vrege,  vimge,
     &                 vreuep, vimuep,
     &                 vredep, vimdep,
     &                 vresep, vimsep,
     &                 vregep, vimgep


      character name*15
      character dir*4
      character name2*19

      logical first
      data first /.true./
      save first

      if (first) then
         first = .false.
         print*,'read vector really'
      else
c         print*,'vector already read'
         return
      endif
  


      dir='grv/'
c      dir='mrs/'

*
* read leading order data files
*
      name   = 'ldamp.dat'
      name2   = dir//name
      call readvect(name2,ldamp_1_dat,ldamp_2_dat)

c      print*,ldamp_1_dat
c      print*,ldamp_2_dat

c      return

      name   = 'ldampe.dat'
      name2   = dir//name
      call readvect(name2,ldampe_1_dat,ldampe_2_dat)

      name   = 'ldamppol.dat'
      name2   = dir//name
      call readvect(name2,ldamppol_1_dat,ldamppol_2_dat)

      name   = 'ldamppole.dat'
      name2   = dir//name
      call readvect(name2,ldamppole_1_dat,ldamppole_2_dat)

      name   = 'lgamp.dat'
      name2   = dir//name
      call readvect(name2,lgamp_1_dat,lgamp_2_dat)

      name   = 'lgampe.dat'
      name2   = dir//name
      call readvect(name2,lgampe_1_dat,lgampe_2_dat)

      name   = 'lgamppol.dat'
      name2   = dir//name
      call readvect(name2,lgamppol_1_dat,lgamppol_2_dat)

      name   = 'lgamppole.dat'
      name2   = dir//name
      call readvect(name2,lgamppole_1_dat,lgamppole_2_dat)

      name   = 'lsamp.dat'
      name2   = dir//name
      call readvect(name2,lsamp_1_dat,lsamp_2_dat)

      name   = 'lsampe.dat'
      name2   = dir//name
      call readvect(name2,lsampe_1_dat,lsampe_2_dat)

      name   = 'lsamppol.dat'
      name2   = dir//name
      call readvect(name2,lsamppol_1_dat,lsamppol_2_dat)

      name   = 'lsamppole.dat'
      name2   = dir//name
      call readvect(name2,lsamppole_1_dat,lsamppole_2_dat)

      name   = 'luamp.dat'
      name2   = dir//name
      call readvect(name2,luamp_1_dat,luamp_2_dat)

      name   = 'luampe.dat'
      name2   = dir//name
      call readvect(name2,luampe_1_dat,luampe_2_dat)

      name   = 'luamppol.dat'
      name2   = dir//name
      call readvect(name2,luamppol_1_dat,luamppol_2_dat)

      name   = 'luamppole.dat'
      name2   = dir//name
      call readvect(name2,luamppole_1_dat,luamppole_2_dat)

*
* read nlo amplitudes
*

      name   = 'nlodamp.dat'
      name2   = dir//name
      call readvect(name2,nlodamp_1_dat,nlodamp_2_dat)

      name   = 'nlodampe.dat'
      name2   = dir//name
      call readvect(name2,nlodampe_1_dat,nlodampe_2_dat)

      name   = 'nlodamppol.dat'
      name2   = dir//name
      call readvect(name2,nlodamppol_1_dat,nlodamppol_2_dat)

      name   = 'nlodamppole.dat'
      name2   = dir//name
      call readvect(name2,nlodamppole_1_dat,nlodamppole_2_dat)

      name   = 'nlogamp.dat'
      name2   = dir//name
      call readvect(name2,nlogamp_1_dat,nlogamp_2_dat)

      name   = 'nlogampe.dat'
      name2   = dir//name
      call readvect(name2,nlogampe_1_dat,nlogampe_2_dat)

      name   = 'nlogamppol.dat'
      name2   = dir//name
      call readvect(name2,nlogamppol_1_dat,nlogamppol_2_dat)

      name   = 'nlogamppole.dat'
      name2   = dir//name
      call readvect(name2,nlogamppole_1_dat,nlogamppole_2_dat)

      name   = 'nlosamp.dat'
      name2   = dir//name
      call readvect(name2,nlosamp_1_dat,nlosamp_2_dat)

      name   = 'nlosampe.dat'
      name2   = dir//name
      call readvect(name2,nlosampe_1_dat,nlosampe_2_dat)

      name   = 'nlosamppol.dat'
      name2   = dir//name
      call readvect(name2,nlosamppol_1_dat,nlosamppol_2_dat)

      name   = 'nlosamppole.dat'
      name2   = dir//name
      call readvect(name2,nlosamppole_1_dat,nlosamppole_2_dat)

      name   = 'nlouamp.dat'
      name2   = dir//name
      call readvect(name2,nlouamp_1_dat,nlouamp_2_dat)

      name   = 'nlouampe.dat'
      name2   = dir//name
      call readvect(name2,nlouampe_1_dat,nlouampe_2_dat)

      name   = 'nlouamppol.dat'
      name2   = dir//name
      call readvect(name2,nlouamppol_1_dat,nlouamppol_2_dat)

      name   = 'nlouamppole.dat'
      name2   = dir//name
      call readvect(name2,nlouamppole_1_dat,nlouamppole_2_dat)

*
* fill amplitude vectors
*
      do ic = 1,10
         do jc = 1,10

            vreu(1,ic,jc) = luamp_2_dat(ic,jc,1)
            vreu(2,ic,jc) = nlouamp_2_dat(ic,jc,1)
            vimu(1,ic,jc) = luamp_2_dat(ic,jc,2)
            vimu(2,ic,jc) = nlouamp_2_dat(ic,jc,2)

            vred(1,ic,jc) = ldamp_2_dat(ic,jc,1)
            vred(2,ic,jc) = nlodamp_2_dat(ic,jc,1)
            vimd(1,ic,jc) = ldamp_2_dat(ic,jc,2)
            vimd(2,ic,jc) = nlodamp_2_dat(ic,jc,2)

            vres(1,ic,jc) = lsamp_2_dat(ic,jc,1)
            vres(2,ic,jc) = nlosamp_2_dat(ic,jc,1)
            vims(1,ic,jc) = lsamp_2_dat(ic,jc,2)
            vims(2,ic,jc) = nlosamp_2_dat(ic,jc,2)

            vreg(1,ic,jc) = lgamp_2_dat(ic,jc,1)
            vreg(2,ic,jc) = nlogamp_2_dat(ic,jc,1)
            vimg(1,ic,jc) = lgamp_2_dat(ic,jc,2)
            vimg(2,ic,jc) = nlogamp_2_dat(ic,jc,2)

            vreup(1,ic,jc) = luamppol_2_dat(ic,jc,1)
            vreup(2,ic,jc) = nlouamppol_2_dat(ic,jc,1)
            vimup(1,ic,jc) = luamppol_2_dat(ic,jc,2)
            vimup(2,ic,jc) = nlouamppol_2_dat(ic,jc,2)

            vredp(1,ic,jc) = ldamppol_2_dat(ic,jc,1)
            vredp(2,ic,jc) = nlodamppol_2_dat(ic,jc,1)
            vimdp(1,ic,jc) = ldamppol_2_dat(ic,jc,2)
            vimdp(2,ic,jc) = nlodamppol_2_dat(ic,jc,2)

            vresp(1,ic,jc) = lsamppol_2_dat(ic,jc,1)
            vresp(2,ic,jc) = nlosamppol_2_dat(ic,jc,1)
            vimsp(1,ic,jc) = lsamppol_2_dat(ic,jc,2)
            vimsp(2,ic,jc) = nlosamppol_2_dat(ic,jc,2)

            vregp(1,ic,jc) = lgamppol_2_dat(ic,jc,1)
            vregp(2,ic,jc) = nlogamppol_2_dat(ic,jc,1)
            vimgp(1,ic,jc) = lgamppol_2_dat(ic,jc,2)
            vimgp(2,ic,jc) = nlogamppol_2_dat(ic,jc,2)

            vreue(1,ic,jc) = luampe_2_dat(ic,jc,1)
            vreue(2,ic,jc) = nlouampe_2_dat(ic,jc,1)
            vimue(1,ic,jc) = luampe_2_dat(ic,jc,2)
            vimue(2,ic,jc) = nlouampe_2_dat(ic,jc,2)

            vrede(1,ic,jc) = ldampe_2_dat(ic,jc,1)
            vrede(2,ic,jc) = nlodampe_2_dat(ic,jc,1)
            vimde(1,ic,jc) = ldampe_2_dat(ic,jc,2)
            vimde(2,ic,jc) = nlodampe_2_dat(ic,jc,2)

            vrese(1,ic,jc) = lsampe_2_dat(ic,jc,1)
            vrese(2,ic,jc) = nlosampe_2_dat(ic,jc,1)
            vimse(1,ic,jc) = lsampe_2_dat(ic,jc,2)
            vimse(2,ic,jc) = nlosampe_2_dat(ic,jc,2)

            vrege(1,ic,jc) = lgampe_2_dat(ic,jc,1)
            vrege(2,ic,jc) = nlogampe_2_dat(ic,jc,1)
            vimge(1,ic,jc) = lgampe_2_dat(ic,jc,2)
            vimge(2,ic,jc) = nlogampe_2_dat(ic,jc,2)

            vreuep(1,ic,jc) = luamppole_2_dat(ic,jc,1)
            vreuep(2,ic,jc) = nlouamppole_2_dat(ic,jc,1)
            vimuep(1,ic,jc) = luamppole_2_dat(ic,jc,2)
            vimuep(2,ic,jc) = nlouamppole_2_dat(ic,jc,2)

            vredep(1,ic,jc) = ldamppole_2_dat(ic,jc,1)
            vredep(2,ic,jc) = nlodamppole_2_dat(ic,jc,1)
            vimdep(1,ic,jc) = ldamppole_2_dat(ic,jc,2)
            vimdep(2,ic,jc) = nlodamppole_2_dat(ic,jc,2)

            vresep(1,ic,jc) = lsamppole_2_dat(ic,jc,1)
            vresep(2,ic,jc) = nlosamppole_2_dat(ic,jc,1)
            vimsep(1,ic,jc) = lsamppole_2_dat(ic,jc,2)
            vimsep(2,ic,jc) = nlosamppole_2_dat(ic,jc,2)

            vregep(1,ic,jc) = lgamppole_2_dat(ic,jc,1)
            vregep(2,ic,jc) = nlogamppole_2_dat(ic,jc,1)
            vimgep(1,ic,jc) = lgamppole_2_dat(ic,jc,2)
            vimgep(2,ic,jc) = nlogamppole_2_dat(ic,jc,2)

         enddo
      enddo




c      do jc=1,10
c         print*,vregep(1,1,jc),vimgep(1,1,jc)
c      enddo

      return
      end



      include 'readvect.f'
