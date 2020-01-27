****************************************************************
*                                                              *      
*                         REU                                  * 
*                                                              *
****************************************************************

      double precision function reu(ord,xi,q2i)

      implicit none

      include 'func.include'      
      integer i,j,k

c      print*,'this is reu'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif

      amp = interpol(ord,vreu,x,q2,xi,q2i)

      reu = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMU                                  * 
*                                                              *
****************************************************************

      double precision function imu(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimu,x,q2,xi,q2i)
      
      imu = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         RED                                  * 
*                                                              *
****************************************************************

      double precision function red(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vred,x,q2,xi,q2i)
      
      red = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMD                                  * 
*                                                              *
****************************************************************

      double precision function imd(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimd,x,q2,xi,q2i)
      
      imd = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         RES                                  * 
*                                                              *
****************************************************************

      double precision function res(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vres,x,q2,xi,q2i)
      
      res = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMS                                  * 
*                                                              *
****************************************************************

      double precision function ims(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vims,x,q2,xi,q2i)
      
      ims = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         REG                                  * 
*                                                              *
****************************************************************

      double precision function reg(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vreg,x,q2,xi,q2i)
      
      reg = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMG                                  * 
*                                                              *
****************************************************************

      double precision function img(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimg,x,q2,xi,q2i)
      
      img = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         REUP                                 * 
*                                                              *
****************************************************************

      double precision function reup(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vreup,x,q2,xi,q2i)
      
      reup = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMUP                                 * 
*                                                              *
****************************************************************

      double precision function imup(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimup,x,q2,xi,q2i)
      
      imup = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         REDP                                 * 
*                                                              *
****************************************************************

      double precision function redp(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vredp,x,q2,xi,q2i)
      
      redp = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMDP                                 * 
*                                                              *
****************************************************************

      double precision function imdp(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimdp,x,q2,xi,q2i)
      
      imdp = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         RESP                                 * 
*                                                              *
****************************************************************

      double precision function resp(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vresp,x,q2,xi,q2i)
      
      resp = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMSP                                 * 
*                                                              *
****************************************************************

      double precision function imsp(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimsp,x,q2,xi,q2i)
      
      imsp = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         REGP                                 * 
*                                                              *
****************************************************************

      double precision function regp(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vregp,x,q2,xi,q2i)
      
      regp = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMGP                                 * 
*                                                              *
****************************************************************

      double precision function imgp(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimgp,x,q2,xi,q2i)
      
      imgp = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         REUE                                 * 
*                                                              *
****************************************************************

      double precision function reue(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vreue,x,q2,xi,q2i)
      
      reue = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMUE                                 * 
*                                                              *
****************************************************************

      double precision function imue(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimue,x,q2,xi,q2i)
      
      imue = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         REDE                                 * 
*                                                              *
****************************************************************

      double precision function rede(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vrede,x,q2,xi,q2i)
      
      rede = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMDE                                 * 
*                                                              *
****************************************************************

      double precision function imde(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimde,x,q2,xi,q2i)
      
      imde = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         RESE                                 * 
*                                                              *
****************************************************************

      double precision function rese(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vrese,x,q2,xi,q2i)
      
      rese = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMSE                                 * 
*                                                              *
****************************************************************

      double precision function imse(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimse,x,q2,xi,q2i)
      
      imse = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         REGE                                 * 
*                                                              *
****************************************************************

      double precision function rege(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vrege,x,q2,xi,q2i)
      
      rege = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMGE                                 * 
*                                                              *
****************************************************************

      double precision function imge(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimge,x,q2,xi,q2i)
      
      imge = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         REUEP                                * 
*                                                              *
****************************************************************

      double precision function reuep(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vreuep,x,q2,xi,q2i)
      
      reuep = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMUEP                                * 
*                                                              *
****************************************************************

      double precision function imuep(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimuep,x,q2,xi,q2i)
      
      imuep = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         REDEP                                * 
*                                                              *
****************************************************************

      double precision function redep(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vredep,x,q2,xi,q2i)
      
      redep = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMDEP                                * 
*                                                              *
****************************************************************

      double precision function imdep(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimdep,x,q2,xi,q2i)
      
      imdep = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         RESEP                                * 
*                                                              *
****************************************************************

      double precision function resep(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vresep,x,q2,xi,q2i)
      
      resep = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMSEP                                * 
*                                                              *
****************************************************************

      double precision function imsep(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimsep,x,q2,xi,q2i)
      
      imsep = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         REGEP                                * 
*                                                              *
****************************************************************

      double precision function regep(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vregep,x,q2,xi,q2i)
      
      regep = dble(amp)

      return
      end

****************************************************************
*                                                              *      
*                         IMGEP                                * 
*                                                              *
****************************************************************

      double precision function imgep(ord,xi,q2i)
      implicit none

      include 'func.include'

      if (first) then
         first = .false.
         call readnlo
c         print*,'read vector'
      endif
      
      amp = interpol(ord,vimgep,x,q2,xi,q2i)
      
      imgep = dble(amp)

      return
      end


      include 'interpol.f'
      include 'readnlo.f'



