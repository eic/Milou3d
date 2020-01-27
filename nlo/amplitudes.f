      subroutine amplitudes
      implicit none

*
* This routine files tables with amplitudes
*
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

      real*8 reu,imu,red,imd,res,ims,reg,img,reup,imup, 
     >     redp,imdp,resp,imsp,regp,imgp,reue,imue,rede,imde,rese,imse
     >     ,rege,imge,reuep,imuep,redep,imdep,resep,imsep,regep,imgep


      character name*17

      double precision x, q2
      double precision xg(10), q2g(10)
      integer i, j, k
      integer ord
      data xg / 0.000288, 0.000355, 0.0005, 0.00075, 0.001,
     &          0.002,    0.00355,  0.005,  0.006,   0.0072 /

      data q2g /  1.9993960, 3.0765160, 4.3848360, 5.9194890,
     &            7.6895290, 9.6907690, 11.923209, 14.386849,
     &           17.073424, 19.998784/ 

      do i = 1,2
         do j= 1, nbins
            do k = 1, nbins
               ord = i
               x   = dble(10**(log10(xg(1))
     &              + real(j)*(log10(xg(10))-log10(xg(1)))/real(nbins)))
               q2  = dble(q2g(1)+real(k)*(q2g(10)-q2g(1))/real(nbins))

c               print*,ord,real(dlog10(x)),real(x),real(q2)

               vectreu(i,j,k)   = real(reu(ord,x,q2))
*               print*,i,j,k,vectreu(i,j,k)

               vectimu(i,j,k)   = real(imu(ord,x,q2))
               vectred(i,j,k)   = real(red(ord,x,q2))
               vectimd(i,j,k)   = real(imd(ord,x,q2))
               vectres(i,j,k)   = real(res(ord,x,q2))
               vectims(i,j,k)   = real(ims(ord,x,q2))
               vectreg(i,j,k)   = real(reg(ord,x,q2))
               vectimg(i,j,k)   = real(img(ord,x,q2))
               vectreup(i,j,k)  = real(reup(ord,x,q2))
               vectimup(i,j,k)  = real(imup(ord,x,q2))
               vectredp(i,j,k)  = real(redp(ord,x,q2))
               vectimdp(i,j,k)  = real(imdp(ord,x,q2))
               vectresp(i,j,k)  = real(resp(ord,x,q2))
               vectimsp(i,j,k)  = real(imsp(ord,x,q2))
               vectregp(i,j,k)  = real(regp(ord,x,q2))
               vectimgp(i,j,k)  = real(imgp(ord,x,q2))
               vectreue(i,j,k)  = real(reue(ord,x,q2))
               vectimue(i,j,k)  = real(imue(ord,x,q2))
               vectrede(i,j,k)  = real(rede(ord,x,q2))
               vectimde(i,j,k)  = real(imde(ord,x,q2))
               vectrese(i,j,k)  = real(rese(ord,x,q2))
               vectimse(i,j,k)  = real(imse(ord,x,q2))
               vectrege(i,j,k)  = real(rege(ord,x,q2))
               vectimge(i,j,k)  = real(imge(ord,x,q2))
               vectreuep(i,j,k) = real(reuep(ord,x,q2))
               vectimuep(i,j,k) = real(imuep(ord,x,q2))
               vectredep(i,j,k) = real(redep(ord,x,q2))
               vectimdep(i,j,k) = real(imdep(ord,x,q2))
               vectresep(i,j,k) = real(resep(ord,x,q2))
               vectimsep(i,j,k) = real(imsep(ord,x,q2))
               vectregep(i,j,k) = real(regep(ord,x,q2))
               vectimgep(i,j,k) = real(imgep(ord,x,q2))
            enddo
         enddo
      enddo


 101  FORMAT(F16.8)


      name  = 'txt/vectreu.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectreu(i,j,k)
*               print*,vectreu(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectimu.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimu(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectred.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectred(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectimd.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimd(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)
      

      name  = 'txt/vectres.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectres(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectims.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectims(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectreg.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectreg(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectimg.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimg(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectreup.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectreup(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectimup.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimup(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectredp.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectredp(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectimdp.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimdp(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectresp.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectresp(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectimsp.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimsp(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectregp.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectregp(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)
      

      name  = 'txt/vectimgp.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimgp(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectreue.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectreue(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectimue.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimue(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectrede.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectrede(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)
      

      name  = 'txt/vectimde.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimde(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectrese.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectrese(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectimse.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimse(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectrege.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectrege(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectimge.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimge(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectreuep.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectreuep(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectimuep.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimuep(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectredep.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectredep(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectimdep.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimdep(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectresep.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectresep(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectimsep.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimsep(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectregep.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectregep(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)


      name  = 'txt/vectimgep.txt'
      OPEN (UNIT=95,FILE=name,STATUS='UNKNOWN')
      do i= 1,2
         do j=1,nbins
            do k = 1,nbins
               write(95,101) vectimgep(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(95)

      return
      end



