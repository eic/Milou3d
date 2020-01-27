*CMZ :          28/11/2003  09.31.39  by  Unknown
*-- Author :    Unknown   28/11/2003


*==========================================
      SUBROUTINE ROTATE(AR,BR,iflag,ARNEW)
*==========================================

* iflag  = -1 :
* Transforms vector AR given in a system
* axis where z = BR(3), back to the
* standard system axis

* iflag = +1 :
* Transforms vector AR given in standard system
* into a system axis where z = BR(3)


      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION AR(4),BR(3)
      dimension arnew(4)

      DIMENSION ROT(3,3)
      DIMENSION APR(3)
      dimension avt(3)
      dimension aa(3)

      dimension zz(3)
      logical do_xaxis
      common/xaxis/do_xaxis,i_xaxis

      zz(1) = 0.0d0
      zz(2) = 0.0d0
      zz(3) = 1.0d0
      psbefore=prodsca3norm(ar,zz)

      do i=1,3
       aa(i) = ar(i)
      enddo

      ANORM = AR(1)**2+AR(2)**2+AR(3)**2
      ANORM = DSQRT(ANORM)


      BNORM = BR(1)**2+BR(2)**2+BR(3)**2
      BNORM = DSQRT(BNORM)


      COSTHETA = BR(3)/BNORM

      SIN2THETA = 1.0d0 - COSTHETA**2
      SINTHETA = DSQRT(SIN2THETA)

      SINPHI = BR(1) / (SINTHETA * BNORM)
      COSPHI = BR(2) / (SINTHETA * BNORM)

c      write(6,*) 'costheta  cosphi  ',costheta, cosphi

      if (iflag.eq.-1) then

      ROT(1,1) = COSPHI
      ROT(1,2) = -SINPHI
      ROT(1,3) = 0.0D0

      ROT(2,1) = COSTHETA * SINPHI
      ROT(2,2) = COSTHETA * COSPHI
      ROT(2,3) = -SINTHETA

      ROT(3,1) = SINTHETA * SINPHI
      ROT(3,2) = SINTHETA * COSPHI
      ROT(3,3) = COSTHETA

      elseif (iflag.eq.1) then

      ROT(1,1) = COSPHI
      ROT(2,1) = -SINPHI
      ROT(3,1) = 0.0D0

      ROT(1,2) = COSTHETA * SINPHI
      ROT(2,2) = COSTHETA * COSPHI
      ROT(3,2) = -SINTHETA

      ROT(1,3) = SINTHETA * SINPHI
      ROT(2,3) = SINTHETA * COSPHI
      ROT(3,3) = COSTHETA

      endif


      do i=1,3
        apr(i) = 0.0d0
        avt(i) = ar(i)
      enddo

      if (iflag.eq.-1) then
       avt(1) = ar(2)
       avt(2) = -ar(1)
       avt(1) = avt(1) * i_xaxis
       avt(2) = avt(2) * i_xaxis
      endif

      do i=1,3
        do j=1,3
        APR(I) = apr(i) +  Avt(J) * ROT(J,I)
        enddo
      ENDDO

      apnorm = apr(1)**2+apr(2)**2+apr(3)**2
      apnorm=dsqrt(apnorm)



      do i=1,3
       arnew(i) = apr(i)
      enddo

* --- rotation angles compatible with Belitsky paper :
      if (iflag.eq.1) then
       arnew(1) = -apr(2)
       arnew(2) = apr(1)
      endif

* --- if we're rotating the incoming electron vector,
*     determine the direction of the x-axis
      if (do_xaxis) then
       i_xaxis = +1
       if (arnew(1).lt.0) i_xaxis = -1
      endif

      if (iflag.eq.1.and.i_xaxis.eq.-1) then
       arnew(1) = -arnew(1)
       arnew(2) = -arnew(2)
      endif
* ---

      arnew(4) = ar(4)

      psafter=prodsca3norm(ar,br)


      RETURN
      END

*----------------------------------------------------
!      real function prodsca3norm*8(a,b)
      double precision function prodsca3norm(a,b)
*----------------------------------------------------

*----------------------------------------------------
*     scalar product with the metrix (1,1,1,0)
*     normalized
*----------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)

      dimension a(3),b(3)

      prodsca3norm = -999.
      prodsca3 = 0.
     &    +a(1)*b(1)
     &    +a(2)*b(2)
     &    +a(3)*b(3)
      pn1 = dsqrt(a(1)**2+a(2)**2+a(3)**2)
      pn2 = dsqrt(b(1)**2+b(2)**2+b(3)**2)
      if (pn1*pn2.gt.0) then
      prodsca3norm = prodsca3 / (pn1*pn2)
      else
       write(6,*) 'prodsca3norm : pn1*pn2 = 0 !'
      endif

      return
      end


C----------------------------------------------------------------------
      SUBROUTINE BOOST (PV,BST)
      IMPLICIT NONE
C----------------------------------------------------------------------
C-
C-   PV --> PV' (BOOSTED BY BST)
C-
C-   NB.  PV   (PX,PY,PZ,E,M)
C-        BST  (BX,BY,BZ)
C-
C
CIMP      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PV(4),BST(3)
C
      real*8 bb2,gam,tm1,tm2
      integer j
      BB2 = BST(1)**2 + BST(2)**2 + BST(3)**2

      if(BB2.ge.1.d0)then
       write(1,*)' warning in boostsusy ',BB2
       bb2=0.99999d0
      endif

      GAM = 1. / DSQRT (1. - BB2)
      TM1 = BST(1) * PV(1) + BST(2) * PV(2) + BST(3) * PV(3)
      TM2 = GAM * TM1 / (GAM + 1.)  + PV(4)
      PV(4) = GAM * (PV(4) + TM1)
      DO 10 J=1,3
   10 PV(J) = PV(J) + GAM * BST(J) * TM2
      END


c -------------------------------------------------------
      subroutine v_print(vec)
c -------------------------------------------------------
      real*8 vec(4)
      real*8 mass2, mass
      mass2 = vec(4)**2-vec(1)**2-vec(2)**2-vec(3)**2
      if (mass2.gt.0) then
       mass = dsqrt(mass2)
      else
       mass = - dsqrt(-mass2)
      endif
      write(6,*) 'vec = ',vec(1),vec(2),vec(3),vec(4),mass
      return
      end


c -------------------------------------------------------
      subroutine v_print4(vec)
c -------------------------------------------------------
      real vec(4)
      real mass2, mass
      mass2 = vec(4)**2-vec(1)**2-vec(2)**2-vec(3)**2
      if (mass2.gt.0) then
       mass = sqrt(mass2)
      else
       mass = - sqrt(-mass2)
      endif
      write(6,*) 'vec = ',vec(1),vec(2),vec(3),vec(4),mass
      return
      end


c -------------------------------------------------------
      subroutine v_copy(a,b)
c -------------------------------------------------------
      real*8 a(4)
      real   b(4)
      do i=1,4
       b(i) = sngl(a(i))
      enddo
      return
      end


c -------------------------------------------------------
      subroutine v_copy8(a,b)
c -------------------------------------------------------
      real*8 a(4)
      real*8 b(4)
      do i=1,4
       b(i) = a(i)
      enddo
      return
      end


c -------------------------------------------------------
      subroutine v_copy48(a,b)
c -------------------------------------------------------
      real*4 a(4)
      real*8 b(4)
      do i=1,4
       b(i) = dble(a(i))
      enddo
      return
      end


c -------------------------------------------------------
      subroutine check(li,pi,lo,po,rg)
c -------------------------------------------------------
      real*4 li(4),pi(4),lo(4),po(4),rg(4)
      px_in=0
      py_in =0
      pz_in=0
      e_in=0
      do i=1,4
       px_in = li(1) + pi(1)
       py_in = li(2) + pi(2)
       pz_in = li(3) + pi(3)
       e_in  = li(4) + pi(4)
      enddo
      px_out = 0
      py_out = 0
      pz_out = 0
      e_out = 0
      do i=1,4
       px_out = lo(1) + po(1) + rg(1)
       py_out = lo(2) + po(2) + rg(2)
       pz_out = lo(3) + po(3) + rg(3)
       e_out  = lo(4) + po(4) + rg(4)
      enddo
      write(6,*) ' --- check kinematics --- '
      write(6,*) 'In : ',px_in,py_in,pz_in,e_in
      write(6,*) 'Out : ',px_out,py_out,pz_out,e_out
      return
      end




