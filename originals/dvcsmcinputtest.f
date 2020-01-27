c      program test

      implicit double precision (a-h,o-z)

C
C     the first common block: iord = 1 (LO), 2 (NLO)
C     the second common block: icountt = 1 (t changed), 2 (t did not change)
C                              icountxq = 1 (x or Q^2 changed, do interpolation), 
C                                         2 (x and Q^2 unchanged, do not interpolate)
C     third common block: set = 1 (most general setting with all polarizations etc.,
C                               2 (all angles integrated out)
C

     
      common /order/ iord

      common /counter1/ icountxq,icountt

      common /setting/ nset
C
C	icount = 1 (read in data files, only on first call of subroutine), 2 (do not read-in on subsequent calls)
C	x = x_bj, q = Q^2, s = center of mass energy, t as in t 
C       phi = azimuthal angle between lepton and final state scattering planes
C       pphi = angle between lepton plane and transverse polarization vector (see Belitsky, Mueller and Kirchner), 
C       theta = angle between long. and trans. polarization, 
C       tchar = +1 positron, -1 electron,
C	xlaml = polarization of lepton, +1 in direction of movement, -1 against direction of movement,
C	xlamp = target polarization, +1(-1) = polarized along(opposite) probe/target beam direction,
C                                     0 = unpolarized
C       Z = charge of nucleus, A= Atomic number, spin = spin of target (0=0, 1=1/2, 2=1 etc.)
C


C
C     data arrays to store results for observables prior to integration over t.
C
C
      parameter(mx=300)
      dimension f(mx),f1(mx),res(100),f2(mx),f3(mx)
      dimension f4(mx),f5(mx),f6(mx),f7(mx),f8(mx),f9(mx),resut(2,mx)
      dimension f10(mx),f11(mx)
      character*2 temp,temp0
      character*78 temp1,temp2
      external SMPSNA
C
C     # of points in x and Q^2 you want to analyze
C
       mx1 = 15
       mq = 12

C
C     initalize value for dvcs subroutine
C
      spin = 1D0
      Z = 1D0
      A = 1D0
C
C     # of points in x and Q^2
C
      nx = 58
      nq = 40
C
C     initializing flags 
C
      icount = 1
      icountxq = 1
      icountt = 1
      nset = 2
C
C     initial x and Q^2 value
C
      x = 0.16D0
      q = 1.0D0
C
C     setting the energy s
C


C An EIC setting (electron = 5 GeV, target = 200 GeV)
c      s = 4.*5.*200.
C HERMES setting
c      s = 2.*27.6*0.938272D0
C HERA setting with 920 GEV protons
c      s = 4.*27.6*920.
C CLAS setting: CLAS I = 4.3 GEV, CLAS II = 5.75 GeV
      s = 2.*4.3D0*0.938272D0
C
C     dummy settings for the angles at the moment
C
      phi = 0.1*3.141592653589793
      pphi = 0.0
      theta = 0.0
C
C     charge, polarization and order in alpha_s
C
      tchar = -1D0
      xlaml = 1D0
      xlamp = 0D0
      iord = 1
C
C     initial value of x, needed to restore later on x loop
C
      x1 = 0.16D0
 
C
C     open files for x vs. Q^2 printout
C
      do n4 = 1,mx1

      write(temp0,333) n4+10 

      temp1 = 'resulttintxvsq'//temp0//'.dat'

      n3 = n4+50
      
      open(unit=n3,file=temp1,status='unknown')
      

      enddo

C
C     Q^2 loop
C
      do n1 = 1,mq

C
C     open file for printout Q^2 vs. x
C
      write(temp,333) n1 

      temp2 = 'resulttintqvsx'//temp//'.dat'

      open(unit=n1,file=temp2,status='unknown')

C
C     x loop
C
      do n2 = 1,mx1

C
C     kinematic variables needed to check wehther chosen values of x, Q^2
C     are in or out of physical region
C
      ep1 = (2.*x*0.938272D0)**2D0/q

      tmin = - q*(2.*(1.-x/A)*(1. - dsqrt(1+ep1)) + ep1)/
     >           (4.*(x/A)*(1.-(x/A)) + ep1)

C
C     starting value for t
C
      t = -0.0001D0
      if (t.gt.tmin) then
         t = 1.01*tmin
         endif
C
C     stepsize in t (necessary for simpson integration below)
C
      dx = (0.25D0+tmin)/mx

C
C     maximum value of t to which to integrate
C      
      tmax = -0.25D0

C
C     LO/NLO loop
C
      do it = 1,2
C
C     loop in t
C
      do n = 1,mx 

C
C      kinematic variables needed to check wehther chosen values of x, Q^2
C     are in or out of physical region
C
         ym = (q+t)/(q+x*t)

         y1 = q/(x*(s - (0.938272D0)**2)) 
         

C
C     if y > y_max, skip to next value in x
C
       if (y1.ge.0.999D0*ym) goto 999

C
C     setting counter to be used later in dvcs subroutine
C
         if (n.gt.1) then

            icount = 2

            endif
C
C     call of DVCS subroutine
C         
            
      call dvcsob(spin,Z,A,nx,nq,icount,x,q,s,t,phi,pphi,
     >     theta,tchar,xlaml,xlamp,iord,res)
c      print *, t,res(1),res(2),res(3),res(4),res(5)/res(6),res(9)
c     > /res(10),res(5),res(9)
c      read (*,*)
C
C     storage of results for CA, SSA and total cross section with and without twist-3
C
      f(n) = res(5)

      f1(n) = res(6)

      f2(n) = res(7)

      f3(n) = res(8)

      f4(n) = res(9)

      f5(n) = res(10)

      f6(n) = res(11)

      f7(n) = res(12)
      
      f8(n) = res(1)

      f9(n) = res(13)

      f10(n) = res(2)

      f11(n) = res(4)
C
C     goto next value in t
C
      t = t - dx

C
C     set counter for interpolation in x and Q^2
C

      icountxq = 2

      icountt = 1
      
      x = x*A

      enddo 
C
C     integration in t
C

      result0 = SMPSNA (mx, dx, f10, ERR)
      result01 = SMPSNA (mx, dx, f11, ERR)
      resut(iord,7) = result0/result01
      result1 = SMPSNA (mx, dx, f, ERR)
      result2 = SMPSNA (mx, dx, f1, ERR1)
      result3 = result1/result2
      resut(iord,1) = result1/result2
      result4 = SMPSNA (mx, dx, f2, ERR2)
      result5 = SMPSNA (mx, dx, f3, ERR3)
      result6 = result4/result5
      resut(iord,2) = result4/result5
      result7 = SMPSNA (mx, dx, f4, ERR4)
      result8 = SMPSNA (mx, dx, f5, ERR5)
      result9 = result7/result8
      resut(iord,3) = result7/result8
      result10 = SMPSNA (mx, dx, f6, ERR6)
      result11 = SMPSNA (mx, dx, f7, ERR7)
      result12 = result10/result11
      resut(iord,4) = result10/result11
      resut(iord,5) = SMPSNA (mx, dx, f8, ERR8)
      resut(iord,6) = SMPSNA (mx, dx, f9, ERR9)
      print *, t, x, q,iord
      print *," SSA wo tw-3=",result12," SSA w tw-3=",result9
      print *," CA wo tw-3=",result6," CA w tw-3=",result3
      print *," Totsig wo tw-3=",resut(iord,6)," Totsig w tw-3="
     > ,resut(iord,5)

C
C     reset t to starting value
C      

      t = -0.0001D0
      if (t.gt.tmin) then
         t = 1.01*tmin
         endif
      
C
C     after LO (iord=1) now do NLO (iord=2)
C
      iord = 2
C
C     reset internal counters for dvcs subroutine
C
      icount = 1
      icountxq = 1
      icountt = 1

      enddo

C
C     summary: t,x,Q^2,CA-w-tw3LO,CA-wo-tw3LO,CA-w-tw3LO/CA-wo-tw3LO,SSA-w-tw3LO,SSA-wo-tw3LO
C     ,SSA-w-tw3LO/SSA-wo-tw3LO,totsig-w-tw3LO,totsig-wo-tw3LO,totsig-w-tw3LO/totsig-wo-tw3LO,
C     CA-w-tw3NLO,CA-wo-tw3NLO,CA-w-tw3NLO/CA-wo-tw3NLO,SSA-w-tw3NLO,SSA-wo-tw3NLO
C     ,SSA-w-tw3NLO/SSA-wo-tw3NLO,totsig-w-tw3NLO,totsig-wo-tw3NLO,totsig-w-tw3NLO/totsig-wo-tw3NLO,
C     CA-wo-tw3LO/CA-wo-tw3NLO,SSA-wo-tw3LO/SSA-wo-tw3NLO,totsig-wo-tw3LO/totsig-wo-tw3NLO
C

C
C     writeout for Q^2 vs. x
C

      write (n1,102) tmax,x,q,resut(1,1),resut(1,2)
     >,resut(1,1)/resut(1,2)
     >,resut(1,3),resut(1,4),resut(1,3)/resut(1,4),resut(1,5),resut(1,6)
     >,resut(1,5)/resut(1,6),resut(2,1),resut(2,2),resut(2,1)/resut(2,2)
     >,resut(2,3),resut(2,4),resut(2,3)/resut(2,4),resut(2,5),resut(2,6)
     >,resut(2,5)/resut(2,6),resut(1,2)/resut(2,2),resut(1,4)/resut(2,4)
     >,resut(1,6)/resut(2,6),resut(1,7),resut(2,7)

C
C     writeout for x vs. Q^2
C
      n7 = n2+50
      
      write (n7,102) tmax,x,q,resut(1,1),resut(1,2)
     >,resut(1,1)/resut(1,2)
     >,resut(1,3),resut(1,4),resut(1,3)/resut(1,4),resut(1,5),resut(1,6)
     >,resut(1,5)/resut(1,6),resut(2,1),resut(2,2),resut(2,1)/resut(2,2)
     >,resut(2,3),resut(2,4),resut(2,3)/resut(2,4),resut(2,5),resut(2,6)
     >,resut(2,5)/resut(2,6),resut(1,2)/resut(2,2),resut(1,4)/resut(2,4)
     >,resut(1,6)/resut(2,6),resut(1,7),resut(2,7)
      

 
 999  iord = 1

      if (x.lt.0.001D0) then

         x = x + 0.0002

         elseif (x.ge.0.001D0.and.x.lt.0.01D0) then

            x = x + 0.002

            elseif (x.ge.0.01D0) then

               x = x + 0.02D0

               endif

      enddo

      close (unit=n1)

      x = x1

      if (s.lt.190D0) then

      q = q + 0.25D0

      else

         q = q + 2

         endif

      enddo

      do n5 = mx1+50,51, -1

         close(unit=n5)

         enddo

 102  FORMAT(26(E12.5,1X))
 333  FORMAT(I2.2)

      end

      FUNCTION SMPSNA (NX, DX, F, ERR)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MAXX = 1000)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DIMENSION F(NX)
      DATA IW1, IW2, TINY / 2*0, 1.E-35 /
C
      IF (DX .LE. 0.) THEN
c        CALL WARNR(IW2,NWRT,'DX cannot be < 0. in SMPSNA', 'DX', DX,
c     >               D0, D1, 0)
        SMPSNA = 0.
        RETURN
      ENDIF

      IF (NX .LE. 0 .OR. NX .GT. MAXX) THEN
c        CALL WARNI(IW1, NWRT, 'NX out of range in SMPSNA', 'NX', NX,
c     >               1, MAXX, 1)
        SIMP = 0.
      ELSEIF (NX .EQ. 1) THEN
        SIMP = 0.
      ELSEIF (NX .EQ. 2) THEN
        SIMP = (F(1) + F(2)) / 2.
        ERRD = (F(1) - F(2)) / 2.
      ELSE
        MS = MOD(NX, 2)

C For odd # of intervels, compute the diff between the Simpson 3/8-rule result
C for the last three bins and the regular Simpson rule result for the next to
C last two bins.  The problem is thereby reduced to a even-bin one in all cases.

        IF (MS .EQ. 0) THEN
          ADD = (9.*F(NX) + 19.*F(NX-1) - 5.*F(NX-2) + F(NX-3)) / 24.
          NZ = NX - 1
        ELSE
          ADD = 0.
          NZ = NX
        ENDIF

        IF (NZ .EQ. 3) THEN
          SIMP = (F(1) + 4.* F(2) + F(3)) / 3.
          TRPZ = (F(1) + 2.* F(2) + F(3)) / 2.
        ELSE
          SE = F(2)
          SO = 0
          NM1 = NZ - 1
          DO 60 I = 4, NM1, 2
            IM1 = I - 1
            SE = SE + F(I)
            SO = SO + F(IM1)
c	  print *, F(I), F(IM1)
   60     CONTINUE
          SIMP = (F(1) + 4.* SE + 2.* SO + F(NZ)) / 3.
          TRPZ = (F(1) + 2.* (SE + SO) + F(NZ)) / 2.
        ENDIF

        ERRD = TRPZ - SIMP
        SIMP = SIMP + ADD

      ENDIF
C
      SMPSNA = SIMP * DX

      IF (ABS(SIMP) .GT. TINY) THEN
        ERR = ERRD / SIMP
      ELSE
        ERR = 0.
      ENDIF
C
      RETURN
C                        ****************************
      END

