C**********************************************************************
C*     !!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!        *
C*                                                                    *
C* If you read this you might have already altered the program in an  *
C* unforseen way!!!! Delete it and download it again!!!!!!            *
C*                                                                    *
C* If you want to study the programm, please copy it to another file  *
C* name and do not use it for any scientific work!!!!!!!!             *
C*                                                                    *
C* Please contact us, if you want to seriously study the code in more * 
C* detail!                                                            *
C*                                                                    *
C**********************************************************************

      BLOCK DATA DATPDF 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C 
      CHARACTER*10 NAMPDF, NMHDRN, NAMQRK 
      CHARACTER*12 MRSFLN 
      LOGICAL LSTX 
C 
      PARAMETER (Z = 1D -10, ZZ = 1D -20) 
 
C              MxF is the maximum number of quark flavors 
C              MxAdF (= 1 gluon + 1 singlet quark + Additional fake flavors 
C                       to save space for storing info on evolution) 
      PARAMETER (MXX = 1050, MxQ = 25, MxF = 6, MxAdF = 6) 
C      PARAMETER (MxPN = MxF * 2 + MxAdF) 

      PARAMETER (MxPN = MxF * 2 + 2) 
      PARAMETER (MxQX= MxQ * MXX,   MxPQX = MxQX * MxPN) 
 
      PARAMETER (MxPDF = 20, MXHDRN = 6) 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
      Common / PdCntrl/ LPrt, LDbg 
 
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX, DEL,IV 
      COMMON / QARAY1 / QINI,QMAX, QV(0:MxQ),TV(0:MxQ), NT,JT,NG 
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx 
      COMMON / PEVLDT / UPD(MXPQX), KF, Nelmt 
      COMMON / PARPRZ / XMNPRZ(9), QMNPRZ(9), ALF(9), 
     >                  NfLPRZ(9), NfAL(9), MORD(9) 
      COMMON / INITIA / AN(-MxF : MxF), AI(-MxF : MxF), 
     >                  BI(-MxF : MxF), CI(-MxF : MxF) 
C 
      COMMON / NAMPDF / NAMPDF(MxPDF), NMHDRN(-1 : MxHDRN), 
     >                  NAMQRK(-MxF : MxF) 
      COMMON / MRSDAT / MRSFLN (3) 
C 
      Data LPrt, LDbg / 1, 0 / 
      DATA QINI, QMAX, XMIN, XCR / 1.9, 1.001D2, 0.999D-3, 0.3 / 
      DATA KF, IKNL, IPD0, IHDN / 10, 1, 1, 1 / 
      DATA NX, NT, JT, LSTX / 40, 6, 1, .FALSE. / 
      DATA (NfLPRZ(I), I=1,9) / 6,4,5,5,5,6,6,6,6/ 
      DATA (XMNPRZ(I), I=1,9) / 3*1.d-4,6*1.D-5/ 
      DATA (QMNPRZ(I), I=1,9) / 2.25, 2.0,3.2, 2*2.25, 4*2.0  / 
C Validity needs to be checked. 
      DATA (MORD(I), I=1,9) / 2*1, 7*2/ 
      DATA (ALF(I),  I=1,9) / 0.2, 0.2, 0.3, 0.19, 0.215, 
     $     0.22, 0.225,0.2,0.2/ 
      DATA (NfAL(I), I=1,9) / 3*4, 6*5 / 
C 
      DATA MRSFLN / 'prmz:hmrse', 'prmz:hmrsb', 'GRID3' / 
      DATA (NAMPDF(I),I=1,11) / 'Ehlq_1    ','Dk-Ow 1   ', 'DFLM_NLLA', 
     >'KMRS B0  ','MRS92_D0  ','MT90_S1  ','CTEQ1M   ', '          ', 
     >'         ', 'QCD_EVL_1','QCD_EVL_2' / 
C 
      DATA (NMHDRN(I), I=-1,5)/ 'AntiProton', 'Neutron', 'Proton', 
     > 'IsoScalar', 'Pi_Plus', 'Pi_Minus', 'K_Plus' / 
C 
      DATA (NAMQRK(I), I=-6,6)  / '-t_Quark', '-b_Quark', '-c_Quark', 
     >    '-s_Quark', '-d_Quark', '-u_Quark',    'Gluon',  'u_Quark', 
     >     'd_Quark',  's_Quark', 'c_Quark',   'b_Quark',  't_Quark' / 
C 
C                               The following data correspond to EHLQ set 1 
      DATA AI / 7 * Z,                      0.5,  0.4,  4 * Z / 
      DATA BI / 6 * Z,                3.5,  2 * 1.51,   4 * 1./ 
      DATA CI / 3 * 5, 3 * 8.54,      5.9,  3.5,  4.5,  4 * 5./ 
      DATA AN / 3 *ZZ, .081, 2 *.182, 2.62, 1.78, 0.67, 4 *ZZ / 
C 
C                        **************************** 
      END 
C
C 
      FUNCTION PDF (Iset, Ihadron, Iparton, X, Q, DEL, IR) 
 
 
C ==================================================================== 
C      FUNCTION PDF (Iset, Ihadron, Iparton, X, Q, IR) 
C     ------------------------------------------------- 
C   Front-end function for parton distributions of pbar, n,  p,  D,  C,  Fe 
C      all tranformed in terms of distribution of proton (PrtPdf) 
C   Also checks to see if pdf < 0; if so issue warning and set pdf=0. 
 
C     Revised 4/4/94 by HLL & WKT:  
C        PDF retains its name and argument list for compatibility with all 
C        existing programs which Call this function; PDF switches between 
C        target hadrons;  
C        PdfPrt is for proton target; it switches between different Iset's. 
C 
C     Ihadron = -1,  0,  1,  2,  4,  6 :  
C              pbar, n,  p,  D,  C,  Fe  
 
C               5 is "isoscalar-corrected iron" hence = D 
 
C         In all cases, adjust the Iparton label 
C         to convert to the corresponding proton distribution which is 
C         given in Function PdfPrt; 
 
C     --------------------------------------------------- 
 
C     If any of the arguments are unphysical, it issues a warning. 
C     Return error code IR: 
C          0  : O.K. 
C          1,2: Error from call to PdfPrt (proton parton distribution) 
C          3  : X < 0 or X > 1 : unphysical: warning + set PDF = 0; 
C          4  : Ihadron out of range; 
C          5  : PDF < 0.,        unphysical: warning + set PDF = 0; 
C          6  : Iparton out of range; warning + set PDF = 0; 
 
C     ------------------------------------- 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      Character Msg*80 
 
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1) 
 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx 
      
 
      DATA 
     1   HUGE, DUM,  D0m,   D1p  
     1 / 1D10, 0.D0, -1D-6, 1.000001D0 / 
     1   IW1,  IW2    / 2*0  / 
C                                 
      Ier = 0 
C                                                -------  Check x-range 
      IF (X.LE.D0m .or. X.GT.D1p) Then 
        Ier = 3 
        Call WARNR(IW1,NWRT,'X out of range in PDF.', 'X', X, 
     >               D0, D1, 1) 
        TEM = DUM 
      EndIf 
C                                              Check bounds on Ihadron 
      If (Ihadron.gt.6 .or. Ihadron.lt.-1 .or. Ihadron.eq.3) then 
         Call WARNI(IW, NWRT, 
     >    'Only Ihardon=-1,0,1,2,4,5,6 (pbar,n,p,D,*,) are active', 
     >    'Ihadron', Ihadron, -1,6,1) 
         Ir = 4 
       PDF=0.D0 
       Return 
      Endif 
C                                              Check bounds on Iparton 
      Jp=Abs(Iparton) 
      Neff = NFL(Q) 
c nfl(q) returns the number of `light' flavors at scale Q - effective 
      If ( Jp .gt. NEFF) then 
C                                   if Jp > Neff, then set PDF=0 and return 
         Call WARNI(IW, NWRT, 
     >    'Iparton out of range', 
     >    'Iparton', Iparton, -Neff, Neff,1) 
         Ier = 6 
           PDF = 0D0 
         Return 
      Endif 
 
C                   --- Conversion of  Ihadron  to proton distributions, 
C                                     if necessary ----  
      If (Jp.eq.1 .or. Jp.eq.2) Then 
C                                              For u and d 
C                                   Use Isospin symmetry n<->p  == u<->d 
         Ipartner=3-Jp 
         If (Iparton.lt.0) Ipartner=-Ipartner 
 
         If (Ihadron.eq.1) then 
C p: 
            Tem= PdfPrt(Iset, Iparton, X, Q, Ir) 
         Elseif (Ihadron.eq.-1) then 
C pbar: 
            Tem= PdfPrt(Iset, -Iparton, X, Q, Ir) 
         Elseif (Ihadron.eq.0) then 
C n: 
            Tem= PdfPrt(Iset, Ipartner, X, Q, Ir) 
         Elseif (Ihadron.eq.2 .or.Ihadron.eq.4 .or.Ihadron.eq.5) then 
C isoscalar: D, C, ... 
            Tem=( PdfPrt(Iset, Iparton, X, Q, Ir) 
     >           +PdfPrt(Iset, Ipartner, X, Q, Ir) )/2.D0 
         Elseif (Ihadron.eq.6) then 
C Fe: 
            Tem=( 26.D0* PdfPrt(Iset, Iparton, X, Q, Ir) 
     >           +30.D0* PdfPrt(Iset, Ipartner, X, Q, Ir) )/56.D0 
         Endif 
      Else 
C                                              For s,c,b,t 
         Tem= PdfPrt(Iset, Iparton, X, Q, Ir) 
      Endif 


C                                      --- Make sure PDF >= 0 -------- 
C                                       (unless Iknl<0 - polarized pdf) 
c      IF (TEM .LT. D0 .and. Iknl .ge. 1 .and. X.GE.DEL) Then 
c        IF (TEM .LT. D0m .AND. X .LE. 0.9D0) Then 
c        Call WARNR(IW2,NWRT,'PDF < 0; Set --> 0', 'PDF',TEM,D0,HUGE,1) 
c        WRITE (NWRT, '(A, 2I5, 2(1PE15.3))') 
c     >      ' ISET, Iparton, X, Q = ', ISET, Iparton, X, Q 
c        Ier = 5 
c        EndIf 
c        TEM = D0 
c      EndIf 
C                        -------- Return function value and error code 
      PDF = TEM 
      IR  = Ier 

      return
C
      end
C
      	      FUNCTION PdfPrt (IPDF, LPRTN, XD, QD, IR) 
C 
C 
C    This routine gives the parton ( IPARTN ) distribution function inside 
C    the proton in a chosen evolved or parametrized form ( Ipdf ) 
 
C  2. For Ipdf = 10 or 11:  steer to results from QCD evolution program; 
C                1  -   9:           CTEQ parametrizations 
C                12 -  20:           CTEQ tables or other forms   
 
C     > 30 : other parametrizations: see codes below for details.              
 
C    It gives the probability distribution, not the momemtum-weighted one. 
C 
C                       Return Code: IER = 0:   No error 
c                                        < 10:   PDF fundtion error 
C                                        = 10:   Ipdf out of range 
c                                        > 10:   multi-error 
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1) 
C 
      PARAMETER (MXPDF = 20, MXHDRN = 6, MXF = 6, MXPN = MXF*2+2) 
C 
      CHARACTER MSG*75 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
C 
      DATA IW1 / 0 / 
 
      IER = 0 
      Irr = 0 
C 
      Iprtn = LPRTN 
      x = XD 
      Q = QD 
C                              ====  Begin of overall Ipdf IfBlock ===== 
      If     (Ipdf .eq.  10) Then 
C                                                             10 evolve 
                            Tmp = ParDis (Iprtn, X, Q) 

      Else 
                            Ier  = 10 
      MSG=  
     >'Ipdf chosen is currently inactive. PdfPrt set equal to zero.' 
      CALL WARNI (IW1, NWRT, MSG, 'Ipdf', Ipdf, 1, 10, 0) 
                            Tmp = 0. 
      EndIf 
 
      PdfPrt = Tmp 
      Ir = Ier + Irr 
 
   10 RETURN 
C                        **************************** 
      END 
  

      FUNCTION PARDIS (IPRTN, XX, QQ) 
C 
C       Given the parton distribution function in the array U in 
C       COMMON / PEVLDT / , this routine fetches u(fl, x, q) at any value of 
C       x and q using Mth-order polinomial (II=0) or rational fraction (II=1)  
C       interpolation for x. It always uses quadratic polinomial interpolation 
C       in ln ln (Q/lambda). 
 
C       The calling program must ensure that 0 =< x =< 1 ; 
C       If 0 =< x < Xmin, extrapolation is used and a warning is given 
 
C       The calling program must ensure that Alambda < Q ; 
C       If (Alambda < Q < Qini .or. Qmax < Q), 
C          extrapolation is used and a warning is given 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      Character Msg*80 
      LOGICAL LSTX 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (MXPN = MXF * 2 + 2) 
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN) 
      PARAMETER (Smll = 1D-9) 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
 
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX, DEL,IV 
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX) 
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG 
      COMMON / QARAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF) 
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx 
      COMMON / PEVLDT / UPD(MXPQX), KF, Nelmt 
C 
      Dimension Fq(5), Df(5) 

C                           M determines the order of the polynomial 
C            II switches between the polint/ratint interpolation routines 
C                                     They are fixed for this version 
      Save  
      Data M , II / 2, 0 / 
      Data Iwrn1, Iwarn2, Iwarn3 / 3*0 / 
 
      X = XX 
      Q = QQ 

      Md = M / 2 
      Amd= Md 

C 
      IF (Q .LT. QINI-Smll) THEN 
         Msg = 'Q less than QINI in PARDIS call; Q SET = QINI.' 
      CALL WARNR (IWRN1, NWRT, Msg, 'Q', Q, QINI, QMAX, 1) 
         Q = QINI 
      ElseIF (Q .GT. QMAX) THEN 
         Msg = 'Q greater than QMAX in PARDIS call; '  
         Msg = Msg // 'Extrapolation will be used.' 
         CALL WARNR(IWRN2, NWRT, Msg, 'Q', Q, QINI, QMAX, 1) 
      EndIf 
C                           Find lower end of interval containing X 
      JL = -1 
      JU = Nx+1 
 11   If (JU-JL .GT. 1) Then 
         JM = (JU+JL) / 2 
         If (X .GT. XV(JM)) Then 
            JL = JM 
         Else 
            JU = JM 
         Endif 
         Goto 11 
      Endif 
 
      Jx = JL - (M-1)/2 
      If     (Jx .LT. 0) Then 
         Jx = 0 
      Elseif (Jx .GT. Nx-M) Then 
         Jx = Nx - M 
      Endif 
C                                    Find the interval where Q lies 
      JL = -1 
      JU = NT+1 
 12   If (JU-JL .GT. 1) Then 
         JM = (JU+JL) / 2 
         If (Q .GT. QV(JM)) Then 
            JL = JM 
         Else 
            JU = JM 
         Endif 
         Goto 12 
       Endif 
 
      Jq = JL - (M-1)/2 
      If     (Jq .LT. 0) Then 
         Jq = 0 
      Elseif (Jq .GT. Nt-M) Then 
         Jq = Nt - M 
      Endif 
 
      SS  = LOG (Q/AL) 
      TT  = LOG (SS) 
C                             Find the off-set in the linear array Upd 
      JFL = IPRTN + NfMx 
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx + M/2
C 
C                    Nt  -| .......................... 
C                         | .......................... 
C                   Jq+M -| .....o ......o............  Iq=M+1 
C                         | .........X................ 
C                    Jq  -| .(J0)o ......o ...........  Iq=1 
C                         | ...... ................... 
C                     0  --------|-------|-----------| 
C                         0     Jx     Jx+M          Nx 
 
 
      Do 21 Iq = 1, M+1 
      J1 = J0 + (Nx+1)*(Iq-1) + 1 

      Fq(Iq) = Upd(J1)

 21   Continue 

C                                                     Interpolate in LnLnQ 
      Call Polint (TV(Jq), Fq(1), M+1, TT, Ftmp, Ddf) 

      PARDIS = Ftmp 
C 
      RETURN 
C                        **************************** 
      END 


       SUBROUTINE EVOLVE (NPO,FINI, IRET) 

C ==================================================================== 
C      SUBROUTINE EVOLVE (NPO,FINI, IRET) 
C 
C               Input argument: FINI is a function 
C 
C                               FINI (LPARTN, X) 
C 
C                     where LPARTN = -6, ... 6 labels the parton flavor: 
C                     t-bar(-6), b-bar(-5), ... gluon(0), u(1), ... t(6) res. 
C 
C               Output parameter: 
C 
C                       Iret = 0  :  normal execution 
C 
C     NPO : 1 = evolution of H and \tilde H
C           2 = evolution of E and \tilde E
C                       ---------------------------- 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      LOGICAL LSTX 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (MXPN = MXF * 2 + 2) 
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN) 
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2) 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX, DEL,IV 
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG 
      COMMON / QARAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF) 
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx 
      COMMON / PEVLDT / UPD(MXPQX), Kf, Nelmt 
C 
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX) 
      COMMON / VARBAB / GB(NDG, NDH, MXX), H(NDH, MXX, M1:M2) 
C 
      DIMENSION QRKP(MXF) 
      DIMENSION JI(-MXF : MXF+1)
      
C 
      EXTERNAL NSRHSP, NSRHSM, FINI 
 
      DATA D0 / 0.0 / 
C 
   11 IRET = 0 
C     Set up number of "valence quarks" 

      IF (IHDN .LE. 4) THEN 
        MXVAL = 2 
      ElseIF (IHDN .LE. 6) THEN 
        MXVAL = 3 
      EndIf 
C                                    Set up X-mesh points and common parameters 
      IF (.NOT. LSTX) CALL XARRAY 
      DLX = 1D0 / NX 
      J10 = NX / 10 

C                                                       Set up Q-related arrays 
      CALL PARPDF (2, 'ALAM', AL, IR) 
C                         Nini and NfMx determined by Qini and Qmax in Qarray 
      CALL QARRAY (NINI) 

C                            Flavor # NFSN is for the singlet quark combination 
      NFSN = NFMX + 1 
C                  Kf is the range of the flavor index for the Upd data array 
C                 = quark flavors + 1 for singlet combination and 1 for gluon 
      KF = 2 * NFMX + 2 
C                                 Total number data points to be stored in Upd 
      Nelmt = Kf * (Nt+1) * (Nx+1) 
C                                                       ---------------------- 
C                                                        Various initiations: 
C 
C      Offset, or Starting address -1, of Pdf(Iflv) in Upd for the first stage. 
C                                                    Will be updated each turn 
      DO 101 IFLV = -NFMX, NFMX+1 
        JFL = NFMX + IFLV 
        JI(IFLV) = JFL * (NT+1) * (NX+1) 
  101 CONTINUE 
C                               Define initial distributions for the evolution. 
C                                        Input functions are defined on (1:Nx). 
C 
C                             Input gluon and quark distributions are inserted 
C                          in the Upd array directly;   Output of each stage of 
C                              evolution serve as the input for the next stage. 
C 
    3 DO 31 IZ = 1, NX 

C                              Gluon:  (Position Ji(0)+1 is reserved for x = 0) 
C 
        UPD(JI(0)+IZ+1) = FINI (NPO,0, XV(IZ)) 

C                                                                 Singlet quark 
C                        Input singlet quark distribution must be initiated for 
C                                                        each stage separately. 
C 
        UPD(JI(NFSN)+IZ+1) = 0 
C                                   If no quark, bypass filling of quark arrays 
        IF (NFMX .EQ. 0) GOTO 31 
C 
        DO 331 IFLV = 1, NINI 
          A = FINI (NPO, IFLV, XV(IZ)) 
          B = FINI (NPO,-IFLV, XV(IZ)) 
          IF (XV(IZ).GT.DEL) THEN
          QRKP (IFLV) = A + B 
          ELSE
          QRKP(IFLV) = A
          ENDIF
C                                               Acculumate singlet distribution 
          UPD(JI(NFSN)+IZ+1) = UPD(JI(NFSN)+IZ+1) + QRKP (IFLV) 

C                                       Initialize the "minus" non-singlet com- 
C                                        bination (Q - Q-bar) in area for Q-bar 
          IF (XV(IZ).GT.DEL) THEN
          UPD(JI(-IFLV)+IZ+1) = A - B
          ELSE
          UPD(JI(-IFLV)+IZ+1) = B
          ENDIF
  331   CONTINUE 
C                            The "plus" non-singlet combination is initialized 
C                           in the array area for the Quark of the same flavor 
        DO 332 IFLV = 1, NINI 
           UPD(JI(IFLV)+IZ+1) = QRKP(IFLV) - UPD(JI(NFSN)+IZ+1)/NINI 
  332   CONTINUE 
C

   31 CONTINUE

C                                                       ----------------------- 
C                                                   Start of the Q2- Evolution: 
C 
C                                 Outer loop is by the effective number of flvr 
      DO 21 NEFF = NINI, NFMX 
C                                    Set up 1st and 2nd order kernel functions

          IF (IKNL .EQ. 2 .OR. IKNL.EQ.-2) CALL STUPKL (NEFF) 
C                                                           Singlet Calculation 
          ICNT = NEFF - NINI + 1

C                                                 Skip if new quark mass = old 
          IF (NTN(ICNT) .EQ. 0) GOTO 21 
C                                       Otherwise, recall iteration parameters 
          NITR = NTN (ICNT) 
          DT   = DTN (ICNT) 
          TIN  = TLN (ICNT)

C                                                            Perform evolution 
          CALL SNEVL (NPO,XV,IV,DEL,IKNL, NX, NITR, JT, DT, TIN, NEFF, 
     >    UPD(JI(NFSN)+2), UPD(JI(0)+2), UPD(JI(NFSN)+1), UPD(JI(0)+1)) 
C 
C                                        Non-singlet sector, one flavor a time. 
C 
C                                        Skip this section if quark flavor = 0.

          IF (NEFF .EQ. 0) GOTO 88 

    5     DO 333 IFLV = 1, NEFF 
C                                    First evolve the (Q+Q-bar) "PLUS NS" part 
           CALL NSEVL (NPO,0,NSRHSP,XV,IV,DEL,IKNL, NX, NITR, JT, DT, 
     >                 TIN, NEFF, UPD(JI( IFLV)+2), UPD(JI( IFLV)+1)) 
C 
C                                  The (Q-Q-bar) "MINUS NS" evolution is needed 
C                                           only for those flavors with valence 
           IF (IFLV .LE. MXVAL) 
     >     CALL NSEVL (NPO,1,NSRHSM,XV,IV,DEL,IKNL, NX, NITR, JT, DT, 
     >                 TIN, NEFF, UPD(JI(-IFLV)+2), UPD(JI(-IFLV)+1)) 
C 
C                            To obtain  the real  quark distribution functions, 
C                      combine the singlet piece to the two non-singlet pieces. 
C                            Enforce positivity conditions also at this stage. 
C 
           DO 55 IS = 0, NITR 
           DO 56 IX = 0, NX 

             TP = UPD (IS*(NX+1) + IX + 1 + JI( IFLV)) 
             TS = UPD (IS*(NX+1) + IX + 1 + JI( NFSN)) / NEFF 
             TP = TP + TS 
             IF (IKNL .GT. 0 .AND. XV(IX).GT.DEL) TP = TP 
 
             IF (IFLV .LE. MXVAL) THEN 
                TM = UPD (IS*(NX+1) + IX + 1 + JI(-IFLV)) 
                IF (IKNL .GT. 0 .AND. XV(IX).GT.DEL) THEN 
                  TM = TM 
c                  TP = MAX (TP, TM) 
                  TP = TP
                EndIf 
             Else 
                TM = 0. 
             EndIf 
 
             IF (XV(IX).GT.DEL) THEN
             UPD (JI( IFLV) + IS*(NX+1) + IX + 1) = (TP + TM)/2. 
             UPD (JI(-IFLV) + IS*(NX+1) + IX + 1) = (TP - TM)/2. 
             ELSE
             UPD (JI( IFLV) + IS*(NX+1) + IX + 1) = TP
             UPD (JI(-IFLV) + IS*(NX+1) + IX + 1) = TM
             ENDIF
 
   56      CONTINUE 

           IF (IS.EQ.0) THEN 

              GOTO 55

           ELSE

              ENDIF

   55      CONTINUE 
333      CONTINUE 

C                   Heavy quarks above current threshold have zero distribution 
C 
        DO 334 IFLV = NEFF + 1, NFMX 
          DO 57 IS = 0, NITR 
          DO 58 IX = 0, NX 
            UPD(JI( IFLV) + IS*(NX+1) + IX + 1) = 0 
            UPD(JI(-IFLV) + IS*(NX+1) + IX + 1) = 0 
   58     CONTINUE 
   57     CONTINUE 
  334   CONTINUE 
   88   CONTINUE 
C                                               ------------------------------ 
C                                      Define initial parameters for next stage 
        IF (NFMX .EQ. NEFF) GOTO 21 
C                                                                   New Offsets 
        DO 335 IFLV = -NFMX, NFMX+1 
           JI(IFLV) = JI(IFLV) + NITR * (NX+1) 
  335   CONTINUE 
C                                                            New distributions: 
C                                           gluon input functions are in place; 
C                       only non-singlet and singlet quark needs re-initiation 
C 
C                             Calculate initial heavy quark distribution due to 
C                           change of renormalization scheme across threshold: 
        CALL HQRK (NX, TT, NEFF+1, UPD(JI(0)+2), UPD(JI(NEFF+1)+2)) 
C 
        DO 32 IZ = 1, NX 

         QRKP (NEFF+1) = 2. * UPD(JI( NEFF+1) + IZ + 1) 
C                                                             New Singlet piece 
         UPD (JI(NFSN)+IZ+1) = UPD (JI(NFSN)+IZ+1)  + QRKP (NEFF+1) 
         VS00 =  UPD (JI(NFSN)+IZ+1) / (NEFF+1) 

C 
C                                        "plus" non-singlet for the new flavor 
         UPD(JI( NEFF+1) + IZ + 1) = QRKP(NEFF+1) - VS00 
C 
C                 Calculate the non-singulet parts of the other quark distr. 
C 
C                               Change from the output of last stage of calcu- 
C                               lation due to two sources of change in Vs/Neff: 
C                               change of Neff and addition of new quark distr. 
         DO 321 IFL = 1, NEFF 
           A = UPD(JI( IFL)+IZ+1) 
           B = UPD(JI(-IFL)+IZ+1) 
           
           IF (XV(IZ).GT.DEL) THEN
             QRKP(IFL) = A + B 
           ELSE
             QRKP(IFL) = A
           ENDIF
C                                             "plus" non-singlet for flavor IFL 
           UPD(JI( IFL)+IZ+1) = QRKP(IFL) - VS00 
C                                 "minus" non-singlet for flavors with valence 
           IF (XV(IZ).GT.DEL) THEN
              IF (IFL .LE. MXVAL)  UPD(JI(-IFL)+IZ+1) = A - B 
           ELSE
              IF (IFL .LE. MXVAL)  UPD(JI(-IFL)+IZ+1) = B 
           ENDIF
  321    CONTINUE 
C 
   32   CONTINUE 
C                               Return of Q-2 evolution loop to the next stage 
   21 CONTINUE 
C                                            Conclusion of the full calculation 
      Return 
C                   ********************** 
      End 

      SUBROUTINE NSEVL (NPO,NSET,RHS,XT,IV,DEL,
     > IKNL,NX,NT,JT,DT,TIN,NEFF,U0,UN)
C
C                               IKNL determines to which order in Alpha is the
C                                               calculation to be carried out.
C
C       Given the non-singlet parton distribution function U0 at some initial
C       QIN (Tt= 0) in the x interval (0, 1) covered by the array Iz = 1, NX,
C       this routine calculates the evoluted function U at Nt values of Tt at
C       intervals of Dt by numerically integrating the A-P equation using the
C       non-singlet kernel.
C
C       Un(Ix, Tt) = Y(x,es) at the sites Ix=0,..,Nx  (x = Ix * Dz);
C       Un(0, Tt) is obtained by quadratic extrapolation from Ix = 1, 2, 3
C                 for each Tt rather then by evolution because of possible
C                                               singular behavior at x = 0.
C
C       Data is stored at Tt = Is * Dt, Is = 0, ... , Nt.
C                       The function at Is = 0 is the input distribution.
C
C       Numerical iteration is performed with finer grain Ddt = Dt/Jt.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      Common / PdCntrl/ LPrt, LDbg
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      
C
      DIMENSION  U0(NX), UN(0:NX, 0:NT)
      DIMENSION  Y0(MXX), Y1(MXX), YP(MXX), F0(MXX), F1(MXX), FP(MXX)
      DIMENSION XT(0:MXX)
      external rhs

C
         If (Ldbg .Eq. 1) Then
            Write (Nwrt, '(A)') ' Non-Singlet:'
            N5 = Nx / 5 + 1
         Endif

      DDT = DT / JT
C
      IF (NX .GT. MXX) THEN
      WRITE (NOUT,*) 'Nx =', NX, ' greater than Max # of pts in NSEVL.'
      STOP 'Program stopped in NSEVL'
      EndIf
C
C                                       Compute an effective first order lamda
C                                       to be used in checking of moment evl.
      TMD = TIN + DT * NT / 2.
      AMU = EXP(TMD)
      TEM = 6./ (33.- 2.* NEFF) / ALPI(AMU)
      TLAM = TMD - TEM
C
C    Fill first rows of output array, for initial value of Q
C
      DO 9 IX = 1, NX
      UN(IX, 0)  = U0(IX)
    9 CONTINUE
      UN(0, 0) = 3D0*U0(1) - 3D0*U0(2) - U0(1)
C
C    Initiation
C
      TT = TIN
      DO 10 IZ = 1, NX
      Y0(IZ)   = U0(IZ)
   10 CONTINUE
C
C    loop in the Tt variable
C
      DO 20 IS = 1, NT
C
C    fine- grained iteration
C
         IPO = JT

         DO 202 JS = 1, JT
C
C    Irnd is the counter of the Q-iteration
C
            IRND = (IS-1) * JT + JS
C
C   Use Runge-Katta the first round
C
            IF (IRND .EQ. 1) THEN
C
                CALL RHS (NPO,TT, Neff, Y0, F0)
C smoothing

                CALL SMOOTH (NSET,IKNL,MXX,DEL,IV,XT,F0)

                DO 250 IZ = 1, NX
                   Y0(IZ) = Y0(IZ) + DDT * F0(IZ)
  250           CONTINUE
C     smoothing

                CALL SMOOTH (NSET,IKNL,MXX,DEL,IV,XT,Y0)
C
                TT = TT + DDT
C
                CALL RHS (NPO,TT, NEFF, Y0, F1)
C     smoothing

                CALL SMOOTH (NSET,IKNL,MXX,DEL,IV,XT,F1)

                DO 251 IZ = 1, NX
                   Y1(IZ) = U0(IZ) + DDT * (F0(IZ) + F1(IZ)) / 2D0
  251           CONTINUE
c     smoothing

                CALL SMOOTH (NSET,IKNL,MXX,DEL,IV,XT,Y1)

C
C    What follows is a combination of the 2-step method
C       plus the Adams Predictor-Corrector Algorithm
C
            Else
C
                CALL RHS (NPO,TT, NEFF, Y1, F1)

c     smoothing

                CALL SMOOTH (NSET,IKNL,MXX,DEL,IV,XT,F1)
C
C    Predictor
C
                DO 252 IZ = 1, NX
                   YP(IZ) = Y1(IZ) + DDT * (3D0 * F1(IZ) - F0(IZ)) / 2D0
  252           CONTINUE
c     smoothing

                CALL SMOOTH (NSET,IKNL,MXX,DEL,IV,XT,YP)

C
C  Increment of Tt at this place is part of the formalism
C
                TT = TT + DDT
C
                CALL RHS (NPO,TT, NEFF, YP, FP)
C     smoothing

                CALL SMOOTH (NSET,IKNL,MXX,DEL,IV,XT,FP)
C
C   Corrector
C
                DO 253 IZ = 1, NX
                   Y1(IZ) = Y1(IZ) + DDT * (FP(IZ) + F1(IZ)) / 2D0
                   F0(IZ) = F1(IZ)
  253           CONTINUE
c     smoothing
               
                CALL SMOOTH (NSET,IKNL,MXX,DEL,IV,XT,Y1)
 
            EndIf
C
  202    CONTINUE

C
C    Fill output array
C
         DO 260 IZ = 1, NX
            UN (IZ, IS) = Y1(IZ)
  260    CONTINUE
C
C               The value of the function at x=0 is obtained by extrapolation
C
         UN(0, IS) = 3D0*Y1(1) - 3D0*Y1(2) + Y1(3)
C
C    Print out for Debugging
C
      If (LDbg .Eq. 1) Then
         Write (Nwrt, '(A, 5(1pE12.3))') '   :', (Un(Iz,Is), Iz=1,Nx,N5)
      Endif
C
      JT = IPO

   20 CONTINUE
C
      RETURN
C                        ****************************
      END
      

 
      SUBROUTINE NSRHSM (NPO,TT, NEFF, FI, FO) 
C 
C       Subroutine to compute the Right-Side of the Altarelli-Parisi Equation 
C                This copy applies to the "NS-minus" piece -- (Qrk - Qrk-bar) 
C 
C       See comments in NSRHSP for details 
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      LOGICAL LSTX 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2) 
      PARAMETER (MXQX = MXX*MXQ) 
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2,  
     >Z3 = 2.404113806319188) 
      PARAMETER (D0 = 0.0, D1 = 1.0)	 
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX) 
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX, DEL, IV 
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX) 
      COMMON / KRNL01 / AGG2 (0:MXX), ANSP (0:MXX), ANSM (0:MXX) 
      COMMON / KRN1ST / FF1(0:MXX,0:MXX), GG1(0:MXX,0:MXX),    
     > FG2(0:MXX,0:MXX),GF2(0:MXX,0:MXX),GG2(0:MXX,0:MXX), 
     > PNSP(0:MXX,0:MXX), PNSM(0:MXX,0:MXX), SFF2(0:MXX,0:MXX), 
     > SFF2BLK(0:MXX,0:MXX), FG2BLK(0:MXX,0:MXX), GF2BLK(0:MXX,0:MXX),  
     > GG2BLK(0:MXX,0:MXX), PNSPBLK(0:MXX,0:MXX), PNSMBLK(0:MXX,0:MXX), 
     > GG2BLC1(0:MXX,0:MXX), PNSPBLC1(0:MXX,0:MXX),FF1BLX(0:MXX,0:MXX),  
     > GG2BLCX(0:MXX,0:MXX), PNSPBLCX(0:MXX,0:MXX),GG1BLX(0:MXX,0:MXX), 
     > PNSMBLC1(0:MXX,0:MXX), PNSMBLCX(0:MXX,0:MXX),FF1BL1(0:MXX,0:MXX), 
     > GG1BL1(0:MXX,0:MXX),SFF2BLKD(0:MXX,0:MXX), FG2BLKD(0:MXX,0:MXX), 
     > GF2BLKD(0:MXX,0:MXX), GG2BLKD(0:MXX,0:MXX),PNSPBLKD(0:MXX,0:MXX), 
     > PNSMBLKD(0:MXX,0:MXX),FG2BLK1(0:MXX,0:MXX),GF2BLK1(0:MXX,0:MXX),
     > SFF2BLK1(0:MXX,0:MXX),PNSPC(0:MXX,0:MXX),PNSMC(0:MXX,0:MXX),
     > GG2C(0:MXX,0:MXX), PNSPDGBL(0:MXX,0:MXX), PNSMDGBL(0:MXX,0:MXX),
     > GG2DGBL(0:MXX,0:MXX)
      COMMON / REGFUNC / GB1(MXX), GBX(MXX), GB21(MXX), GB510(MXX),  
     > GB2X(MXX), GB31(MXX), GB41(MXX), GB5X(MXX), GB51(MXX),  
     > GB3X(MXX), GB4X(MXX), GB210(MXX), GB10(MXX), GB310(MXX), 
     > G1(MXX), G2(MXX), G3(MXX), G4(MXX), G5(MXX), GB410(MXX)

      COMMON / DIFFD / DIFFDEL
      
c      COMMON / SECDERIV / DGBX(3,MXX), DGB2X(3,MXX), DGB1(3,MXX), 
c     > DGB21(3,MXX),
c     > DGB01(3,MXX), DGB41(3,MXX), DGB10(3,MXX), DGB210(3,MXX), 
c     > DG1(3,MXX), DG2(3,MXX)
 
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx 
 
      DIMENSION  FI(NX), FO(NX), FN(MXX), FFN(3*MXX),FQ(MXX)
      DIMENSION W0(MXX), W1(MXX), WH(MXX), W2(MXX), W01(MXX) 
      DIMENSION W3(MXX), W5(MXX), W22(MXX), FFI(3*MXX), WH1(MXX) 
      DIMENSION X2(10), X1(10), FX(10), FX1(10), AF(MXX),  AF1(MXX), 
     > AF2(MXX), TQQ1(MXX), TQQX(MXX),TQQ11(MXX),TQQD(MXX)
      DIMENSION TEMPI1(0:MXX),TEMPI2(0:MXX)
       
	EXTERNAL SMPSNA,SMPNOL,SMPNOR
 
	SAVE 
	DATA AERR, RERR / 0.00001, 0.0002/ 
C 
      S = EXP(TT) 
      Q = AL * EXP (S) 
      CPL = ALPI(Q) 
      CPL2= CPL ** 2 / 2. * S 
      CPL = CPL * S 
C 
C	Set up Quark pdf for integration in BL region! 
C 

	DZ = 1./(NX-1) 

	FN(NX) = 0.0
	FQ(NX) = 0.0


	DO 229 IY = 1, NX-1
           FQ(IY) = FI(IY)*XA(IY,1)
           TMX = ABS(XV(IY)-DEL/2.)
        IF (TMX.GT.1E-10) then
	FN(IY) = FI(IY)        
        ELSE
           IF (IKNL.LT.0) then
              FN(IY) = 0.0
           ELSE
             FN(IY) = FI(IY) 
           ENDIF
        ENDIF

 229    CONTINUE 

C 
C	LO convolution of kernel with pdf! First unpol. then pol. case! 
C 
      CALL NEWARRAY (NX, DEL, FQ, FFI) 
      CALL NEWARRAY (NX, DEL, FN, FFN) 
      CALL INTEGR (NX, 0, FQ, W0, IR1) 
      CALL INTEGR (NX, 0, FN, W01, IR12) 
      CALL NINTEGR (NX, DEL, 1, FFI, FQ, W1, IR11) 
      CALL NINTEGR (NX, DEL, 15, FFN, FN, W2, IR01) 
      CALL NINTEGR (NX, DEL, 8, FFI, FQ, W22, IR02) 
      CALL NINTEGR (NX, DEL, 7, FFI, FQ, W3, IR22) 
      CALL NINTEGR (NX, DEL, 12, FFN, FN, W5, IR3) 
      CALL NINTEGR (NX, DEL, 16, FFN, FN, AF, IR4)
      CALL NINTEGR (NX, DEL, 20, FFN, FN, AF2, IR4)
      CALL HINTEGN (2,NX,IV,DEL,FFI,FQ, WH)
      CALL HINTEGN (1,NX,IV,DEL,FFN,FN, WH1) 	
C 
     
	DIM1 = 1D0/DEL 
 
      DO 230 IZ = 2, NX 
C
C     First treatment of ERBL region
C        
      IF (XA(IZ,1) .LE. DEL) THEN 

	AA1 = 1D0 - XA(IZ,1)*DIM1 
C
C     Special treatment of the point y=del and y= 0
C
        IF (IZ.EQ.IV) THEN
C
C     For H and \tilde H use normal extrapolation and then symmetry properties
C     for E and \tilde E set those points to zero!
C
        IF (NPO.EQ.1) THEN

          IF (IKNL.GT.0) THEN

      DO 11 IL = 1,7
         X2(IL) = 1D0*IL
         FX(IL) = FO(IL+1)
         
 11      CONTINUE

         XX = 0D0

         CALL RATINT(X2,FX,7,XX,TEM,ERR)

         FO(1) = TEM 

      FO(IZ) = FO(1)

        ELSE

      DO 1 IL = 1,7
c         X2(IL) = XV(IL+1)
         X2(IL) = 1D0*IL
         FX(IL) = FO(IL+1)
         
 1      CONTINUE

c         XX = XV(1)
        XX = 0D0

         CALL RATINT(X2,FX,7,XX,TEM,ERR)

         FO(1) = TEM

      FO(IZ) = -FO(1)

        ENDIF 
 
        ELSE
C
C     now E and \tilde E
C
        FO(1) = 0D0
        FO(IV) = 0D0

        ENDIF

           ELSE
C
C     separate treatment for E and \tilde E once more, but first do H and \tilde H
C

       IF (NPO.EQ.1) THEN
C
C     treatment of H and \tilde H
C

        IF (IKNL.GT.0) THEN

      FO(IZ) = 2.*FN(IZ) + 4./3.*(WH1(IZ) - AF(IZ)* 
     > AA1 + AA1*W5(IZ) + W2(IZ) + (AF2(IZ) - (1D0-AA1)*AF(IV)))

        ELSE

      FO(IZ) = 2.*FN(IZ) + 4./3.*(WH1(IZ) - AF(IZ)* 
     > AA1 + AA1*W5(IZ) + W2(IZ) - (AF2(IZ) - (1D0-AA1)*AF(IV)))

      TMX = ABS(XV(IZ)-DEL/2.)
        IF (TMX.LT.1E-10) FO(IZ) = 0.0

        ENDIF

        ELSE
C
C now E and \tilde E
C        

       IF (IKNL.GT.0) THEN

      FO(IZ) = 2.*FN(IZ) + 4./3.*(WH1(IZ) - AF(IZ)* 
     > AA1 + AA1*W5(IZ) + W2(IZ))

        ELSE

      FO(IZ) = 2.*FN(IZ) + 4./3.*(WH1(IZ) - AF(IZ)* 
     > AA1 + AA1*W5(IZ) + W2(IZ))

      TMX = ABS(XV(IZ)-DEL/2.)
        IF (TMX.LT.1E-10) FO(IZ) = 0.0

         ENDIF

        ENDIF

        FO(IZ) = CPL*FO(IZ)

       ENDIF

      ELSE

C
C     set pdf = 0 for DGLAP region if E or \tilde E
C
      IF (NPO.EQ.1) THEN
      
      FO(IZ) = 2.* FQ(IZ) + 4./3.* (2D0*WH(IZ) - W0(IZ)  
     >         - W1(IZ) - W3(IZ) + W22(IZ))

      FO(IZ) = CPL * FO(IZ)*XA(IZ,-1) 

         ELSE

         FO(IZ) = 0D0

         ENDIF

      ENDIF
     
 230  CONTINUE 

C
C     Now NLO treatment
C

      IF (IKNL .EQ. 2 .OR. IKNL.EQ.-2) THEN

           IF (IKNL.EQ.2) THEN

         ET = 1D0
	ET1 = 1D0
        ET4 = 1D0
        ET2 = 1D0  

        ELSE

         ET = 1D0
	ET1 = 1D0
        ET4 = 1D0
        ET2 = 1D0  

        ENDIF   

      DO 21 IX = 2, NX-1

         IF (IX.EQ.IV) GOTO 21

        X = XV(IX)

        NR = IX - IV
C
C       Evaluate the integrand for the I-th integral
C
C	First the BL integrands and then the DGLAP integrands.
C
C

        GB10(NX-IV+1) = 0.0

	IF (IX.LE.IV) THEN

	IQ = 1

	IQ1 = NX

	ELSE

	IQ = IX

	IQ1 = NX

	ENDIF
	
        DO 31 KZ = IQ, IQ1

	IF (IX.EQ.KZ) THEN

           IF (IX.LE.IV) THEN
C
C     Values of the integrand at y=x. Dummy values only for book keeping purposes.
C
        GBX(1) = 0.0

        GB1(1) = 0.0

        ENDIF
           
	GOTO 31

	ELSE

	IF (IX.GT.KZ .AND. IX.LE.IV) THEN

	GBX(KZ) = DXTZ(KZ)*PNSMBLK(IX,KZ)*(FN(KZ) - FN(IX))/ET

	ELSEIF (KZ.GT.IX .AND. IX.LE.IV .AND. KZ.LE.IV) THEN

        NT = KZ-IX+1
 
	GB1(NT) = DXTZ(KZ)*(PNSMBLK(IX,KZ)*(FN(KZ) - FN(IX))/ET)

	ENDIF

        IF (NPO.EQ.1) THEN
C
C     These terms are only for H and \tilde H nonzero
C

	IF (IX.LE.IV .AND. KZ.GT.IV) THEN
C
C     The relative signs, - for the quark kernel und + for the gluon switch for the
C     polarized case, as well as the signs infront of the QG and GQ and quark 
C     singlet piece.
C 

	IK = KZ - IV + 1

         IF (IKNL.EQ.2) THEN

	GB10(IK) = DXTZ(KZ)*(PNSMBLK(IX,KZ)+PNSMBLKD(IX,KZ+1))*FN(KZ)/ET
        
         ELSE

        GB10(IK) = DXTZ(KZ)*(PNSMBLK(IX,KZ)-PNSMBLKD(IX,KZ+1))*FN(KZ)/ET
        
         ENDIF

        ELSEIF (IX.LE.IV.AND.KZ.EQ.IV) THEN

        IK = KZ - IV + 1

         IF (IKNL.EQ.2) THEN

	GB10(IK) = DIFFDEL*(PNSMBLK(IX,KZ)+PNSMBLKD(IX,KZ+1))*FN(KZ)/ET
        
         ELSE

        GB10(IK) = DIFFDEL*(PNSMBLK(IX,KZ)-PNSMBLKD(IX,KZ+1))*FN(KZ)/ET
        
         ENDIF
		
	ELSEIF (IX.GT.IV .AND. KZ.GT.IX .AND. IX.LT.NX-1) THEN

	IW1 = KZ - IV
        IW = KZ - IX + 1

        G1(IW) = DXTZ(KZ)*PNSM(NR,IW1)*(FQ(KZ) - FQ(IX))/ET4 

	  ENDIF

	 ENDIF

        ENDIF

   31   CONTINUE


	IF (IX.LE.IV) THEN

C
C     Specify number of points in the various integration regions.
C

           NC = IX
           NC1 = IV-IX+1
           NC2 = NX-IV+1

C
C     Do the integration from 0..x, x..d, d..1 using Simpson
C

        IF (IX.LE.6) THEN

        TBLQQX = 0.0

	TBLQQ1 = SMPNOL(NC1,DZ,GB1,ERR)

        ELSEIF (IX.GE.IV-5) THEN

	TBLQQX = SMPNOR(NC,DZ,GBX,ERR)

	TBLQQ1 = 0.0

        ELSE

        TBLQQX = SMPNOR(NC,DZ,GBX,ERR)

	TBLQQ1 = SMPNOL(NC1,DZ,GB1,ERR)

        ENDIF

        IF (NPO.EQ.1) THEN

           TBLQQD1 = SMPSNA(NC2,DZ,GB10,ERR)

        ELSE

           TBLQQD1 = 0D0

        ENDIF

C
C	Save integration results in temp. array!
C

	TQQX(IX) = TBLQQX

	TQQ1(IX) = TBLQQ1

	TQQ11(IX) = TBLQQD1

	ELSE

      IF (NPO.EQ.1) THEN

C
C     Number of points in the integral
C

        NC3 = NX-IX
C
C     Value of the integrand at y=x. Only dummy value for book keeping purposes.
C

        G1(1) = 0.0
        G1(NX-IX) = -DXTZ(NX)*PNSM(NR,NX-IV)*FQ(IX)/ET4

C
C     Do integration in DGLAP region using Simpson
C
        IF (IX.LT.NX-3) THEN

	TEM = SMPNOL(NC3,DZ,G1,ERR)
      
        ELSE

           TEM = 0.0

        ENDIF

C	Save results in temp. array!

	TQQD(IX) = TEM

         ELSE

        TQQD(IX) = 0D0

         ENDIF

	ENDIF
	
   21 CONTINUE
C

	IY = 6

	DO 9999 IM = 1,NX-1

C
C	Successive extrapolation of missing integration results for the 0..x integral 
C	(first) and the x..del integral (second). In both cases 6 bins. are missing. 
C	Extrapolation uses 8th. order polynomial extrapolation.
C

	IF (IM.LE.6) THEN
           
           KY = 8

	DO 998 IE = 1,KY
	X2(IE) = XV(IY+IE)
 998	CONTINUE
	XX = XV(IY)

	DO 997 IE = 1,KY
	FX(IE) = TQQX(IY+IE)
 997	CONTINUE

	CALL RATINT(X2,FX,KY,XX,TEM,ERR)

	TQQX(IY) = TEM

	IY = IY - 1

	ELSEIF(IM.GE.IV-5 .AND.IM.LE.IV) THEN

           KU = 8

	DO 831 IE = 1,KU

	X2(IE) = XV(IM-IE)

 831	CONTINUE

	XX = XV(IM)

	DO 992 IE = 1,KU
	FX(IE) = TQQ1(IM-IE)
 992	CONTINUE

	CALL RATINT(X2,FX,KU,XX,TEM,ERR)

	TQQ1(IM) = TEM

	ENDIF

C
C	Extrapolation of last 5 missing bins for DGLAP region, again using 10th. order 
C	polynomial approximation.
C
	IF (IM.GE.NX-3) THEN

           KT = 7

	DO 990 IH = 1,KT
	X2(IH) = XV(IM-IH)
 990	CONTINUE
	XX = XV(IM)

	DO 989 IH = 1,KT
	FX(IH) = TQQD(IM-IH)
 989	CONTINUE
	
	CALL RATINT(X2,FX,KT,XX,TEM,ERR)

	TQQD(IM) = TEM

	ENDIF

 9999   CONTINUE

        DO 999 IM= 1,NX-1

C	2nd. order results. BL region

	IF (IM.LE.IV) THEN

           IF (IKNL.EQ.2) THEN
	
	TMP = ((TQQX(IM)+TQQ1(IM)+TQQ11(IM))*DIM1
     >  + FN(IM)*ANSM(IM)) * CPL2

        ELSE

       TMP = ((TQQX(IM)+TQQ1(IM)+TQQ11(IM))*DIM1
     >  + FN(IM)*(ANSM(IM)-2./9.*(13./2. - PI2 + 2.*Z3))) * CPL2

       TMX = ABS(XV(IM)-DEL/2.)
        IF (TMX.LT.1E-10) TMP = 0.0        

        ENDIF

C	1st. + 2nd. order

        TEMPI2(IM) = FO(IM)
	FO(IM) = FO(IM) + TMP
	
	ELSE

C	2nd. order results DGLAP region.

        IF (NPO.EQ.1) THEN

        IF (IKNL.EQ.2) THEN

	TMP1 = (TQQD(IM) - FQ(IM)*ANSM(IM))*CPL2 

        TEMPI1(IM) = TMP1*XA(IM,-1)  

        ELSE

        TMP1 = (TQQD(IM) - FQ(IM)*(ANSM(IM) +
     >  2./9.*(13./2. - PI2 + 2.*Z3)))*CPL2 

        TEMPI1(IM) = TMP1*XA(IM,-1)

        ENDIF

C	1st. + 2nd. order

	FO(IM) = FO(IM) + TMP1*XA(IM,-1)

        ELSE

        FO(IM) = 0D0

        ENDIF

	ENDIF

 999	CONTINUE

        IF (NPO.EQ.1) THEN

            KT = 5

         DO 1973 IP = 2,0,-1

          DO 1972 IE = 1,KT
      
           X2(IE) = 1D0*(IE+1)

           FX1(IE) = FO(IV+IE+IP)

 1972     CONTINUE

         XX = 1D0
	
         CALL POLINT(X2,FX1,KT,XX,TEM,ERR)

        	FO(IV+IP) = TEM

 1973    CONTINUE

         KT2 = 1

        DO 2011 IK = KT2,1,-1

        IF (IKNL.EQ.2) THEN

        FO(IK) = FO(IV-IK+1)

        ELSE

        FO(IK) = -FO(IV-IK+1)

        ENDIF
         
 2011   CONTINUE

        ELSE

        FO(1) = 0D0
        FO(IV) = 0D0

        ENDIF

      ENDIF	

      RETURN 
C                        **************************** 
      END 
 
      SUBROUTINE NSRHSP (NPO,TT, NEFF, FI, FO) 
C 
C       Subroutine to compute the Right-Side of the Altarelli-Parisi Equation 
C       This copy applies to the "NS-plus" piece involving (Qrk + Qrk-bar) 
C 
C       IKNL = 1,  2 : 1st & 2nd order evolution of the unpolarized case; 
C             -1, -2 :             .....                  polarized  .. 
C 
C       Nx is the number of mesh-points, Tt is the Log Q variable. 
C       Y is the input distribution function defined on the mesh points 
C       F is the RHS value (which is also = dY/dt) defined on the same mesh 
C  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      LOGICAL LSTX 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2) 
      PARAMETER (MXQX = MXX*MXQ) 
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2,  
     >Z3 = 2.404113806319188) 
      PARAMETER (D0 = 0.0, D1 = 1.0)	 
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX) 
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX, DEL, IV 
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX) 
      COMMON / KRNL01 / AGG2 (0:MXX), ANSP (0:MXX), ANSM (0:MXX) 
      COMMON / KRN1ST / FF1(0:MXX,0:MXX), GG1(0:MXX,0:MXX),    
     > FG2(0:MXX,0:MXX),GF2(0:MXX,0:MXX),GG2(0:MXX,0:MXX), 
     > PNSP(0:MXX,0:MXX), PNSM(0:MXX,0:MXX), SFF2(0:MXX,0:MXX), 
     > SFF2BLK(0:MXX,0:MXX), FG2BLK(0:MXX,0:MXX), GF2BLK(0:MXX,0:MXX),  
     > GG2BLK(0:MXX,0:MXX), PNSPBLK(0:MXX,0:MXX), PNSMBLK(0:MXX,0:MXX), 
     > GG2BLC1(0:MXX,0:MXX), PNSPBLC1(0:MXX,0:MXX),FF1BLX(0:MXX,0:MXX),  
     > GG2BLCX(0:MXX,0:MXX), PNSPBLCX(0:MXX,0:MXX),GG1BLX(0:MXX,0:MXX), 
     > PNSMBLC1(0:MXX,0:MXX), PNSMBLCX(0:MXX,0:MXX),FF1BL1(0:MXX,0:MXX), 
     > GG1BL1(0:MXX,0:MXX),SFF2BLKD(0:MXX,0:MXX), FG2BLKD(0:MXX,0:MXX), 
     > GF2BLKD(0:MXX,0:MXX), GG2BLKD(0:MXX,0:MXX),PNSPBLKD(0:MXX,0:MXX), 
     > PNSMBLKD(0:MXX,0:MXX),FG2BLK1(0:MXX,0:MXX),GF2BLK1(0:MXX,0:MXX),
     > SFF2BLK1(0:MXX,0:MXX),PNSPC(0:MXX,0:MXX),PNSMC(0:MXX,0:MXX),
     > GG2C(0:MXX,0:MXX), PNSPDGBL(0:MXX,0:MXX), PNSMDGBL(0:MXX,0:MXX),
     > GG2DGBL(0:MXX,0:MXX)

      COMMON / REGFUNC / GB1(MXX), GBX(MXX), GB21(MXX), GB510(MXX),  
     > GB2X(MXX), GB31(MXX), GB41(MXX), GB5X(MXX), GB51(MXX),  
     > GB3X(MXX), GB4X(MXX), GB210(MXX), GB10(MXX), GB310(MXX), 
     > G1(MXX), G2(MXX), G3(MXX), G4(MXX), G5(MXX), GB410(MXX)

      COMMON / DIFFD / DIFFDEL
      
c      COMMON / SECDERIV / DGBX(3,MXX), DGB2X(3,MXX), DGB1(3,MXX), 
c     > DGB21(3,MXX),
c     > DGB01(3,MXX), DGB41(3,MXX), DGB10(3,MXX), DGB210(3,MXX), 
c     > DG1(3,MXX), DG2(3,MXX)

      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx

      DIMENSION FI(NX), FO(NX), FN(MXX), FFN(3*MXX),FQ(MXX)
      DIMENSION W0(MXX), W1(MXX), WH(MXX), W2(MXX), W01(MXX)
      DIMENSION W3(MXX), W5(MXX), W22(MXX), FFI(3*MXX), WH1(MXX)
      DIMENSION X2(10), X1(10), FX(10), FX1(10), AF2(MXX), AF(MXX)
     > ,TQQ1(MXX), TQQX(MXX),TQQ11(MXX),TQQD(MXX)
c     > T3(0:MXX),T4(0:MXX)
      DIMENSION TEMPI1(0:MXX),TEMPI2(0:MXX)

	EXTERNAL SMPNOL,SMPNOR,SMPSNA

	SAVE 
	DATA AERR, RERR / 0.00001, 0.0002/ 


C
      S = EXP(TT)
      Q = AL * EXP (S)
      CPL = ALPI(Q)
      CPL2= CPL ** 2 / 2. * S
      CPL = CPL * S
C

      DZ = 1./(NX - 1)
	

        FN(NX) = 0.0
	FQ(NX) = 0.0

	DO 229 IY = 1, NX-1
        FQ(IY) = FI(IY)*XA(IY,1)
        TMX = ABS(XV(IY)-DEL/2.)
        IF (TMX.GT.1E-10) then
	FN(IY) = FI(IY)        
        ELSE
           IF (IKNL.GT.0) then
              FN(IY) = 0.0
           ELSE
             FN(IY) = FI(IY) 
           ENDIF
        ENDIF   

 229    CONTINUE 

C 
C	LO convolution of kernel with pdf! First unpol. then pol. case! 
C 
      CALL NEWARRAY (NX, DEL, FQ, FFI) 
      CALL NEWARRAY (NX, DEL, FN, FFN) 
      CALL INTEGR (NX, 0, FQ, W0, IR1) 
      CALL INTEGR (NX, 0, FN, W01, IR12) 
      CALL NINTEGR (NX, DEL, 1, FFI, FQ, W1, IR11) 
      CALL NINTEGR (NX, DEL, 15, FFN, FN, W2, IR01) 
      CALL NINTEGR (NX, DEL, 8, FFI, FQ, W22, IR02) 
      CALL NINTEGR (NX, DEL, 7, FFI, FQ, W3, IR22) 
      CALL NINTEGR (NX, DEL, 12, FFN, FN, W5, IR3) 
      CALL NINTEGR (NX, DEL, 16, FFN, FN, AF, IR4)
      CALL NINTEGR (NX, DEL, 20, FFN, FN, AF2, IR5)
      CALL HINTEGN (2,NX,IV,DEL,FFI,FQ, WH) 
      CALL HINTEGN (1,NX,IV,DEL,FFN,FN, WH1) 	

C 
	DIM1 = 1D0/DEL 

      DO 230 IZ = 2, NX 
      
      IF (XA(IZ,1) .LE. DEL) THEN 

	AA1 = 1D0 - XA(IZ,1)*DIM1 
         
        IF (IZ.EQ.IV) THEN

          IF (NPO.EQ.1) THEN

           IF (IKNL.GT.0) THEN

         DO 11 IL = 1,7
c         X2(IL) = XV(IL+1)
         X2(IL) = 1D0*IL
         FX(IL) = FO(IL+1)
         
 11      CONTINUE

c         XX = XV(1)
         XX = 0D0

         CALL RATINT(X2,FX,7,XX,TEM,ERR)

         FO(1) = TEM

        FO(IZ) = -FO(1) 

        ELSE

          DO 1 IL = 1,7
c         X2(IL) = XV(IL+1)
         X2(IL) = 1D0*IL
         FX(IL) = FO(IL+1)
         
 1       CONTINUE

c         XX = XV(1)

         XX = 0D0

         CALL RATINT(X2,FX,7,XX,TEM,ERR)

         FO(1) = TEM 

        FO(IZ) = FO(1)

        ENDIF
           
         ELSE

          FO(1) = 0D0
          FO(IV) = 0D0

         ENDIF

           ELSE

       IF (NPO.EQ.1) THEN

        IF (IKNL.GT.0) THEN

      FO(IZ) = 2.*FN(IZ) + 4./3.*(WH1(IZ) - AF(IZ)* 
     > AA1 + AA1*W5(IZ) + W2(IZ) - (AF2(IZ) - (1D0-AA1)*AF(IV)))

      TMX = ABS(XV(IZ)-DEL/2.)
        IF (TMX.LT.1E-10) FO(IZ) = 0.0

        ELSE

      FO(IZ) = 2.*FN(IZ) + 4./3.*(WH1(IZ) - AF(IZ)* 
     > AA1 + AA1*W5(IZ) + W2(IZ) + (AF2(IZ) - (1D0-AA1)*AF(IV)))

        ENDIF
        
        ELSE

        IF (IKNL.GT.0) THEN

      FO(IZ) = 2.*FN(IZ) + 4./3.*(WH1(IZ) - AF(IZ)* 
     > AA1 + AA1*W5(IZ) + W2(IZ))

      TMX = ABS(XV(IZ)-DEL/2.)
        IF (TMX.LT.1E-10) FO(IZ) = 0.0

        ELSE

      FO(IZ) = 2.*FN(IZ) + 4./3.*(WH1(IZ) - AF(IZ)* 
     > AA1 + AA1*W5(IZ) + W2(IZ))

        ENDIF

        ENDIF

         FO(IZ) = CPL*FO(IZ)
      
        ENDIF
      
      ELSE

      IF (NPO.EQ.1) THEN
       
      FO(IZ) = 2.* FQ(IZ) + 4./3.* (2D0*WH(IZ) - W0(IZ)  
     >         - W1(IZ) - W3(IZ) + W22(IZ))

      FO(IZ) = CPL * FO(IZ)*XA(IZ,-1) 

      ELSE

       FO(IZ) = 0D0

       ENDIF

	ENDIF 

  230 CONTINUE

C
C     Now NLO!
C

      IF (IKNL .EQ. 2 .OR. IKNL.EQ.-2) THEN

           IF (IKNL.EQ.2) THEN

         ET = 1D0
	ET1 = 1D0
        ET4 = 1D0
        ET2 = 1D0  

        ELSE

         ET = 1D0
	ET1 = 1D0
        ET4 = 1D0
        ET2 = 1D0  

        ENDIF

      DO 21 IX = 2, NX-1

        IF (IX.EQ.IV) GOTO 21

        X = XV(IX)

        NR = IX - IV

C
C       Evaluate the integrand for the I-th integral
C
C	First the BL integrands and then the DGLAP integrands.
C
C

        GB10(NX-IV+1) = 0.0

	IF (IX.LE.IV) THEN

	IQ = 1

	IQ1 = NX

	ELSE

	IQ = IX

	IQ1 = NX

	ENDIF
	
        DO 31 KZ = IQ, IQ1

	IF (IX.EQ.KZ) THEN

           IF (IX.LE.IV) THEN
C
C     Values of the integrand at y=x. Dummy values only for book keeping purposes.
C
        GBX(1) = 0.0

        GB1(1) = 0.0

        ENDIF
           
	GOTO 31

	ELSE

	IF (IX.GT.KZ .AND. IX.LE.IV) THEN

	GBX(KZ) = DXTZ(KZ)*PNSPBLK(IX,KZ)*(FN(KZ) - FN(IX))/ET

	ELSEIF (KZ.GT.IX .AND. IX.LE.IV .AND. KZ.LE.IV) THEN

        NT = KZ-IX+1
 
	GB1(NT) = DXTZ(KZ)*(PNSPBLK(IX,KZ)*(FN(KZ) - FN(IX))/ET)

	ENDIF

        IF (NPO.EQ.1) THEN

	IF (IX.LE.IV .AND. KZ.GT.IV) THEN
C
C     The relative signs, - for the quark kernel und + for the gluon switch for the
C     polarized case, as well as the signs infront of the QG and GQ and quark 
C     singlet piece.
C 

	IK = KZ - IV + 1

        IF (IKNL.EQ.2) THEN

	GB10(IK) = DXTZ(KZ)*(PNSPBLK(IX,KZ)-PNSPBLKD(IX,KZ+1))*FN(KZ)/ET
        
        ELSE

        GB10(IK) = DXTZ(KZ)*(PNSPBLK(IX,KZ)+PNSPBLKD(IX,KZ+1))*FN(KZ)/ET
        
        ENDIF

        ELSEIF (IX.LE.IV.AND.KZ.EQ.IV) THEN

        IK = KZ - IV + 1

        IF (IKNL.EQ.2) THEN

	GB10(IK) = DIFFDEL*(PNSPBLK(IX,KZ)-PNSPBLKD(IX,KZ+1))*FN(KZ)/ET
        
        ELSE

        GB10(IK) = DIFFDEL*(PNSPBLK(IX,KZ)+PNSPBLKD(IX,KZ+1))*FN(KZ)/ET
        
        ENDIF
		
	ELSEIF (IX.GT.IV .AND. KZ.GT.IX .AND. IX.LT.NX-1) THEN

	IW1 = KZ - IV
        IW = KZ - IX + 1

        G1(IW) = DXTZ(KZ)*PNSP(NR,IW1)*(FQ(KZ) - FQ(IX))/ET4 

	ENDIF

        ENDIF

	ENDIF

   31   CONTINUE


	IF (IX.LE.IV) THEN

C
C     Specify number of points in the various integration regions.
C

           NC = IX
           NC1 = IV-IX+1
           NC2 = NX-IV+1

C
C     Do the integration from 0..x, x..d, d..1 using Simpson
C

        IF (IX.LE.6) THEN

        TBLQQX = 0.0

	TBLQQ1 = SMPNOL(NC1,DZ,GB1,ERR)

        ELSEIF (IX.GE.IV-5) THEN

	TBLQQX = SMPNOR(NC,DZ,GBX,ERR)

	TBLQQ1 = 0.0

        ELSE

        TBLQQX = SMPNOR(NC,DZ,GBX,ERR)

	TBLQQ1 = SMPNOL(NC1,DZ,GB1,ERR)

        ENDIF

        IF (NPO.EQ.1) THEN

        TBLQQD1 = SMPSNA(NC2,DZ,GB10,ERR)

        ELSE

        TBLQQD1 = 0D0

        ENDIF
C
C	Save integration results in temp. array!
C

	TQQX(IX) = TBLQQX

	TQQ1(IX) = TBLQQ1

	TQQ11(IX) = TBLQQD1


	ELSE

      IF (NPO.EQ.1) THEN

C
C     Number of points in the integral
C

        NC3 = NX-IX
C
C     Value of the integrand at y=x. Only dummy value for book keeping purposes.
C

        G1(1) = 0.0
        G1(NX-IX) = -DXTZ(NX)*PNSP(NR,NX-IV)*FQ(IX)/ET4 
C
C     Do integration in DGLAP region using Simpson
C
        IF (IX.LT.NX-5) THEN

	TEM = SMPNOL(NC3,DZ,G1,ERR)

        ELSE

        TEM = 0.0

        ENDIF
      
C	Save results in temp. array!

	TQQD(IX) = TEM

        ELSE

        TQQD(IX) = 0D0

        ENDIF

	ENDIF
	
   21 CONTINUE
C

	IY = 6
c        IY1 = 5

	DO 9999 IM = 1,NX-1

C
C	Successive extrapolation of missing integration results for the 0..x integral 
C	(first) and the x..del integral (second). In both cases 6 bins. are missing. 
C	Extrapolation uses 8th. order polynomial extrapolation.
C

	IF (IM.LE.6) THEN
           
           KY = 8

	DO 998 IE = 1,KY
	X2(IE) = XV(IY+IE)
 998	CONTINUE
	XX = XV(IY)

	DO 997 IE = 1,KY
	FX(IE) = TQQX(IY+IE)
 997	CONTINUE

	CALL RATINT(X2,FX,KY,XX,TEM,ERR)

	TQQX(IY) = TEM

	IY = IY - 1

	ELSEIF(IM.GE.IV-5 .AND.IM.LE.IV) THEN

           KU = 8

	DO 831 IE = 1,KU
	X2(IE) = XV(IM-IE)
 831	CONTINUE
	XX = XV(IM)

	DO 992 IE = 1,KU
	FX(IE) = TQQ1(IM-IE)
 992	CONTINUE

	CALL RATINT(X2,FX,KU,XX,TEM,ERR)

	TQQ1(IM) = TEM

	ENDIF
C
C	Extrapolation of last 5 missing bins for DGLAP region, again using 7th. order 
C	polynomial approximation.
C
	IF (IM.GE.NX-3) THEN

           KT = 7

	DO 990 IH = 1,KT
	X2(IH) = XV(IM-IH)
 990	CONTINUE
	XX = XV(IM)

	DO 989 IH = 1,KT
	FX(IH) = TQQD(IM-IH)
 989	CONTINUE
	
	CALL RATINT(X2,FX,KT,XX,TEM,ERR)

	TQQD(IM) = TEM

	ENDIF

 9999   CONTINUE

C	2nd. order results. BL region

        DO 999 IM = 1,NX-1

	IF (IM.LE.IV) THEN

           IF (IKNL.EQ.2) THEN
	
	TMP = ((TQQX(IM)+TQQ1(IM)+TQQ11(IM))*DIM1
     >  - FN(IM)*(2./9.*(13./2. - PI2 + 2.*Z3)-ANSP(IM))) * CPL2

        TMX = ABS(XV(IM)-DEL/2.)
        IF (TMX.LT.1E-10) TMP = 0.0

        ELSE

       TMP = ((TQQX(IM)+TQQ1(IM)+TQQ11(IM))*DIM1
     >  + FN(IM)*ANSP(IM)) * CPL2

        ENDIF

C	1st. + 2nd. order

        TEMPI2(IM) = FO(IM)

	FO(IM) = FO(IM) + TMP
	
	ELSE

        IF (NPO.EQ.1) THEN

C	2nd. order results DGLAP region.

        IF (IKNL.EQ.2) THEN

	TMP1 = (TQQD(IM) - FQ(IM)*(ANSP(IM) +
     >  2./9.*(13./2. - PI2 + 2.*Z3)))*CPL2 

        TEMPI1(IM) = TMP1*XA(IM,-1)

        ELSE

        TMP1 = (TQQD(IM) - FQ(IM)*ANSP(IM))*CPL2 

         TEMPI1(IM) = TMP1*XA(IM,-1)

        ENDIF

C	1st. + 2nd. order

	FO(IM) = FO(IM) + TMP1*XA(IM,-1)

        ELSE

         FO(IM) = 0D0

        ENDIF

	ENDIF
         
 999	CONTINUE

        IF (NPO.EQ.1) THEN

            KT = 5

         DO 1973 IP = 2,0,-1

          DO 1972 IE = 1,KT
      
           X2(IE) = 1D0*(IE+1)

           FX1(IE) = FO(IV+IE+IP)

 1972     CONTINUE

         XX = 1D0
	
         CALL POLINT(X2,FX1,KT,XX,TEM,ERR)

        	FO(IV+IP) = TEM

 1973    CONTINUE



       KT2 = 1

       DO 2011 IK = KT2,1,-1

        IF (IKNL.EQ.2) THEN

        FO(IK) = -FO(IV-IK+1)

        ELSE

        FO(IK) = FO(IV-IK+1)

        ENDIF
         
 2011   CONTINUE

        ELSE

        FO(1) = 0D0
        FO(IV) = 0D0

        ENDIF

        ENDIF	

      RETURN 
C                        **************************** 
      END 

      SUBROUTINE SNEVL(NPO,XT,IV,DEL,IKNL, NX, NT, JT, DT, 
     >TIN, NEFF, UI, GI, US, GS)
C
C       This is the singlet counter-part of the NSEVL subroutine. Refer to
C                       comments at the beginning of that program section.
C
C     Input parton distributions are Gi (for gluon) and Ui (for singlet quark)

C                               at Tt = 0; output distributions are Gs and Us
C                                       at Tt = IS*dt with IS = 1, 2, ... , Nt.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6)
      PARAMETER (MXQX= MXQ * MXX)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      Common / PdCntrl/ LPrt, LDbg
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)      
      
C
      DIMENSION UI(NX), US(0:NX, 0:NT)
      DIMENSION GI(NX), GS(0:NX, 0:NT)
      DIMENSION Y0(MXX), Y1(MXX), YP(MXX), F0(MXX), F1(MXX), FP(MXX)
      DIMENSION Z0(MXX), Z1(MXX), ZP(MXX), G0(MXX), G1(MXX), GP(MXX)
      DIMENSION XT(0:MXX)

      DATA D0 / 0.0 /

         If (Ldbg .Eq. 1) Then
            Write (Nwrt, '(A)') 'Singlet:'
            N5 = Nx / 5 + 1
         Endif

C                                       Faster evolution of the singlet sector
C                                       requires finer iteration to achieve the
C                                       same accuracy  as in the non-singlet.
      JTT = 2 * JT
      DDT = DT / JTT
C
      IF (NX .GT. MXX) THEN
      WRITE (NOUT,*) 'Nx =', NX, ' greater than Max # of pts in SNEVL.'
      STOP 'Program stopped in SNEVL'
      EndIf
C                                       Compute an effective first order lamda
C                                       to be used in checking of moment evl.
      TMD = TIN + DT * NT / 2.
      AMU = EXP(TMD)
      TEM = 6./ (33.- 2.* NEFF) / ALPI(AMU)
      TLAM = TMD - TEM

C                                      Initialization (see previous subroutine)
      DO 9 IX = 1, NX
      US (IX, 0) = UI(IX)
      GS (IX, 0) = GI(IX)
    9 CONTINUE
      US ( 0, 0) = (UI(1) - UI(2))* 3D0 + UI(3)
      GS ( 0, 0) = (GI(1) - GI(2))* 3D0 + GI(3)
C
      TT = TIN
      DO 10 IZ = 1, NX
      Y0(IZ) = UI(IZ)
      Z0(IZ) = GI(IZ)
   10 CONTINUE

C     loop in the Tt variable


      DO 20 IS = 1, NT

         IPO = JTT
C                                                       fine- grained iteration
         DO 202 JS = 1, JTT
C                                          Irnd is the counter for Q-iterations
            IRND = (IS-1) * JTT + JS
C                                               Use Runge-Katta the first round

            IF (IRND .EQ. 1) THEN
C
                CALL SNRHS (NPO,TT, NEFF, Y0,Z0,  F0,G0)
C
c     smoothing as in NSEVL
C
                CALL SMOOTH (0,IKNL,MXX,DEL,IV,XT,F0)
                CALL SMOOTH (1,IKNL,MXX,DEL,IV,XT,G0)

                DO 250 IZ = 1, NX
                   Y0(IZ) = Y0(IZ) + DDT * F0(IZ)
                   Z0(IZ) = Z0(IZ) + DDT * G0(IZ)
  250           CONTINUE

C
c     smoothing

              CALL SMOOTH (0,IKNL,MXX,DEL,IV,XT,Y0)
              CALL SMOOTH (1,IKNL,MXX,DEL,IV,XT,Z0)  

                TT = TT + DDT
C
                CALL SNRHS (NPO,TT, NEFF, Y0, Z0,  F1, G1)
C
c     smoothing
              CALL SMOOTH (0,IKNL,MXX,DEL,IV,XT,F1)
              CALL SMOOTH (1,IKNL,MXX,DEL,IV,XT,G1) 
                
                DO 251 IZ = 1, NX
                   Y1(IZ) = UI(IZ) + DDT * (F0(IZ) + F1(IZ)) / 2D0
                   Z1(IZ) = GI(IZ) + DDT * (G0(IZ) + G1(IZ)) / 2D0
  251           CONTINUE
C     smoothing

              CALL SMOOTH (0,IKNL,MXX,DEL,IV,XT,Y1)
              CALL SMOOTH (1,IKNL,MXX,DEL,IV,XT,Z1) 


C                            What follows is a combination of the 2-step method
C                                   and the Adams Predictor-Corrector Algorithm
            Else
C
                CALL SNRHS (NPO,TT, NEFF, Y1, Z1,  F1, G1)
C     smoothing

              CALL SMOOTH (0,IKNL,MXX,DEL,IV,XT,F1)
              CALL SMOOTH (1,IKNL,MXX,DEL,IV,XT,G1)
C                                                                     Predictor

                DO 252 IZ = 1, NX
                   YP(IZ) = Y1(IZ) + DDT * (3D0 * F1(IZ) - F0(IZ)) / 2D0
                   ZP(IZ) = Z1(IZ) + DDT * (3D0 * G1(IZ) - G0(IZ)) / 2D0
  252           CONTINUE
c     smoothing

              CALL SMOOTH (0,IKNL,MXX,DEL,IV,XT,YP)
              CALL SMOOTH (1,IKNL,MXX,DEL,IV,XT,ZP)

C                        Increment of Tt at this place is part of the formalism
                TT = TT + DDT
C
                CALL SNRHS (NPO,TT, NEFF, YP, ZP,  FP, GP)
C
c                smoothing

              CALL SMOOTH (0,IKNL,MXX,DEL,IV,XT,FP)
              CALL SMOOTH (1,IKNL,MXX,DEL,IV,XT,GP)

C                                                                     Corrector
                DO 253 IZ = 1, NX
                   Y1(IZ) = Y1(IZ) + DDT * (FP(IZ) + F1(IZ)) / 2D0
                   Z1(IZ) = Z1(IZ) + DDT * (GP(IZ) + G1(IZ)) / 2D0
                   F0(IZ) = F1(IZ)
                   G0(IZ) = G1(IZ)
  253           CONTINUE
C     smoothing

              CALL SMOOTH (0,IKNL,MXX,DEL,IV,XT,Y1)
              CALL SMOOTH (1,IKNL,MXX,DEL,IV,XT,Z1)

            EndIf
C
  202    CONTINUE
C                       Fill output array and restore factor of X, if necessary
C                                    For spin-averaged case, enforce positivity
         DO 260 IX = 1, NX
           IF (IKNL .GT. 0 .AND. XA(IX,1).GT.DEL) THEN
c            US (IX, IS) = MAX(Y1(IX), D0)
c            GS (IX, IS) = MAX(Z1(IX), D0)
            US (IX, IS) = Y1(IX)
            GS (IX, IS) = Z1(IX)
           Else
            US (IX, IS) = Y1(IX)
            GS (IX, IS) = Z1(IX)
           EndIf
  260    CONTINUE
C
C               The value of the function at x=0 is obtained by extrapolation
C
         US(0, IS) = 3D0*Y1(1) - 3D0*Y1(2) + Y1(3)
         GS(0, IS) = 3D0*Z1(1) - 3D0*Z1(2) + Z1(3)
C
C                                                    Print out for Debugging
      If (LDbg .Eq. 1) Then
         Write (Nwrt, '(A, 5(1pE12.3))') ' SQ:',(Us(Iz,Is), Iz=1,Nx,N5) 
         Write (Nwrt, '(A, 5(1pE12.3))') '  G:',(Gs(Iz,Is), Iz=1,Nx,N5)
      Endif

      JTT = IPO

   20 CONTINUE
C
      RETURN
C                        ****************************
      END
      
	SUBROUTINE SNRHS (NPO,TT, NEFF, FI, GI,  FO, GO) 
C 
C       Subroutine to compute the Right-Side of the Altarelli-Parisi Equation 
C                                                       for the Singlet sector: 
C       See comments in NSRHSP for notes on IKNL, and programm comments! 
C 
C       FI, Z are the input distributions for quark and gluon respectively; 
C       FO, G are the output dY/dt, dZ/dt. 
C       Nx is the number of mesh-points, Tt is the Log Q variable. 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      LOGICAL LSTX 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2) 
      PARAMETER (MXQX = MXX*MXQ) 
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2,  
     >Z3 = 2.404113806319188) 
      PARAMETER (D0 = 0.0, D1 = 1.0)	 
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX) 
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX, DEL, IV 
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX) 
      COMMON / KRNL01 / AGG2 (0:MXX), ANSP (0:MXX), ANSM (0:MXX) 
      COMMON / KRN1ST / FF1(0:MXX,0:MXX), GG1(0:MXX,0:MXX),    
     > FG2(0:MXX,0:MXX),GF2(0:MXX,0:MXX),GG2(0:MXX,0:MXX), 
     > PNSP(0:MXX,0:MXX), PNSM(0:MXX,0:MXX), SFF2(0:MXX,0:MXX), 
     > SFF2BLK(0:MXX,0:MXX), FG2BLK(0:MXX,0:MXX), GF2BLK(0:MXX,0:MXX),  
     > GG2BLK(0:MXX,0:MXX), PNSPBLK(0:MXX,0:MXX), PNSMBLK(0:MXX,0:MXX), 
     > GG2BLC1(0:MXX,0:MXX), PNSPBLC1(0:MXX,0:MXX),FF1BLX(0:MXX,0:MXX),  
     > GG2BLCX(0:MXX,0:MXX), PNSPBLCX(0:MXX,0:MXX),GG1BLX(0:MXX,0:MXX), 
     > PNSMBLC1(0:MXX,0:MXX), PNSMBLCX(0:MXX,0:MXX),FF1BL1(0:MXX,0:MXX), 
     > GG1BL1(0:MXX,0:MXX),SFF2BLKD(0:MXX,0:MXX), FG2BLKD(0:MXX,0:MXX), 
     > GF2BLKD(0:MXX,0:MXX), GG2BLKD(0:MXX,0:MXX),PNSPBLKD(0:MXX,0:MXX), 
     > PNSMBLKD(0:MXX,0:MXX),FG2BLK1(0:MXX,0:MXX),GF2BLK1(0:MXX,0:MXX),
     > SFF2BLK1(0:MXX,0:MXX),PNSPC(0:MXX,0:MXX),PNSMC(0:MXX,0:MXX),
     > GG2C(0:MXX,0:MXX), PNSPDGBL(0:MXX,0:MXX), PNSMDGBL(0:MXX,0:MXX),
     > GG2DGBL(0:MXX,0:MXX)
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx 
      COMMON / REGFUNC / GB1(MXX), GBX(MXX), GB21(MXX), GB510(MXX),  
     > GB2X(MXX), GB31(MXX), GB41(MXX), GB5X(MXX), GB51(MXX),  
     > GB3X(MXX), GB4X(MXX), GB210(MXX), GB10(MXX), GB310(MXX), 
     > G1(MXX), G2(MXX), G3(MXX), G4(MXX), G5(MXX), GB410(MXX)

      COMMON / DIFFD / DIFFDEL

c      COMMON / SECDERIV / DGBX(3,MXX), DGB2X(3,MXX), DGB1(3,MXX), 
c     > DGB21(3,MXX),
c     > DGB01(3,MXX), DGB41(3,MXX), DGB10(3,MXX), DGB210(3,MXX), 
c     > DG1(3,MXX), DG2(3,MXX)

      DIMENSION GI(NX), GO(NX)
      DIMENSION FI(NX), FO(NX), W0(MXX), W1(MXX), WH(MXX), R4(MXX) 
      DIMENSION R0(MXX), R1(MXX), RH(MXX),W2(MXX), R3(MXX) 
      DIMENSION FII(3*MXX), GII(3*MXX), R5(MXX), R6(MXX), RH1(MXX) 
      DIMENSION WH1(MXX), WH2(MXX), W21(MXX), R51(MXX), R61(MXX) 
      DIMENSION W3(MXX), W7(MXX), W8(MXX), W9(MXX), W11(MXX) 
      DIMENSION W22(MXX), RH3(MXX), FN(MXX), FIN(3*MXX) 
      DIMENSION W01(MXX), WH0(MXX), FQ(MXX),RHX(MXX)
      DIMENSION X2(10), FX(10), GX(10),TQQX(MXX)
     > , TGGX(MXX), TQQ1(MXX), TGG1(MXX), TQQD(MXX) 
     > , TGGD(MXX), TQQ11(MXX),FX1(10)
     > , TGG11(MXX), AF1(MXX), AF2(MXX), AG1(MXX)
     > , AG2(MXX), AG3(MXX), AG4(MXX), AF(MXX), AF0(MXX)
     
      DIMENSION TEMPI1(MXX),TEMPI2(0:MXX),TEMPI3(0:MXX),TEMPI4(0:MXX)
      
      EXTERNAL SMPSNA,SMPNOR,SMPNOL
 
	SAVE 
	DATA AERR, RERR / 0.0001, 0.0005/ 
 
C	Exponentiated running alpha_s 
      S = EXP(TT) 
      Q = AL * EXP (S) 
      CPL = ALPI(Q) 
      CPL2= CPL ** 2 / 2. * S 
      CPL = CPL * S 
      DIM = 1D0/DEL
C   
 
	DZ = 1./(NX -1) 

C	Initialize quark and gluon pdfs at x=1. 
	FN(NX) = 0.0 
	FI(NX)= 0.0 
	GI(NX) = 0.0 
        ETC = 1D0

C	convert Q to q via q=Q*x.

	DO 228 IY = 1, NX-1 
        FQ(IY) = FI(IY)*XA(IY,1)   
        TMX = ABS(XV(IY)-DEL/2.)
        IF (TMX.GT.1E-10) then
	FN(IY) = FI(IY)/ETC        
        ELSE
           IF (IKNL.GT.0) then
              FN(IY) = 0.0
           ELSE
             GI(IY) = 0.0
             FN(IY) = FI(IY)/ETC 
           ENDIF
        ENDIF

 228    CONTINUE 

C
C	Start LO integration.
C 
 
      CALL NEWARRAY (NX, DEL, FQ, FII) 
      CALL NEWARRAY (NX, DEL, FN, FIN) 
      CALL NEWARRAY (NX, DEL, GI, GII) 
C 
      CALL NINTEGR (NX, DEL, 1, FII, FQ, W1, IR11)       
      CALL NINTEGR (NX, DEL, 8, FII, FQ, WH1, IR22) 
      CALL NINTEGR (NX, DEL, 15, FIN, FN, W22, IR20) 
      CALL NINTEGR (NX, DEL, 7, FII, FQ, WH2, IR21) 
      CALL NINTEGR (NX, DEL, 22, FIN, FN, AF1, IR58) 
C 
      CALL INTEGR (NX, 0, FQ, W0, IR2) 
      CALL INTEGR (NX, 0, FN, W01, IR19)  
      CALL INTEGR (NX, 0, GI, R0, IR5) 
      CALL INTEGR (NX, 1, GI, R1, IR6) 
 
C 
      CALL NINTEGR (NX, DEL, 8, GII, GI, RH1, IR43) 
      CALL NINTEGR (NX, DEL, 15, GII, GI, RH3, IR39) 
      CALL NINTEGR (NX, DEL, 14, FIN, FN, W11, IR48) 
      CALL NINTEGR (NX, DEL, 12, GII, GI, W8, IR46) 
      CALL NINTEGR (NX, DEL, 13, GII, GI, W7, IR47) 
      CALL NINTEGR (NX, DEL, 12, FIN, FN, W9, IR49) 
      CALL NINTEGR (NX, DEL, 16, FIN, FN, AF, IR50) 
      CALL NINTEGR (NX, DEL, 16, GII, GI, AG1, IR51) 
      CALL NINTEGR (NX, DEL, 17, GII, GI, AG2, IR52) 
      CALL NINTEGR (NX, DEL, 20, FIN, FN, AF2, IR54) 
      CALL NINTEGR (NX, DEL, 19, GII, GI, AG4, IR56) 
      CALL NINTEGR (NX, DEL, 21, FIN, FN, AF0, IR57) 
      
C
      CALL HINTEGN (2,NX,IV,DEL,FII,FQ, WH) 
      CALL HINTEGN (1,NX,IV,DEL,FIN,FN, WH0) 
      CALL HINTEGN (2,NX,IV,DEL,GII,GI, RH) 
      CALL HINTEGN (1,NX,IV,DEL,GII,GI, RHX) 
C 
       A = 1D0 - DEL/2D0

      IF (IKNL .GT. 0) THEN 
       
      CALL NINTEGR (NX, DEL, 2, GII, GI, R3, IR13) 
      CALL NINTEGR (NX, DEL, 3, GII, GI, R4, IR23) 
      CALL NINTEGR (NX, DEL, 6, FII, FQ, W2, IR01) 
      CALL NINTEGR (NX, DEL, 4, GII, GI, R5, IR33) 
      CALL NINTEGR (NX, DEL, 5, GII, GI, R6, IR42) 
      
     

      DO 230 IZ = 2, NX 
 
	IF (XA(IZ,1) .LE. DEL) THEN 
 
	AA1 = XA(IZ,1)*DIM - 1D0 
      
      IF (IZ.EQ.IV) THEN

C
C     set first and last point in ERBL region (NPO=1 evol. in ERBL+DGLAP, NPO=2
C     evol. only in ERBL.)
C
      IF (NPO.EQ.1) THEN

      DO 1 IL = 1,7
         X2(IL) = XV(IL+1)
         FX(IL) = FO(IL+1)
         GX(IL) = GO(IL+1)
 1       CONTINUE

         XX = XV(1)

         CALL RATINT(X2,FX,7,XX,TEM,ERR)

         FO(1) = TEM

         CALL RATINT(X2,GX,7,XX,TEM,ERR)

         GO(1) = TEM

      FO(IZ) = - FO(1)

      GO(IZ) = GO(1) 

      ELSE

      FO(1) = 0D0
      GO(1) = 0D0
      FO(IV) = 0D0
      GO(IV) = 0D0

      ENDIF

           ELSE
C
C     result of convolution of kernel with GPD. again split according to NPO value.
C     ERBL region.
C

      IF (NPO.EQ.1) THEN

      FO(IZ) = (2.*FN(IZ) + 4./3.*(WH0(IZ) +
     > AF(IZ)*AA1 - AA1*W9(IZ) + W22(IZ) 
     > - (AF2(IZ) - (AA1+1D0)*AF(IV))))
     > + NEFF*(AG2(IZ)*(2.*AA1 + 1.)*XA(IZ,1) 
     > + 4.*(AA1+1.)*AA1*AG1(IZ) - AA1* 
     > (2.*AA1 +1.)*W7(IZ) - 4.*(AA1 + 1.)
     > *AA1*W8(IZ) - DEL*AA1*(2D0*AA1+1D0)*AG2(IV)
     > - 4D0*AA1*(AA1+1D0)*AG1(IV))*DIM*A


      GO(IZ) = (33D0 - 2D0 * NEFF)/6D0*GI(IZ) + 3D0
     > *(RHX(IZ) + AG1(IZ)*(2D0*(AA1+1D0)**2*(1D0 - 2D0
     > *AA1) - 1D0)-AG2(IZ)*(AA1**2+(AA1+1D0)**2)*XA(IZ,1)+RH3(IZ)
     > + 2D0*W8(IZ)*AA1**2*(3D0+2D0*AA1) + ((AA1+1D0)**2
     > +AA1**2)*AA1*W7(IZ)
     > + AG4(IZ) + DEL*AA1*(AA1**2 + (AA1+1D0)**2)*AG2(IV)
     > + 2D0*AA1**2*(2D0*AA1 + 3D0)*AG1(IV)) 
     > + 4./3.*(-(1D0+AA1)*AF(IZ)*XA(IZ,1)
     > + W9(IZ)*DEL*AA1**2 - 2D0*W11(IZ)*AA1**2
     > *(3D0 + 2D0*AA1) - DEL*AA1**2*AF(IV) 
     > + 2D0*(AA1+1D0)**2*(1D0 - 2D0*AA1)*AF1(IZ) + 2D0*AA1**2*
     > (2D0*AA1+3D0)*AF0(IV))/A 

      ELSE

       FO(IZ) = (2.*FN(IZ) + 4./3.*(WH0(IZ) +
     > AF(IZ)*AA1 - AA1*W9(IZ) + W22(IZ)))
     > + NEFF*(AG2(IZ)*(2.*AA1 + 1.)*XA(IZ,1) 
     > + 4.*(AA1+1.)*AA1*AG1(IZ) - AA1* 
     > (2.*AA1 +1.)*W7(IZ) - 4.*(AA1 + 1.)
     > *AA1*W8(IZ))*DIM*A


      GO(IZ) = (33D0 - 2D0 * NEFF)/6D0*GI(IZ) + 3D0
     > *(RHX(IZ) + AG1(IZ)*(2D0*(AA1+1D0)**2*(1D0 - 2D0
     > *AA1) - 1D0)-AG2(IZ)*(AA1**2+(AA1+1D0)**2)*XA(IZ,1)+RH3(IZ)
     > + 2D0*W8(IZ)*AA1**2*(3D0+2D0*AA1) + ((AA1+1D0)**2
     > +AA1**2)*AA1*W7(IZ)) 
     > + 4./3.*(-(1D0+AA1)*AF(IZ)*XA(IZ,1)
     > + W9(IZ)*DEL*AA1**2 - 2D0*W11(IZ)*AA1**2
     > *(3D0 + 2D0*AA1) 
     > + 2D0*(AA1+1D0)**2*(1D0 - 2D0*AA1)*AF1(IZ))/A

      ENDIF

         FO(IZ) = CPL*FO(IZ)

         TMX = ABS(XV(IZ)-DEL/2.)
         IF (TMX.LT.1E-10) FO(IZ) = 0.0

      	 GO(IZ) = GO(IZ) * CPL   

        ENDIF

      else 

C
C     as above but now DGLAP region!
C
      IF (NPO.EQ.1) THEN

      FO(IZ) = 2D0*FQ(IZ) + 4D0/3D0*(2D0*WH(IZ)-W0(IZ)  
     > + WH1(IZ) - WH2(IZ) - W1(IZ))+NEFF*(R3(IZ) + R4(IZ))*A 

      FO(IZ) = XA(IZ,-1) * FO(IZ) * CPL 
      
      GO(IZ) = 4D0 / 3D0 * (W0(IZ) + W2(IZ))/A
     >   + (33D0 - 2D0 * NEFF) / 6D0 * GI(IZ) + 3D0  
     >   * (2D0*RH(IZ)  
     >   - R0(IZ) - R1(IZ) - 2D0*R4(IZ)   
     >   + RH1(IZ) + R6(IZ) + 2D0*R5(IZ)) 
      GO(IZ) = GO(IZ) * CPL 

      ELSE

      FO(IZ) = 0D0
      GO(IZ) = 0D0

      ENDIF

	ENDIF 


  230 CONTINUE 

      Else 

C
C     polarized case now. analogous to unpol.
C

       
      CALL NINTEGR (NX, DEL, 9, GII, GI, W21, IR99) 
      CALL NINTEGR (NX, DEL, 10, FII, FQ, R51, IR98) 
      CALL NINTEGR (NX, DEL, 11, GII, GI, R61, IR97) 
       
      DO 240 IZ = 2, NX 
 
	IF (XA(IZ,1) .LE. DEL) THEN 
 
	AA1 = XA(IZ,1)*DIM - 1D0  

      IF (IZ.EQ.IV) THEN

       IF (NPO.EQ.1) THEN

        DO 11 IL = 1,7
         X2(IL) = XV(IL+1)
         FX(IL) = FO(IL+1)
         GX(IL) = GO(IL+1)
 11   CONTINUE

         XX = XV(1)

         CALL RATINT(X2,FX,7,XX,TEM,ERR)

         FO(1) = TEM

         CALL RATINT(X2,GX,7,XX,TEM,ERR)

         GO(1) = TEM
      
       FO(IZ) = FO(1)

       GO(IZ) = -GO(1)

       ELSE

       FO(1) = 0D0
       FO(IV) = 0D0
       GO(1) = 0D0
       GO(IV) = 0D0 

       ENDIF

       ELSE

       IF (NPO.EQ.1) THEN
 
      FO(IZ) = 2.*FN(IZ) + 4./3.*(WH0(IZ) + AF(IZ)*
     > AA1 - AA1*W9(IZ) + W22(IZ)
     > + (AF2(IZ) - (AA1+1D0)*AF(IV)))
     > + NEFF*(-AG2(IZ)*XA(IZ,1) - W7(IZ)*AA1 + DEL*AA1*AG2(IV) 
     > )*DIM*A
 
      GO(IZ) = (33D0 - 2D0 * NEFF)/6D0*GI(IZ) + 3D0
     > *(RHX(IZ) - AG1(IZ) + AG2(IZ)*(2.*AA1 + 1.)*XA(IZ,1) 
     > + RH3(IZ) + (1. + 2.*AA1)*AA1*W7(IZ) - AG4(IZ) 
     > - DEL*AA1*(2D0*AA1+1D0)*AG2(IV))
     > + 4./3.*((1. + AA1)*AF(IZ)*XA(IZ,1) 
     > - W9(IZ)*DEL*AA1**2 - DEL*AA1**2*AF(IV))/A
 
      ELSE

       FO(IZ) = 2.*FN(IZ) + 4./3.*(WH0(IZ) + AF(IZ)*
     > AA1 - AA1*W9(IZ) + W22(IZ))
     > + NEFF*(-AG2(IZ)*XA(IZ,1) - W7(IZ)*AA1)*DIM*A
 
      GO(IZ) = (33D0 - 2D0 * NEFF)/6D0*GI(IZ) + 3D0
     > *(RHX(IZ) - AG1(IZ) + AG2(IZ)*(2.*AA1 + 1.)*XA(IZ,1) 
     > + RH3(IZ) + (1. + 2.*AA1)*AA1*W7(IZ))
     > + 4./3.*((1. + AA1)*AF(IZ)*XA(IZ,1) 
     > - W9(IZ)*DEL*AA1**2)/A

      ENDIF

       FO(IZ) = CPL*FO(IZ)

       GO(IZ) = GO(IZ) * CPL 

        TMX = ABS(XV(IZ)-DEL/2.)
        IF (TMX.LT.1E-10) GO(IZ) = 0.0

      ENDIF
      
      else 

      IF (NPO.EQ.1) THEN
 
      FO(IZ) = NEFF * W21(IZ)*A 
     > + 2.* FQ(IZ) + 4./ 3.* (2D0*WH(IZ)-W0(IZ)  
     >      + WH1(IZ) - WH2(IZ) - W1(IZ) ) 

      FO(IZ) = XA(IZ,-1) * FO(IZ) * CPL 
 
      GO(IZ) = 4./3. * R51(IZ)/A 
     >+ (33.- 2.* NEFF) / 6.* GI(IZ) + 3.*(2D0*RH(IZ)  
     >   - R0(IZ) - R1(IZ) + RH1(IZ) + R61(IZ)) 

      GO(IZ) = GO(IZ) * CPL 

      ELSE

      FO(IZ) = 0D0
      GO(IZ) = 0D0

      ENDIF

	ENDIF 
 
  240 CONTINUE 
 
      EndIf 
C
C     NOW NLO!
C

 	IF (IKNL .EQ. 2 .OR. IKNL.EQ.-2) THEN


           IF (IKNL.EQ.2) THEN

         ET = 1D0
	ET1 = 1D0
        ET4 = 1D0
        ET2 = 1D0  
        
        ELSE

         ET = 1D0
	ET1 = 1D0
        ET4 = 1D0
        ET2 = 1D0  

        ENDIF

      DO 21 IX = 2, NX-1


        IF (IX.EQ.IV) GOTO 21

        X = XV(IX)

        NR = IX - IV
C
C       Evaluate the integrand for the I-th integral
C
C	First the BL integrands and then the DGLAP integrands.
C
C
        GB10(NX-IV+1) = 0.0
        
	GB210(NX-IV+1) = 0.0
       


	IF (IX.LE.IV) THEN

	IQ = 1

	IQ1 = NX

	ELSE

	IQ = IX

	IQ1 = NX

	ENDIF
	
        DO 31 KZ = IQ, IQ1

	IF (IX.EQ.KZ) THEN

           IF (IX.LE.IV) THEN
C
C     Values of the integrand at y=x. Dummy values only for book keeping purposes.
C
        GBX(IX) = 0.0

        GB2X(IX) = 0.0

        GB1(1) = 0.0

	GB21(1) = 0.0

        ENDIF
           
	GOTO 31

	ELSE

	IF (IX.GT.KZ .AND. IX.LE.IV) THEN

	GBX(KZ) = DXTZ(KZ)*((PNSPBLK(IX,KZ)*(FN(KZ) - FN(IX))
     > + SFF2BLK(IX,KZ)*FN(KZ))*DIM + FG2BLK(IX,KZ)*GI(KZ)*DIM**2*A)
 
	GB2X(KZ) = DXTZ(KZ)*(GG2BLK(IX,KZ)*(GI(KZ) - GI(IX))*DIM
     > + GF2BLK(IX,KZ)*FN(KZ)/A) 

	ELSEIF (KZ.GT.IX .AND. IX.LE.IV .AND. KZ.LE.IV) THEN

        NT1 = KZ-IX+1
 
	GB1(NT1) = DXTZ(KZ)*((PNSPBLK(IX,KZ)*(FN(KZ) - FN(IX))
     > +SFF2BLK1(IX,KZ)*FN(KZ))*DIM + FG2BLK1(IX,KZ)*GI(KZ)*DIM**2*A) 

	GB21(NT1) = DXTZ(KZ)*(GG2BLK(IX,KZ)*(GI(KZ) - GI(IX))*DIM
     > + GF2BLK1(IX,KZ)*FN(KZ)/A)

	ENDIF

       IF (NPO.EQ.1) THEN

	IF (IX.LE.IV .AND. KZ.GT.IV) THEN
C
C     The relative signs, - for the quark kernel und + for the gluon switch for the
C     polarized case, as well as the signs infront of the QG and GQ and quark singlet 
C     piece.
C 

	IK = KZ - IV + 1

        IF (IKNL.EQ.2) THEN

	GB10(IK)=DXTZ(KZ)*(((PNSPBLK(IX,KZ)-PNSPBLKD(IX,KZ+1))*FN(KZ)
     >  + (SFF2BLK1(IX,KZ) - SFF2BLKD(IX,KZ+1))*FN(KZ))*DIM
     >  + (FG2BLK1(IX,KZ) - FG2BLKD(IX,KZ+1))*GI(KZ)*DIM**2*A)       

	GB210(IK) = DXTZ(KZ)*((GG2BLK(IX,KZ)+GG2BLKD(IX,KZ+1))*GI(KZ)
     > *DIM + (GF2BLK1(IX,KZ)+GF2BLKD(IX,KZ+1))*FN(KZ)/A) 

        ELSE

        GB10(IK)=DXTZ(KZ)*(((PNSPBLK(IX,KZ)+PNSPBLKD(IX,KZ+1))*FN(KZ)
     >  + (SFF2BLK1(IX,KZ) + SFF2BLKD(IX,KZ+1))*FN(KZ))*DIM
     >  + (FG2BLK1(IX,KZ) + FG2BLKD(IX,KZ+1))*GI(KZ)*DIM**2*A)
        
	GB210(IK) = DXTZ(KZ)*((GG2BLK(IX,KZ)-GG2BLKD(IX,KZ+1))*GI(KZ)
     > *DIM + (GF2BLK1(IX,KZ)-GF2BLKD(IX,KZ+1))*FN(KZ)/A) 


        ENDIF

        ELSEIF (IX.LE.IV.AND.KZ.EQ.IV) THEN

        IK = KZ - IV + 1

        IF (IKNL.EQ.2) THEN

	GB10(IK)=DIFFDEL*(((PNSPBLK(IX,KZ)-PNSPBLKD(IX,KZ+1))*FN(KZ)
     >  + (SFF2BLK1(IX,KZ) - SFF2BLKD(IX,KZ+1))*FN(KZ))*DIM
     >  + (FG2BLK1(IX,KZ) - FG2BLKD(IX,KZ+1))*GI(KZ)*DIM**2*A)
        
	GB210(IK) = DIFFDEL*((GG2BLK(IX,KZ)+GG2BLKD(IX,KZ+1))*GI(KZ)
     > *DIM + (GF2BLK1(IX,KZ)+GF2BLKD(IX,KZ+1))*FN(KZ)/A) 

        ELSE

        GB10(IK)=DIFFDEL*(((PNSPBLK(IX,KZ)+PNSPBLKD(IX,KZ+1))*FN(KZ)
     >  + (SFF2BLK1(IX,KZ) + SFF2BLKD(IX,KZ+1))*FN(KZ))*DIM
     >  + (FG2BLK1(IX,KZ) + FG2BLKD(IX,KZ+1))*GI(KZ)*DIM**2*A)

	GB210(IK) = DIFFDEL*((GG2BLK(IX,KZ)-GG2BLKD(IX,KZ+1))*GI(KZ)
     > *DIM + (GF2BLK1(IX,KZ)-GF2BLKD(IX,KZ+1))*FN(KZ)/A) 

        ENDIF
		
	ELSEIF (IX.GT.IV .AND. KZ.GT.IX .AND. IX.LT.NX-1) THEN

	XY = X/XV(KZ)

	IW1 = KZ - IV
        IW = KZ - IX + 1

        G1(IW) = DXTZ(KZ)*(PNSP(NR,IW1)*(FQ(KZ) - FQ(IX))
     > + SFF2(NR,IW1)*FQ(KZ)+ A*FG2(NR,IW1)*GI(KZ))
     
	G2(IW) = DXTZ(KZ)*(GG2(NR,IW1)*(GI(KZ) - XY*GI(IX))
     > +  GF2(NR,IW1)*FQ(KZ)/A)

	ENDIF

        ENDIF

	ENDIF

   31   CONTINUE

        GB10(NX-IV+1) = 0.0
	GB210(NX-IV+1) = 0.0

	IF (IX.LE.IV) THEN

C
C     Specify number of points in the various integration regions.
C

           NC = IX
           NC1 = IV-IX+1
           NC2 = NX-IV+1

C
C     Do the integration from 0..x, x..d, d..1 using Simpson
C

           IF (IX.LE.6) THEN

        TBLQQX = 0.0

        TBLGGX = 0.0
        
        ELSE

        TBLQQX = SMPNOR(NC,DZ,GBX,ERR)
	TBLGGX = SMPNOR(NC,DZ,GB2X,ERR) 

        ENDIF

        IF (IX.GE.IV-5) THEN

        TBLQQ1 = 0.0

        TBLGG1 = 0.0

        ELSE

      TBLQQ1 = SMPNOL(NC1,DZ,GB1,ERR)
      TBLGG1 = SMPNOL(NC1,DZ,GB21,ERR)


        ENDIF

        IF (NPO.EQ.1) THEN

      TBLQQD1 =SMPSNA(NC2,DZ,GB10,ERR)
      TBLGGD1 =SMPSNA(NC2,DZ,GB210,ERR)

        ELSE

       TBLQQD1 = 0D0
       TBLGGD1 = 0D0

        ENDIF

C
C	Save integration results in temp. array!
C

	TQQX(IX) = TBLQQX
	TGGX(IX) = TBLGGX
	TQQ1(IX) = TBLQQ1
	TGG1(IX) = TBLGG1
	TQQ11(IX) = TBLQQD1
	TGG11(IX) = TBLGGD1

	ELSE

      IF (NPO.EQ.1) THEN

C
C     Number of points in the integral in the DGLAP region
C

        NC3 = NX-IX
C
C     Value of the integrand at y=x. Only dummy value for book keeping purposes.
C

        G1(1) = 0.0
        G2(1) = 0.0
	G1(NX-IX) = -DXTZ(NX)*PNSP(NR,NX-IV)*FQ(IX)
	G2(NX-IX) = -DXTZ(NX)*GG2(NR,NX-IV)*GI(IX)*X 
       


C
C     Do integration in DGLAP region using Simpson
C

        IF (IX.LT.NX-3) THEN

	TEM = SMPNOL(NC3,DZ,G1,ERR)

	TEM1 = SMPNOL(NC3,DZ,G2,ERR) 

      
        ELSE

        TEM = 0.0

        TEM1 = 0.0

        ENDIF


C	Save results in temp. array!

	TQQD(IX) = TEM
	TGGD(IX) = TEM1

        ELSE

        TQQD(IX) = 0D0
        TGGD(IX) = 0D0

        ENDIF
   
	ENDIF

   21 CONTINUE
C

	IY = 6

	DO 9999 IM = 1,NX-1

C
C	Successive extrapolation of missing integration results for the 0..x integral 
C	(first) and the x..del integral (second). In both cases 6 bins. are missing. 
C	Extrapolation uses 8th. order polynomial extrapolation.
C

	IF (IM.LE.6) THEN
           
           KY = 8

	DO 998 IE = 1,KY
	X2(IE) = XV(IY+IE)
 998	CONTINUE
	XX = XV(IY)

	DO 997 IE = 1,KY
	FX(IE) = TQQX(IY+IE)
 997	CONTINUE

	CALL RATINT(X2,FX,KY,XX,TEM,ERR)

	TQQX(IY) = TEM


	DO 995 IE = 1,KY
	FX(IE) = TGGX(IY+IE)
 995	CONTINUE

	CALL RATINT(X2,FX,KY,XX,TEM,ERR)

	TGGX(IY) = TEM

	IY = IY - 1

	ELSEIF(IM.GE.IV-5 .AND.IM.LE.IV) THEN

           KU = 8

	DO 831 IE = 1,KU
	X2(IE) = XV(IM-IE)
 831	CONTINUE
	XX = XV(IM)

	DO 992 IE = 1,KU
	FX(IE) = TQQ1(IM-IE)
 992	CONTINUE

	CALL RATINT(X2,FX,KU,XX,TEM,ERR)

	TQQ1(IM) = TEM

	DO 991 IE = 1,KU
	FX(IE) = TGG1(IM-IE)
 991	CONTINUE

	CALL RATINT(X2,FX,KU,XX,TEM,ERR)

	TGG1(IM) = TEM

	ENDIF

C
C	Extrapolation of last 5 missing bins for DGLAP region, again using 9th. order 
C	polynomial approximation.
C
	IF (IM.GE.NX-3) THEN

           KT = 9

	DO 990 IH = 1,KT
	X2(IH) = XV(IM-IH)
 990	CONTINUE
	XX = XV(IM)

	DO 989 IH = 1,KT
	FX(IH) = TQQD(IM-IH)
 989	CONTINUE
	
	CALL RATINT(X2,FX,KT,XX,TEM,ERR)

	TQQD(IM) = TEM

	DO 987 IH = 1,KT
	FX(IH) = TGGD(IM-IH)
 987	CONTINUE
	
	CALL RATINT(X2,FX,KT,XX,TEM,ERR)

	TGGD(IM) = TEM

        ENDIF

 9999   CONTINUE


        DO 999 IM = 1,NX-1


C	2nd. order results. ERBL region

	IF (IM.LE.IV) THEN

C
C     first unpol. then pol.!
C

           IF (IKNL.EQ.2) THEN
	
	TMP = (TQQX(IM)+TQQ1(IM)+TQQ11(IM)
     >  - FN(IM)*(2./9.*(13./2. - PI2 + 2.*Z3)-ANSP(IM)))*CPL2

        TMX = ABS(XV(IM)-DEL/2.)
        IF (TMX.LT.1E-10) TMP = 0.0

	TMPO = (TGGX(IM)+TGG1(IM)+TGG11(IM) 
     >  - GI(IM)*(NEFF/108.*(105. + 296./3.) - AGG2(IM)))*CPL2

        ELSE

        TMP = (TQQX(IM)+TQQ1(IM)+TQQ11(IM)
     >  + FN(IM)*ANSP(IM)) * CPL2

	TMPO = (TGGX(IM)+TGG1(IM)+TGG11(IM)
     >  + GI(IM)*(9.*(95./27.-7.*PI2/9. 
     >  + Z3) + NEFF/54.*(87. - 112./3.) + AGG2(IM)))*CPL2

         TMX = ABS(XV(IM)-DEL/2.)
        IF (TMX.LT.1E-10) TMPO = 0.0

        ENDIF

C	1st. + 2nd. order

	FO(IM) = FO(IM) + TMP

	GO(IM) = GO(IM) + TMPO

	ELSE

        IF (NPO.EQ.1) THEN

C	2nd. order results DGLAP region.
C
C     first unpol. then pol.!
C

        IF (IKNL.EQ.2) THEN

	TMP1 = (TQQD(IM) -FQ(IM)*(ANSP(IM) +
     >  2./9.*(13./2. - PI2 + 2.*Z3)))*CPL2 

	TMP2 = (TGGD(IM)-GI(IM)*(AGG2(IM)+NEFF/108.*(105.
     >  + 296./3.)))* CPL2

        ELSE

        TMP1 = (TQQD(IM)-FQ(IM)*ANSP(IM))*CPL2 

	TMP2 = (TGGD(IM) + GI(IM)*(-AGG2(IM) + 9.*(95./27. 
     >  - 7.*PI2/9. + Z3) + NEFF/54.*(87.-112./3.)))*CPL2

        ENDIF

C	1st. + 2nd. order 

	FO(IM) = FO(IM) + TMP1*XA(IM,-1)

	GO(IM) = GO(IM) + TMP2

        ELSE

           FO(IM) = 0D0
           GO(IM) = 0D0

        ENDIF

	ENDIF

 999    CONTINUE

        IF (NPO.EQ.1) THEN

        KT = 10
        KT1 = 5

        IF (IKNL.EQ.2) THEN

           KT2 = 4
           KT21 = 1
           TIM = 0.1

           ELSE

           KT2 = 1
           KT21 = 1
           TIM = 0.01

        ENDIF

        DO 2011 IK = KT2+1,1,-1

        DO 1971 IE = 1,KT

        X2(IE) = 1D0*(IE+1)   
	FX(IE) = GO(IV-IE-IK+1)

 1971	CONTINUE


        XX = 1D0

        IF (DEL.LT.TIM) THEN

	CALL POLINT(X2,FX,KT1,XX,TEM,ERR1)

        ELSE

        KT1 = 8

        CALL RATINT(X2,FX,KT1,XX,TEM,ERR1)

        ENDIF
        
       	GO(IV-IK+1) = TEM

        IF (IKNL.EQ.2) THEN

           GO(IK) = GO(IV-IK+1)
           
           ELSE

           GO(IK) = -GO(IV-IK+1)
           
        ENDIF
         
 2011   CONTINUE

        DO 2012 IK = KT21+1,1,-1

        DO 1972 IE = 1,KT

        X2(IE) = 1D0*(IE+1)   
	
        FX1(IE) = FO(IV-IE-IK+1)

 1972   CONTINUE

        XX = 1D0

        IF (DEL.LT.TIM) THEN

        CALL RATINT(X2,FX1,KT,XX,TEM,ERR)

        ELSE

        KT = 5

        CALL POLINT(X2,FX1,KT,XX,TEM,ERR)

        ENDIF

       	FO(IV-IK+1) = TEM

        IF (IKNL.EQ.2) THEN

           FO(IK) = -FO(IV-IK+1)

           ELSE

           FO(IK) = FO(IV-IK+1) 

        ENDIF
         
 2012   CONTINUE

C
C     smooth quark in last two bins in ERBL region
C

        CALL SMOOTH1 (IKNL,MXX,DEL,IV,XV,FO)

        ELSE

        FO(1) = 0D0
        GO(1) = 0D0
        GO(IV) = 0D0
        FO(IV) = 0D0

        ENDIF

        IF (DEL.GT.0.01) THEN

           KT = 5

           DO 1973 IP = 2,1,-1

        DO 1977 IE = 1,KT
        X2(IE) = 1D0*(IE+1)
	FX(IE) = GO(IV+IE+IP)
        FX1(IE) = FO(IV+IE+IP)
 1977   CONTINUE

        XX = 1D0

	CALL POLINT(X2,FX,KT,XX,TEM,ERR1)

       	GO(IV+IP) = TEM

        CALL POLINT(X2,FX1,KT,XX,TEM,ERR)

       	FO(IV+IP) = TEM
           
 1973   CONTINUE

        ENDIF
       
      ENDIF

      RETURN
C                         ************************* 
      END

C	---------------------------------------------- 
 
      FUNCTION FF1F (XX,YY) 
 
	IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
 
 
C	FIRST LO DGLAP KERNELS AS GIVEN IN hep-ph/9912379. 
C	------------- 
 
	X = XX 
 
	D = YY 
 
	XM1 = 1.- X 
 
	DM1 = 1. - D 
 
	X2 = X**2 
 
	FF1F = (1. + X2 - D*(1.+ X))/(XM1*DM1) 
 
	RETURN 
 
C	------------- 
 
	ENTRY FG1F (XX,YY) 
 
	X = XX 
 
	D = YY 
 
	XM1 = 1.- X 
 
	DM1 = 1. - D 
 
	X2 = X**2 
 
	FG1F = 2.*(X2 + XM1**2 - X*D)/DM1**2 
 
	RETURN 
 
C	------------- 
 
	ENTRY FG1FA (XX,YY) 
 
	X = XX 
 
	D = YY 
 
	XM1 = 1.- X 
 
	DM1 = 1. - D 
 
	FG1FA = 2.*(2.*X - 1. - X*D)/DM1**2 
 
	RETURN	 
 
C	------------- 
 
	ENTRY GF1F (XX,YY) 
 
	X = XX 
 
	D = YY 
 
	XM1 = 1.- X 
 
	DM1 = 1. - D 
 
	GF1F = (1. + XM1**2 - D)/DM1 
 
	RETURN 
 
C	-------------- 
 
	ENTRY GF1FA (XX,YY) 
 
	X = XX 
 
	D = YY 
 
	DM1 = 1. - D 
 
	GF1FA = (X*(2.- X) - D)/DM1 
 
	RETURN 
 
C	-------------- 
 
	ENTRY GG1F (XX,YY) 
 
	X = XX 
 
	D = YY 
 
	XM1 = 1.- X 
 
	DM1 = 1. - D 
 
	X2 = X**2 
 
	XMDM1 = 1. - XM1/DM1 
 
	GG1F = (X2 + XMDM1)/XM1 + 2.*(1.- XMDM1)**2  
     >  + 2.*(1./2. - X2)*XMDM1/DM1 
 
	RETURN 
 
C       ------------ 
 
	ENTRY GG1FA (XX,YY) 
 
	X = XX 
 
	D = YY 
 
	XM1 = 1.- X 
 
	DM1 = 1. - D 
 
	X2 = X**2 
 
	GG1FA = (D**2*(1.+ X2) + 2.*X*(X*DM1 - D)) 
     >  /(DM1**2*XM1) 
     >  + (4.*X*XM1 - 2.*D*(1.- X2))/DM1**2 
 
	RETURN 
 
C	--------------- 
 
	ENTRY FGC1F (XX,YY) 
 
	X = XX 
 
	D = YY 
 
	XM1 = 1.- X 
 
	DM1 = 1. - D 
 
	XMDM1 = 1. - XM1/DM1 
 
	FGC1F = 2.*(1.- XMDM1)**2 
 
	RETURN 
 
C	----------------- 
 
	ENTRY GGA1F (XX,YY) 
 
	X = XX 
 
	D = YY 
 
	XM1 = 1.- X 
 
	DM1 = 1. - D 
 
	X2 = X**2 
 
	GGA1F = (2.*X*XM1 - D*(1.- X2))/DM1**2 
 
	RETURN 
 
C	------------------ 
 
	ENTRY GGC1F (XX,YY) 
 
	X = XX 
 
	D = YY 
 
	XM1 = 1.- X 
 
	DM1 = 1. - D 
 
	XMDM = XM1/DM1 
 
	XMDM2 = XMDM**2  
 
	GGC1F = XM1*XMDM2 
 
	RETURN 
 
C	------------------ 
C	NOW THE BL KERNELS 
 
	ENTRY FF1BLF (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XYM1 = 1. - XY  
 
	FF1BLF = XY*(1. + 1./(Y*XYM1)) 
 
	RETURN 
 
C	----------------- 
 
	ENTRY FG1BLF (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XM1 = 1.- X 
 
	Y2 = Y**2 
 
	FG1BLF = - X*(1. + 2.*(2.*XM1*Y - X))/Y2 
 
	RETURN 
 
C	----------- 
 
	ENTRY FG1BLFA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	Y2 = Y**2 
 
	FG1BLFA = - X/Y2 
 
	RETURN 
 
C	------------- 
 
	ENTRY GF1BLF (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	X2 = X**2 
 
	GF1BLF = X2*(1. + 2.*(2.*XM1*Y - YM1))/Y 
 
	RETURN 
 
C	-------------- 
 
	ENTRY GF1BLFA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	X2 = X**2 
 
	GF1BLFA = X2/Y 
 
	RETURN 
 
C	---------------- 
 
	ENTRY GG1BLF (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY  
 
	GG1BLF = XY2*(2.+ 1./(Y*XYM1) + 2.*(2.*XM1*Y + Y - X))  
 
	RETURN 
 
C	--------------------- 
 
	ENTRY GG1BLFA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY  
 
	GG1BLFA = XY2*(2.+ 1./(Y*XYM1)) 
 
	RETURN 
 
C	------------- 
 
	ENTRY GG1BLCF (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XM1 = 1.- X 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	GG1BLCF = XY2*(2.*XM1*Y + Y - X) 
 
	RETURN  
 
C	----------------	 
 
	END 
 
C	------------------------- 
 
	FUNCTION HFFBL (XX,YY) 
 
	IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
 
	PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2) 
        PARAMETER (D0 = 0.0, D1 = 1.0, TINY = 1E-10) 
 
	EXTERNAL FF1F, FF1BLF, FG1F, FG1FA, FG1BLF, FG1BLFA, GF1F, GF1FA,  
     >  GF1BLF, GF1BLFA, GG1F, GG1FA, GG1BLF, GG1BLFA, GGA1F 
C	Setup of the h functions for BL and DGLAP kernels, First BL, then DGLAP 
 
	X = XX 
 
	Y = YY 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XLNY = LOG(Y) 
 
	XLNXM1 = LOG(XM1)	 
 
	SPENX = SPENCE(X) 
 
	SPENYM1 = SPENCE(YM1) 
 
	FF1 = FF1BLF(X,Y) 
 
	FF1B = FF1BLF(XM1,YM1) 

	HFFBL = 2.*FF1B*XLNXM1*XLNY-2.*FF1*(SPENX+SPENYM1) 
 
	RETURN 
 
C	------------------- 
 
	ENTRY HFFBLB (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNY2 = XLNY**2 
		  
	SPENXM1 = SPENCE(XM1) 
 
	SPENYM1 = SPENCE(YM1) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	FF1 = FF1BLF(X,Y) 
 
	FF1B = FF1BLF(XM1,YM1) 

        IF (ABS(X-Y).LE.TINY) THEN

        HFFBLB = 2.*(2.*SPENYM1 + 2./Y - XLNY2 + (XLNY*(1. - 
     >  XLNY) - SPENYM1 + 3./2.*XLNY2)/(Y*YM1)) 

        ELSE

	HFFBLB = (FF1 - FF1B)*(2.*SPENXYM1 + XLNY2) + 2.*FF1* 
     >  (SPENYM1 - XLNX*XLNY) + 2.*FF1B*SPENXM1 

        ENDIF

	RETURN 

C	--------------------- 
 
	ENTRY HFGBL (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XLNY = LOG(Y) 
 
	XLNXM1 = LOG(XM1)	 
 
	SPENX = SPENCE(X) 
 
	SPENYM1 = SPENCE(YM1) 
 
	FG1 = FG1BLF(X,Y) 
 
	FG1B = FG1BLF(XM1,YM1) 
 
	HFGBL = -2.*FG1B*XLNXM1*XLNY-2.*FG1*(SPENX+SPENYM1) 

	RETURN 
 
C	--------------------- 
 
	ENTRY HFGBLB (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNY2 = XLNY**2 
		  
	SPENXM1 = SPENCE(XM1) 
 
	SPENYM1 = SPENCE(YM1) 
 
	SPENXYM1 = SPENCE(XYM1)	 
 
	FG1 = FG1BLF(X,Y) 
 
	FG1B = FG1BLF(XM1,YM1) 
 
	HFGBLB = (FG1 + FG1B)*(2.*SPENXYM1 + XLNY2) + 2.*FG1* 
     >  (SPENYM1 - XLNX*XLNY) - 2.*FG1B*SPENXM1 

	RETURN 
 
C	--------------------- 
 
	ENTRY HGFBL (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XLNY = LOG(Y) 
 
	XLNXM1 = LOG(XM1)	 
 
	SPENX = SPENCE(X) 
 
	SPENYM1 = SPENCE(YM1) 
 
	GF1 = GF1BLF(X,Y) 
 
	GF1B = GF1BLF(XM1,YM1) 
 
	HGFBL = -2.*GF1B*XLNXM1*XLNY-2.*GF1*(SPENX+SPENYM1) 
 
	RETURN 
 
C	--------------------- 
 
	ENTRY HGFBLB (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNY2 = XLNY**2 
		 
	SPENXM1 = SPENCE(XM1) 
 
	SPENYM1 = SPENCE(YM1) 
 
	SPENXYM1 = SPENCE(XYM1)	 
 
	GF1 = GF1BLF(X,Y) 
 
	GF1B = GF1BLF(XM1,YM1) 
 
	HGFBLB = (GF1 + GF1B)*(2.*SPENXYM1 + XLNY2) + 2.*GF1* 
     >  (SPENYM1 - XLNX*XLNY) - 2.*GF1B*SPENXM1 

	RETURN 
 
C	--------------------- 
 
	ENTRY HGGBL (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XLNY = LOG(Y) 
 
	XLNXM1 = LOG(XM1)	 
 
	SPENX = SPENCE(X) 
 
	SPENYM1 = SPENCE(YM1) 
 
	GG1 = GG1BLF(X,Y) 
 
	GG1B = GG1BLF(XM1,YM1) 
 
	HGGBL = 2.*GG1B*XLNXM1*XLNY-2.*GG1*(SPENX+SPENYM1) 
 
	RETURN 
 
C	------------------- 
 
	ENTRY HGGBLB (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNY2 = XLNY**2 
		  
	SPENXM1 = SPENCE(XM1) 
 
	SPENYM1 = SPENCE(YM1) 
 
	SPENXYM1 = SPENCE(XYM1)	 
	 
	GG1 = GG1BLF(X,Y) 
 
	GG1B = GG1BLF(XM1,YM1)

        IF (ABS(X-Y).LE.TINY) THEN

        HGGBLB = 2.*((4.*SPENYM1 - 2.*XLNY2)*(1. + 2.*Y*YM1) + 2./Y
     >  + (XLNY*(1. -  XLNY) - 2.*SPENYM1 + 2.*XLNY2)/(Y*YM1))

        ELSE

	HGGBLB = (GG1 - GG1B)*(2.*SPENXYM1 + XLNY2) + 2.*GG1* 
     >  (SPENYM1 - XLNX*XLNY) + 2.*GG1B*SPENXM1 

        ENDIF

	RETURN 
 
C	--------------------- 
 
	ENTRY HFGBLA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XLNY = LOG(Y)	 
 
	XLNXM1 = LOG(XM1) 
 
	SPENX = SPENCE(X) 
 
	SPENYM1 = SPENCE(YM1) 
 
	FG1 = FG1BLFA(X,Y) 
 
	FG1B = FG1BLFA(XM1,YM1) 
 
	HFGBLA = -2.*FG1B*XLNXM1*XLNY-2.*FG1*(SPENX+SPENYM1) 
 
	RETURN 
 
C	--------------------- 
 
	ENTRY HFGBLBA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNY2 = XLNY**2 
		  
	SPENXM1 = SPENCE(XM1) 
 
	SPENYM1 = SPENCE(YM1) 
 
	SPENXYM1 = SPENCE(XYM1)	 
 
	FG1 = FG1BLFA(X,Y) 
 
	FG1B = FG1BLFA(XM1,YM1) 
 
	HFGBLBA = (FG1 + FG1B)*(2.*SPENXYM1 + XLNY2) + 2.*FG1* 
     >  (SPENYM1 - XLNX*XLNY) - 2.*FG1B*SPENXM1 
 
	RETURN 
 
C	--------------------- 
 
	ENTRY HGFBLA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XLNY = LOG(Y) 
 
	XLNXM1 = LOG(XM1)	 
 
	SPENX = SPENCE(X) 
 
	SPENYM1 = SPENCE(YM1) 
 
	GF1 = GF1BLFA(X,Y) 
 
	GF1B = GF1BLFA(XM1,YM1) 
 
	HGFBLA = -2.*GF1B*XLNXM1*XLNY-2.*GF1*(SPENX+SPENYM1) 
 
	RETURN 
 
C	--------------------- 
 
	ENTRY HGFBLBA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNY2 = XLNY**2 
		  
	SPENXM1 = SPENCE(XM1) 
 
	SPENYM1 = SPENCE(YM1) 
 
	SPENXYM1 = SPENCE(XYM1)	 
 
	GF1 = GF1BLFA(X,Y) 
 
	GF1B = GF1BLFA(XM1,YM1) 
 
	HGFBLBA = (GF1 + GF1B)*(2.*SPENXYM1 + XLNY2) + 2.*GF1* 
     >  (SPENYM1 - XLNX*XLNY) - 2.*GF1B*SPENXM1 
 
	RETURN 
 
C	--------------------- 
 
	ENTRY HGGBLA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XLNY = LOG(Y) 
 
	XLNXM1 = LOG(XM1)	 
 
	SPENX = SPENCE(X) 
 
	SPENYM1 = SPENCE(YM1) 
 
	GG1 = GG1BLFA(X,Y) 
 
	GG1B = GG1BLFA(XM1,YM1) 
 
	HGGBLA = 2.*GG1B*XLNXM1*XLNY-2.*GG1*(SPENX+SPENYM1) 
 
	RETURN 
 
C	------------------- 
 
	ENTRY HGGBLBA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNY2 = XLNY**2 
		  
	SPENXM1 = SPENCE(XM1) 
 
	SPENYM1 = SPENCE(YM1) 
 
	SPENXYM1 = SPENCE(XYM1)	 
 
	GG1 = GG1BLFA(X,Y) 
 
	GG1B = GG1BLFA(XM1,YM1)

        IF (ABS(X-Y).LE.TINY) THEN

        HGGBLBA = 2.*(4.*SPENYM1 - 2.*XLNY2 + 2./Y
     >  + (XLNY*(1. -  XLNY) - 2.*SPENYM1 + 2.*XLNY2)/(Y*YM1))

        ELSE

	HGGBLBA = (GG1 - GG1B)*(2.*SPENXYM1 + XLNY2) + 2.*GG1* 
     >  (SPENYM1 - XLNX*XLNY) + 2.*GG1B*SPENXM1 

        ENDIF

	RETURN 
 
C	---------------------------- 
C	Now DGLAP 
C	---------------------------- 
 
	ENTRY HFFD (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNXY = LOG(XY) 
 
	XLNXMYM1 = LOG(XMYM1) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	FF1 = FF1F(X,Y) 
 
	XK1 = XY + X/XM1 
 
	HFFD = 2.*FF1*(XLNY*XLNXMYM1 + SPEN1MY - SPENXYM1 - PI2/6.) 
     >  + 2.*XK1*(XLNXMYM1*(XLNXY - XLNY) -2.*(SPEN1MY - SPENXYM1)) 

	RETURN 
 
C	-------------------------- 
 
	ENTRY HFFDB (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXY = LOG(XY) 
 
	XLNXMYM1 = LOG(XMYM1) 
 
	XLNYM1Y = LOG(YM1/Y) 
 
	XLNYM1X = LOG(YM1+X) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	SPENYMX = SPENCE(Y-X) 
 
	SPENMXBY = SPENCE(-X/YM1) 
 
	FF1 = FF1F(Y-X,Y) 
 
	XK2 = XY*(1. - Y/(1.-Y+X))/YM1 
 
	HFFDB = 2.*FF1*(XLNY*XLNXMYM1 + 1./2.*XLNYM1*(XLNXY + XLNX 
     >  - XLNYM1Y) + SPEN1MY - SPENXYM1 - PI2/6.  
     >  - XLNYM1X*(XLNXMYM1 + XLNX) - SPENYMX - SPENMXBY) 
     >  + 2.*XK2*(XLNXMYM1*(XLNXY - XLNY) -2.*(SPEN1MY - SPENXYM1)) 

	RETURN 
 
C	-------------------------- 
 
	ENTRY HFGD (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNXY = LOG(XY) 
 
	XLNXMYM1 = LOG(XMYM1) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	FG1 = FG1F(X,Y) 
 
	XK1 = -2.*(X - 2.*XY2*Y + 4.*XYM1*XY)/Y 
 
	HFGD = 2.*FG1*(XLNY*XLNXMYM1 + SPEN1MY - SPENXYM1 - PI2/6.) 
     >  + 2.*XK1*(XLNXMYM1*(XLNXY - XLNY) -2.*(SPEN1MY - SPENXYM1)) 

	RETURN 
 
C	-------------------------- 
 
	ENTRY HFGDB (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXY = LOG(XY) 
 
	XLNXMYM1 = LOG(XMYM1) 
 
	XLNYM1Y = LOG(YM1/Y) 
 
	XLNYM1X = LOG(YM1+X) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	SPENYMX = SPENCE(Y-X) 
 
	SPENMXBY = SPENCE(-X/YM1) 
 
	FG1 = FG1F(Y-X,Y) 
 
	XK2 = 2.*XY*(1. + XYM1*(4. - 6.*Y)/Y)/YM1**2 
 
	HFGDB = 2.*FG1*(XLNY*XLNXMYM1 + 1./2.*XLNYM1*(XLNXY + XLNX 
     >  - XLNYM1Y) + SPEN1MY - SPENXYM1 - PI2/6.  
     >  - XLNYM1X*(XLNXMYM1 + XLNX) - SPENYMX - SPENMXBY) 
     >  + 2.*XK2*(XLNXMYM1*(XLNXY - XLNY) -2.*(SPEN1MY - SPENXYM1)) 

	RETURN 
 
C	-------------------------- 
 
	ENTRY HGFD (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	XYM1 = 1. - XY 
 
	XLNY = -LOG(1./Y) 
 
	XLNXY = LOG(XY/Y) 
 
	XLNXMYM1 = LOG(XMYM1) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	GF1 = GF1F(X,Y) 
 
	XK1 = - X*XY*(1.- (6.- 4.*XY)/Y) 
  
	HGFD = 2.*GF1*(XLNY*XLNXMYM1 + SPEN1MY - SPENXYM1 - PI2/6.) 
     >  + 2.*XK1*(XLNXMYM1*(XLNXY) -2.*(SPEN1MY - SPENXYM1)) 

	RETURN 
 
C	-------------------------- 
 
	ENTRY HGFDB (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	XYM1 = 1. - XY 
 
	RE = X**2/YM1 
 
	XLNY = -LOG(1./Y) 
 
	XLNX = -LOG(1./X) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXY = LOG(XY/Y) 
 
	XLNXMYM1 = LOG(XMYM1) 
 
	XLNYM1Y = LOG(RE) 
 
	XLNYM1X = LOG(YM1+X) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	SPENYMX = SPENCE(Y-X) 
 
	SPENMXBY = SPENCE(-X/YM1) 
 
	GF1 = GF1F(Y-X,Y) 
 
	XK2 = X*XY*(1. + YM1*(6. - 4.*XY)/Y)/YM1 
 
	HGFDB = 2.*GF1*(XLNY*XLNXMYM1 + 1./2.*XLNYM1*( 
     >  XLNYM1Y) + SPEN1MY - SPENXYM1 - PI2/6.  
     >  - XLNYM1X*(XLNXMYM1 + XLNX) - SPENYMX - SPENMXBY) 
     >  + 2.*XK2*(XLNXMYM1*(XLNXY) -2.*(SPEN1MY - SPENXYM1)) 

	RETURN 
 
C	-------------------------- 
 
	ENTRY HGGD (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY 
 
	XLNY = -LOG(1./Y) 
 
	XLNXY = LOG(XY) 
 
	XLNXMYM1 = -LOG(1./XMYM1) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	GG1 = GG1F(X,Y) 
 
	XK1 = (X2/XM1 + 2.*XY2 + 2.*X*XYM1*XY*(1.+ 2./Y)) 
  
	TEM1 = XLNY*XLNXMYM1 + SPEN1MY - SPENXYM1 - PI2/6. 
	HGGD1 = 2.*GG1*TEM1 
 
	TEM2 = XLNXMYM1*XLNXY - 2.*SPEN1MY 
     >  - XLNXMYM1*XLNY + 2.*SPENXYM1 
	HGGD2 = 2.*XK1*TEM2 
 
	HGGD = HGGD1 + HGGD2 

	RETURN 
 
C	-------------------------- 
 
	ENTRY HGGDB (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	XYM1 = 1. - XY 
 
	XLNY = -LOG(1./Y) 
 
	XLNX = -LOG(1./X) 
 
	XLNXM1 = -LOG(1./XM1) 
 
	XLNYM1 = -LOG(1./YM1) 
 
	XLNXY = LOG(XY) 
 
	XLNXMYM1 = -LOG(1./XMYM1) 
 
	XLNYM1Y = -LOG(Y/YM1) 
 
	XLNYM1X = -LOG(1./(YM1+X)) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	SPENYMX = SPENCE(Y-X) 
 
	SPENMXBY = SPENCE(-X/YM1) 
 
	GG1 = GG1F(Y-X,Y) 
 
	XK2 = (XY/YM1)**2*(10.+ 6.*X - 9.*Y +  
     >  (Y*(1.+ X)*(1.- 4.*XY/Y) - 4.*YM1)/(YM1+X)) 
 
	TEM1 = XLNY*XLNXMYM1 + 1./2.*XLNYM1*XLNXY 
     >  + 1./2.*XLNYM1*XLNX - 1./2.*XLNYM1*XLNYM1Y  
     >  - XLNYM1X*XLNXMYM1 - XLNYM1X*XLNX  
 
	TEM11 = SPEN1MY - SPENXYM1 - PI2/6. - SPENYMX - SPENMXBY 
 
	HGGDB1 = 2.*GG1*(TEM1+TEM11) 
 
	TEM2 =  XLNXMYM1*(XLNXY - XLNY)  
        TEM22 = - 2.*(SPEN1MY - SPENXYM1) 
 
	HGGDB2 = 2.*XK2*TEM2+2.*XK2*TEM22 
 
	HGGDB = HGGDB1 + HGGDB2 

	RETURN 
 
C	-------------------------- 
 
	ENTRY HFGDA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNXY = LOG(XY) 
 
	XLNXMYM1 = LOG(XMYM1) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	FG1 = FG1FA(X,Y) 
 
	XK1 = -2.*XY  
 
	HFGDA = 2.*FG1*(XLNY*XLNXMYM1 + SPEN1MY - SPENXYM1 - PI2/6.) 
     >  + 2.*XK1*(XLNXMYM1*(XLNXY - XLNY) -2.*(SPEN1MY - SPENXYM1)) 
 
	RETURN 
 
C	-------------------------- 
 
	ENTRY HFGDBA (XX,YY) 
 
        X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	XYM1 = 1. - XY 
 
	RE = X**2/YM1 
 
	XLNY = -LOG(1./Y) 
 
	XLNX = -LOG(1./X) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXY = LOG(XY/Y) 
 
	XLNXMYM1 = LOG(XMYM1) 
 
	XLNYM1Y = LOG(RE) 
 
	XLNYM1X = LOG(YM1+X) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	SPENYMX = SPENCE(Y-X) 
 
	SPENMXBY = SPENCE(-X/YM1) 

	FG1 = FG1FA(Y-X,Y) 
 
	XK2 = - 2.*XY/YM1**2 
 
	HFGDBA = 2.*FG1*(XLNY*XLNXMYM1 + 1./2.*XLNYM1*( 
     >  XLNYM1Y) + SPEN1MY - SPENXYM1 - PI2/6.  
     >  - XLNYM1X*(XLNXMYM1 + XLNX) - SPENYMX - SPENMXBY) 
     >  + 2.*XK2*(XLNXMYM1*(XLNXY) -2.*(SPEN1MY - SPENXYM1)) 
 
	RETURN 
 
C	-------------------------- 
 
	ENTRY HGFDA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNXY = LOG(XY) 
 
	XLNXMYM1 = LOG(XMYM1) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	GF1 = GF1FA(X,Y) 
 
	XK1 = X*XY 
  
	HGFDA = 2.*GF1*(XLNY*XLNXMYM1 + SPEN1MY - SPENXYM1 - PI2/6.) 
     >  + 2.*XK1*(XLNXMYM1*(XLNXY - XLNY) -2.*(SPEN1MY - SPENXYM1)) 
 
	RETURN 
 
C	-------------------------- 
 
	ENTRY HGFDBA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXY = LOG(XY) 
 
	XLNXMYM1 = LOG(XMYM1) 
 
	XLNYM1Y = LOG(YM1/Y) 
 
	XLNYM1X = LOG(YM1+X) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	SPENYMX = SPENCE(Y-X) 
 
	SPENMXBY = SPENCE(-X/YM1) 
 
	GF1 = GF1FA(Y-X,Y) 
 
	XK2 = - X*XY/YM1 
 
	HGFDBA = 2.*GF1*(XLNY*XLNXMYM1 + 1./2.*XLNYM1*(XLNXY + XLNX 
     >  - XLNYM1Y) + SPEN1MY - SPENXYM1 - PI2/6.  
     >  - XLNYM1X*(XLNXMYM1 + XLNX) - SPENYMX - SPENMXBY) 
     >  + 2.*XK2*(XLNXMYM1*(XLNXY - XLNY) -2.*(SPEN1MY - SPENXYM1)) 
 
	RETURN 
 
C	-------------------------- 
 
	ENTRY HGGDA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	X2 = X**2 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNXY = LOG(XY) 
 
	XLNXMYM1 = LOG(XMYM1) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	GG1 = GG1FA(X,Y) 
 
	XK1 = X2/XM1 + 2.*X*XY 
  
	HGGDA = 2.*GG1*(XLNY*XLNXMYM1 + SPEN1MY - SPENXYM1 - PI2/6.) 
     >  + 2.*XK1*(XLNXMYM1*(XLNXY - XLNY) -2.*(SPEN1MY - SPENXYM1)) 
 
	RETURN 
 
C	-------------------------- 
 
	ENTRY HGGDBA (XX,YY) 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XMYM1 = 1. - XM1/YM1 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXY = LOG(XY) 
 
	XLNXMYM1 = LOG(XMYM1) 
 
	XLNYM1Y = LOG(YM1/Y) 
 
	XLNYM1X = LOG(YM1+X) 
 
	SPEN1MY = SPENCE(1.- 1./Y) 
 
	SPENXYM1 = SPENCE(XYM1) 
 
	SPENYMX = SPENCE(Y-X) 
 
	SPENMXBY = SPENCE(-X/YM1) 
 
	GG1 = GG1FA(Y-X,Y) 
 
	XK2 = -Y*(XY/YM1)**2*(2.- Y/(YM1+X)) 
 
	HGGDBA = 2.*GG1*(XLNY*XLNXMYM1 + 1./2.*XLNYM1*(XLNXY + XLNX 
     >  - XLNYM1Y) + SPEN1MY - SPENXYM1 - PI2/6.  
     >  - XLNYM1X*(XLNXMYM1 + XLNX) - SPENYMX - SPENMXBY) 
     >  + 2.*XK2*(XLNXMYM1*(XLNXY - XLNY) -2.*(SPEN1MY - SPENXYM1)) 
 
	RETURN 
 
C	-------------------------- 
 
	END 
 
C	-------------------------- 
 
	FUNCTION PNSMBL (N,XX,YY) 
 
	IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
 
	PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2) 
        PARAMETER (D0 = 0.0, D1 = 1.0) 
 
	DATA CF, CG, TF / 1.333333333333, 3.0, 0.5/ 
 
	EXTERNAL FF1F, FF1BLF, FG1F, FG1FA, FG1BLF, FG1BLFA, GF1F,   
     >  GF1BLF, GF1BLFA, GG1F, GG1FA, GG1BLF, GG1BLFA, GGA1F, FGC1F,  
     >  GG1BLCF, HFFBL, HFFBLB, HFGBL, HFGBLB, HFGBLA, HFGBLBA, HGFBL,  
     >  HGFBLB, HGFBLA, HGFBLBA, HGGBL, HGGBLB, HGGBLA, HGGBLBA, HFFD,   
     >  HFFDB, HFGD, HFGDB, HFGDA, HFGDBA, HGFD, HGFDB, HGFDA, HGFDBA,   
     >  HGGD, HGGDB, HGGDA, HGGDBA, GF1FA, GGC1F 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
C	The NLO kernels, first for the ERBL case. As in hep-ph/9912379. 
C	Also remember that for the polarized case PNSMBLA = PNSPBL and  
C	PNSPBLA = PNSMBL!  
C	------------------------- 
 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	YX = Y/X 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XOYM1 = X/YM1 
 
	XM1OY = XM1/Y 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXY = LOG(XY) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXYM1 = LOG(XYM1) 
 
	XLNYXM1 = LOG(YX - 1.) 
 
	XLNY2 = XLNY**2 
 
	XLNX2 = XLNX**2 
 
	XLNXY2 = XLNXY**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNYM12 = XLNYM1**2 
 
	XLNXYM12 = XLNXYM1**2 
 
	XLNYXM12 = XLNYXM1**2 
 
	FF1 = FF1BLF(X,Y) 
 
	FF1B = FF1BLF(XM1,YM1) 
 
	HFF = HFFBL(X,Y) 
 
	HFFB = HFFBLB(XM1,Y) 
 
	CF2 = (4./3. - PI2/3.)*FF1 + 3.*XY - (3./2.*FF1 - XOYM1/2.) 
     >  *XLNXY - (FF1 - FF1B)*XLNXY*XLNXYM1 + (FF1 + XOYM1/2.)*XLNXY2 
     >  - XOYM1/2.*XLNX*(1. + XLNX - 2.*XLNXM1) - XM1OY/2.*XLNXM1* 
     >  (1. + XLNXM1 - 2.*XLNX) 
 
	CFBO = 5./3.*FF1 + XY + FF1*XLNXY 
 
	CFCFMCA = 4./3.*FF1 + 2.*XY + HFF + HFFB 
 
	PNSMBL = CF**2*CF2 - CF*BO*CFBO/2. - CF*(CF-CG/2.)*CFCFMCA 

	RETURN 
 
C	-------------------- 
 
	ENTRY PNSPBL (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	YX = Y/X 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XOYM1 = X/YM1 
 
	XM1OY = XM1/Y 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXY = LOG(XY) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXYM1 = LOG(XYM1) 
 
	XLNYXM1 = LOG(YX - 1.) 
 
	XLNY2 = XLNY**2 
 
	XLNX2 = XLNX**2 
 
	XLNXY2 = XLNXY**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNYM12 = XLNYM1**2 
 
	XLNXYM12 = XLNXYM1**2 
 
	XLNYXM12 = XLNYXM1**2 
 
	FF1 = FF1BLF(X,Y) 
 
	FF1B = FF1BLF(XM1,YM1) 
 
	HFF = HFFBL(X,Y) 
 
	HFFB = HFFBLB(XM1,Y) 
 
	CF2 = (4./3. - PI2/3.)*FF1 + 3.*XY - (3./2.*FF1 - XOYM1/2.) 
     >  *XLNXY - (FF1 - FF1B)*XLNXY*XLNXYM1 + (FF1 + XOYM1/2.)*XLNXY2 
     >  - XOYM1/2.*XLNX*(1. + XLNX - 2.*XLNXM1) - XM1OY/2.*XLNXM1* 
     >  (1. + XLNXM1 - 2.*XLNX) 
 
	CFBO = 5./3.*FF1 + XY + FF1*XLNXY 
 
	CFCFMCA = 4./3.*FF1 + 2.*XY + HFF - HFFB 
 
	PNSPBL = CF**2*CF2 - CF*BO*CFBO/2. - CF*(CF-CG/2.)*CFCFMCA 

	RETURN 
 
C	-------------------- 
 
	ENTRY SFF2BLA (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	YX = Y/X 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XOYM1 = X/YM1 
 
	XM1OY = XM1/Y 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXY = LOG(XY) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXYM1 = LOG(XYM1) 
 
	XLNYXM1 = LOG(YX - 1.) 
 
	CFNFTF = XOYM1*XLNXY*(1. - XLNXY) 
     >  - XOYM1*XLNX*(1.- XLNX) 
     >  - XM1OY*XLNXM1*(1.- XLNXM1) 
 
	SFF2BLA = - 2.*CF*TFNF*CFNFTF 
 
	RETURN 
 
C	---------------------------- 
 
	ENTRY FG2BLA (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	YX = Y/X 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XOYM1 = X/YM1 
 
	XM1OY = XM1/Y 
 
	XMYM = XM1/YM1 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXY = LOG(XY) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXYM1 = LOG(XYM1) 
 
	XLNYXM1 = LOG(YX - 1.) 
 
	XLNY2 = XLNY**2 
 
	XLNX2 = XLNX**2 
 
	XLNXY2 = XLNXY**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNYM12 = XLNYM1**2 
 
	XLNXYM12 = XLNXYM1**2 
 
	XLNYXM12 = XLNYXM1**2 
 
	FG1 = FG1BLFA(X,Y) 
 
	FG1B = FG1BLFA(XM1,YM1) 
 
	HFG = HFGBLA(X,Y) 
 
	HFGB = HFGBLBA(XM1,Y) 
 
	CFNFTF = 2.*(5. - PI2/3.)*FG1 + 5.*XY/YM1 + (2.*FG1  
     >  - 3.*XOYM1/YM1 + 2.*XY/YM1)*XLNXY - 2.*(FG1 + FG1B  
     >  + 1./(Y*YM1))*XLNXYM1  - FG1B*XLNXY2 
     >  + (FG1 + FG1B)*XLNYXM12 - FG1*XLNX*(2. + 2.*XLNXM1 - XLNX)  
     >  + 3.*XOYM1*XLNX/YM1 + FG1B*XLNXM1*(2. + 2.*XLNX - XLNXM1)  
     >  - 3.*XM1OY*XLNXM1/Y 
 
	CGNFTF = -4.*FG1 - 6.*XY/YM1 - 2.*(FG1 - 2.*XOYM1/YM1)*XLNXY  
     >  + 2.*(FG1 + FG1B + 1./(Y*YM1))*XLNXYM1 + (FG1 - 3.*XOYM1/YM1) 
     >  *XLNXY2 - (FG1 + FG1B)*XLNXYM12 - HFG - HFGB 
     >   + (6.*XY/YM1 - 4.*XOYM1/YM1)*XLNX + FG1*XLNX*(2. - 
     >  XLNX) + 3.*XOYM1*XLNX2/YM1 - (6.*XMYM/Y - 4.*XM1OY/Y)*XLNXM1  
     >  - FG1B*XLNXM1*(2. - XLNXM1) - 3.*XM1OY*XLNXM12/Y 
 
	FG2BLA = CF*TFNF*CFNFTF + CG*TFNF*CGNFTF 
 
        RETURN
 
C	------------------------ 
 
	ENTRY GF2BLA(N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	YX = Y/X 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XOYM1 = X/YM1 
 
	XM1OY = XM1/Y 
 
	XMYM = XM1/YM1 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXY = LOG(XY) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXYM1 = LOG(XYM1) 
 
	XLNYXM1 = LOG(YX - 1.) 
 
	XLNY2 = XLNY**2 
 
	XLNX2 = XLNX**2 
 
	XLNXY2 = XLNXY**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNYM12 = XLNYM1**2 
 
	XLNXYM12 = XLNXYM1**2 
 
	XLNYXM12 = XLNYXM1**2 
 
	GF1 = GF1BLFA(X,Y) 
 
	GF1B = GF1BLFA(XM1,YM1) 
 
	HGF = HGFBLA(X,Y) 
 
	HGFB = HGFBLBA(XM1,Y) 
 
	CF2 = -3.*XY*(1.+ X) + (X*XOYM1 + 2.*XY*XM1 - 2.*XY*X)*XLNXY  
     >  - (XY*XMYM + XY + XMYM)*XLNXYM1 + GF1*XLNXY2   
     >  - (GF1 + GF1B)*XLNXYM12 -(X*XOYM1 + 1./Y - XM1*XM1OY)*XLNX 
     >   - GF1*XLNX2 +(XM1*XM1OY + 1./YM1 - X*XOYM1)*XLNXM1  
     >  + GF1B*XLNXM12  
 
	CFBO = - 2./3.*GF1 + 2.*XY + (GF1 + GF1B)*XLNXYM1 + GF1*XLNX  
     >  - GF1B*XLNXM1 
 
	CFCG = (13./3. - PI2/3.)*GF1 - XY - (XY*(XM1 - X) + 5.*X*XMYM) 
     >  *XLNXY + (XMYM - 2.*XM1*XMYM + XY - 2.*X*XY)*XLNXYM1  
     >  - 1./2.*(GF1B + 3.*X*XOYM1)*XLNXY2 
     >  + 1./2.*(GF1 + GF1B)*XLNYXM12 - 1./2.*HGF + 1./2.*HGFB   
     >  + XLNX*(5.*X*XMYM + XY + GF1*(- 2. + 1./2.*XLNX - XLNXM1) 
     >  + 3./2.*X*XOYM1*XLNX) - XLNXM1*(5.*XM1*XY + XMYM  
     >  + GF1B*(- 2.+ 1./2.*XLNXM1 - XLNX) + 3./2.*XM1*XM1OY*XLNXM1)   
 
	GF2BLA = CF**2/2.*CF2 - CF*BO/2.*CFBO + CF*CG*CFCG 
 
	RETURN 
 
C	----------------------- 
 
	ENTRY GG2BLA (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	YX = Y/X 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XOYM1 = X/YM1 
 
	XM1OY = XM1/Y 
 
	XMYM = XM1/YM1 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXY = LOG(XY) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXYM1 = LOG(XYM1) 
 
	XLNYXM1 = LOG(YX - 1.) 
 
	XLNY2 = XLNY**2 
 
	XLNX2 = XLNX**2 
 
	XLNXY2 = XLNXY**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNYM12 = XLNYM1**2 
 
	XLNXYM12 = XLNXYM1**2 
 
	XLNYXM12 = XLNYXM1**2 
 
	GG1 = GG1BLFA(X,Y) 
 
	GG1B = GG1BLFA(XM1,YM1) 
 
	HGG = HGGBLA(X,Y) 
 
	HGGB = HGGBLBA(XM1,Y) 
 
	CA2 = (2./3. - PI2/3.)*GG1 - 15./4.*XY2 - 6.*XY*XMYM + XY2/YM1  
     >  + XY/Y - (2.*XY/YM1 - 6.*XOYM1*XMYM + XOYM1**2)*XLNXY  
     >  + (GG1 + 2.*XOYM1**2)*XLNXY2 
     >  - (GG1 - GG1B)*XLNXY*XLNXYM1 - 1./2.*HGG + 1./2.*HGGB  
     >  + XLNX*(3.*XY*XMYM -3.*XY*XOYM1 - 6.*XMYM*XOYM1 + XOYM1**2  
     >  - 2.*XOYM1**2*(XLNX + XLNXM1) + (3.*XOYM1/YM1 - Y/YM1**2) 
     >  *XLNXM1) + XLNXM1*(3.*XY*XMYM - 3.*XMYM*XM1OY - 6.*XY*XM1OY  
     >  + XM1OY**2 
     >  - 2.*XM1OY**2*(XLNX + XLNXM1) + (3.*XM1OY/Y - YM1/Y2)*XLNX) 
     
	CGBO = - 5./3.*GG1 - 5./2.*XY2 + XY/Y + XY*XOYM1  
     >  + XOYM1**2*XLNY  
     >  + XM1OY**2*XLNXM1 
 
	CFNFTF = 2.*XY/(Y*YM1) - 4.*XY2 + 2.*XY*XOYM1*XLNXY  
     >  - XOYM1*XLNXY*(2. + X*XLNXY)/YM1 + XOYM1*XLNX* 
     >  (2. + X*XLNX)/YM1 + XM1OY*XLNXM1*(2. + XM1*XLNXM1)/Y  
      
 
	GG2BLA = CG**2*CA2 + CG*BO/2.*CGBO + CF*TFNF*CFNFTF 
 
	RETURN 
 
C	------------------------------ 
 
	ENTRY SFF2BL (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	YX = Y/X 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XOYM1 = X/YM1 
 
	XM1OY = XM1/Y 
 
	XMYM = XM1/YM1 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXY = LOG(XY) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXYM1 = LOG(XYM1) 
 
	XLNYXM1 = LOG(YX - 1.) 
 
	XLNY2 = XLNY**2 
 
	XLNX2 = XLNX**2 
 
	XLNXY2 = XLNXY**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNYM12 = XLNYM1**2 
 
	XLNXYM12 = XLNXYM1**2 
 
	XLNYXM12 = XLNYXM1**2 
		  
	SPENXM1 = SPENCE(XM1) 
 
	SPENYM1 = SPENCE(YM1) 
 
	CFTFNF = 3.*XY - 8.*XY*X*YM1 + (5.*XOYM1 - 8.*X*XOYM1)*XLNXY  
     >  + (XOYM1 - 4.*X*XM1)*XLNXY2 + 8.*X*XM1*(SPENXM1 - SPENYM1 
     >  + XLNXM1*XLNY) - 290./9.*X*XM1 - (5.*XOYM1 - 8.*X*XOYM1   
     >  + 6.*X - 38./3.*XM1*X)*XLNX 
     >  - (XOYM1 - 4.*X*XM1)*XLNX2 + 4.*X*XM1*XLNX*XLNXM1   
     >  - 290./9.*X*XM1 - (5.*XM1OY - 8.*XM1*XM1OY + 6.*XM1  
     >  - 38./3.*X*XM1)*XLNXM1 - (XM1OY - 4.*X*XM1)*XLNXM12  
     >  + 4.*X*XM1*XLNX*XLNXM1 
 
	SFF2BL = 2.*CF*TFNF*CFTFNF 

	RETURN 
 
C	------------------------------------- 
 
	ENTRY FG2BL (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	YX = Y/X 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XOYM1 = X/YM1 
 
	XM1OY = XM1/Y 
 
	XMYM = XM1/YM1 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXY = LOG(XY) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXYM1 = LOG(XYM1) 
 
	XLNYXM1 = LOG(YX - 1.) 
 
	XLNY2 = XLNY**2 
 
	XLNX2 = XLNX**2 
 
	XLNXY2 = XLNXY**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNYM12 = XLNYM1**2 
 
	XLNXYM12 = XLNXYM1**2 
 
	XLNYXM12 = XLNYXM1**2 
 
	FG1A = FG1BLFA(X,Y) 
 
	FG1BA = FG1BLFA(XM1,YM1) 
 
	FG1 = FG1BLF(X,Y) 
 
	FG1B = FG1BLF(XM1,YM1) 
 
	HFG = HFGBL(X,Y) 
 
	HFGB = HFGBLB(XM1,Y) 
 
	CFNFTF = 2.*(5. - PI2/3.)*FG1 - 4.*FG1A + XY/YM1  
     >  + 2.*(FG1 + FG1B  
     >  + XY/YM1**2 - 1./YM1**2 + 1./2.*XOYM1/YM1)*XLNXY  
     >  - 2.*(FG1 + FG1B 
     >  + 1./(Y*YM1))*XLNXYM1 + FG1BA*XLNXY2  
     >  + (FG1 + FG1B)*XLNYXM12  
     >  - 2.*(FG1 + FG1B + FG1BA + 1./2.*XOYM1/YM1  
     >  + 2.*XY)*XLNX  
     >  - 2.*FG1*XLNX*XLNXM1 + (FG1 - FG1B - FG1BA)*XLNX2   
     >  + 2.*(FG1B + FG1 + FG1A + 1./2.*XM1OY/Y   
     >  + 2.*XMYM)*XLNXM1  
     >  + 2.*FG1B*XLNX*XLNXM1 - (FG1B - FG1 - FG1A)*XLNXM12 
 
	CGNFTF = 6.*XOYM1*(3. - 4.*X) - 2.*XY*(11. - 16.*X)  
     >  + 4.*XY*(1.  
     >  - 3.*X)/Y - 2.*(FG1 - 5.*FG1B + (5.-7.*X)/YM1**2  
     >  + 2.*XOYM1*( 
     >  3./Y - 2.*XY - 12.*XM1))*XLNXY + 2.*(FG1 + FG1B  
     >  + 1/(Y*YM1)) 
     >  *XLNXYM1 + (FG1 - FG1B + (1. - 4.*X)/YM1**2)*XLNXY2  
     >  - (FG1 + FG1B)*XLNXYM12 - HFG + HFGB + 2.*(FG1  
     >  - 5.*FG1B   
     >  + (5. - 7.*X)/YM1**2 + XOYM1*(5./Y - 6.*XY - 22.  
     >  + 28.*X))*XLNX  
     >  - (FG1 - FG1B + (1. - 4.*X)/YM1**2)*XLNX2  
     >  - 2.*(FG1B - 5.*FG1 + (5. - 7.*XM1)/Y**2  
     >  + XM1OY*(5./YM1 - 6.*XMYM - 22. + 28.*XM1))*XLNXM1  
     >  + (FG1B - FG1 + (1. - 4.*XM1)/Y**2)*XLNXM12 
 
	FG2BL = CF*TFNF*CFNFTF + CG*TFNF*CGNFTF  

	RETURN 
 
C	----------------------------------------- 
 
	ENTRY GF2BL (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	YX = Y/X 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XOYM1 = X/YM1 
 
	XM1OY = XM1/Y 
 
	XMYM = XM1/YM1 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXY = LOG(XY) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXYM1 = LOG(XYM1) 
 
	XLNYXM1 = LOG(YX - 1.) 
 
	XLNY2 = XLNY**2 
 
	XLNX2 = XLNX**2 
 
	XLNXY2 = XLNXY**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNYM12 = XLNYM1**2 
 
	XLNXYM12 = XLNXYM1**2 
 
	XLNYXM12 = XLNYXM1**2 
 
	GF1A = GF1BLFA(X,Y) 
 
	GF1BA = GF1BLFA(XM1,YM1) 
 
	GF1 = GF1BLF(X,Y) 
 
	GF1B = GF1BLF(XM1,YM1) 
 
	HGF = HGFBL(X,Y) 
 
	HGFB = HGFBLB(XM1,Y)  
 
	CF2 = 3.*GF1A - 9./2.*XY - (1./2.*XOYM1*(2.+ X) - XY*XMYM)*XLNXY 
     > - (3./2.*GF1 + 3./2.*GF1B + XY*XMYM)*XLNXYM1 - 1./2.*GF1A*XLNXY2 
     > - 1./2.*(GF1 + GF1B)*XLNXYM12 - 1./2.*(XY*(2.- 5.*X)  
     > - 3.*X*XOYM1)*XLNX + 1./2.*GF1A*XLNX2 + 1./2.*(XMYM*(2. - 5.*XM1)  
     > - 3.*XM1*XM1OY)*XLNXM1 - 1./2.*GF1BA*XLNXM12  
 
	CFBO = 10./3.*GF1 + 2.*XY*XM1 + (GF1 + GF1B)*XLNXYM1 - GF1A*XLNX 
     > + GF1BA*XLNXM1 
 
	CFCG = (4./9. - PI2/3.)*GF1 - 26./9.*X*XY + XY*(3. + 4.*X2) 
     > - X*(105. - 246.*X + 188.*X2)/9. - (XY*XMYM - XOYM1*(6. - 11.*X  
     > + 8.*X2) + 4.*X*XM1*(XM1 - X))*XLNXY + XY*XMYM*XLNXYM1  
     > + 1./2.*(2.*GF1 + 2.*GF1A + GF1BA - 3.*X*XOYM1 - 2.)*XLNXY2 
     >  + 1./2.*(2. - GF1A - GF1BA)*XLNYXM12 
     > - 1./2.*HGF - 1./2.*HGFB + (XY*XMYM - XOYM1*(6. - 11.*X + 8.*X2) 
     > + X*(18. - 63.*X + 62.*X2)/3.)*XLNX - GF1*XLNX*XLNXM1 - XY*XOYM1* 
     > (1. - 4.*Y)*XLNX2/2. -(XY*XMYM - XM1OY*(6.- 11.*XM1 + 8.*XM1**2) 
     > + XM1*(18. - 63.*XM1 + 62.*XM1**2)/3.)*XLNXM1 + GF1B*XLNX*XLNXM1  
     > + XMYM*XM1OY*(1. - 4.*YM1)*XLNXM12/2.  
 
	GF2BL = CF**2*CF2 - CF*BO/2.*CFBO + CF*CG*CFCG 

	RETURN 
 
C	--------------------------------- 
 
	ENTRY GG2BL (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	Y = YY 
 
	XY = X/Y 
 
	YX = Y/X 
 
	XM1 = 1.- X 
 
	YM1 = 1.- Y 
 
	XOYM1 = X/YM1 
 
	XM1OY = XM1/Y 
 
	XMYM = XM1/YM1 
 
	X2 = X**2 
 
	Y2 = Y**2 
 
	XY2 = X2/Y2 
 
	XYM1 = 1. - XY 
 
	XLNY = LOG(Y) 
 
	XLNX = LOG(X) 
 
	XLNXY = LOG(XY) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNYM1 = LOG(YM1) 
 
	XLNXYM1 = LOG(XYM1) 
 
	XLNYXM1 = LOG(YX - 1.) 
 
	XLNY2 = XLNY**2 
 
	XLNX2 = XLNX**2 
 
	XLNXY2 = XLNXY**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNYM12 = XLNYM1**2 
 
	XLNXYM12 = XLNXYM1**2 
 
	XLNYXM12 = XLNYXM1**2 
 
	GG1A = GG1BLFA(X,Y) 
 
	GG1BA = GG1BLFA(XM1,YM1) 
 
	GG1 = GG1BLF(X,Y) 
 
	GG1B = GG1BLF(XM1,YM1) 
 
	GG1C = GG1BLCF(X,Y) 
 
	HGG = HGGBL(X,Y) 
 
	HGGB = HGGBLB(XM1,Y)  
 
	CG2 = (2./3. - PI2/3.)*GG1 - 1./2.*GG1C + XY*(1./Y + 5./4.*XY)  
     > - XY*(10. - 19.*X + 16.*X2) + XOYM1*(6. - 13.*X + 8.*X2) +  
     > (2.*XY*(1.- 2.*XM1*(X*Y + XM1*YM1))/YM1 + XOYM1*(2.- 3.*X  
     > - 8.*XM1**2)/YM1)*XLNXY + (GG1 + 2.*XOYM1**2)*XLNXY2  
     > - (GG1 - GG1B)*XLNXY*XLNXYM1 
     > - 1./2.*HGG - 1./2.*HGGB - 2.*XOYM1**2*XLNX2   
     > + (XMYM*(X2 + XM1**2 - 2.*XM1*(1. + 2.*X)*YM1)/YM1 + 1./YM1) 
     > *XLNX*XLNXM1   
     > + (XY*(3.*(X2 + XM1**2)*(YM1 - Y) - 4.*X*Y)/YM1  
     > - XOYM1*(2. - 3.*X - 8.*XM1**2)/YM1)*XLNX 
     > - 2.*XM1OY**2*XLNXM12 + (XY*(XM1**2 + X2 -2.*X*(1.+2.*XM1)*Y)/Y  
     > + 1./Y)*XLNX*XLNXM1 + (XMYM*(3.*(X2+ XM1**2)*(Y - YM1)  
     > - 4.*XM1*YM1)/Y - XM1OY*(2. - 3.*XM1 - 8.*X2)/Y)*XLNXM1 
      
 
	CGBO = - 5./3.*GG1 - 13./3.*GG1C - 11./2.*XY2 + XY2/YM1 + XY/Y 
     > + XOYM1**2*XLNY + XM1OY**2*XLNXM1 
 
	CFTFNF = - 20./3.*GG1C - 12.*XY2 + 2.*XY*(2. + Y + 3.*X*Y - 4.*Y2) 
     > /(Y*YM1) - 2.*XY*(X - Y + 2.*X*Y)*XLNXY/YM1**2 - XOYM1**2*XLNXY2 
     > - 2.*XOYM1*(1. - 3.*X - 2.*X*YM1)*XLNX/YM1 + XOYM1**2*XLNX2 
     > - 2.*XM1OY*(1. - 3.*XM1 - 2.*XM1*Y)*XLNXM1/Y + XM1OY**2*XLNXM12 
 
	GG2BL = CG**2*CG2 + CG*BO/2.*CGBO + CF*TFNF*CFTFNF 

	RETURN 
 
C	------------------------- 
C	Now the DGLAP kernels! 
C	------------------------- 
 
	ENTRY PNSMD (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	D = YY 
 
	XD = X/D 
 
	DX = D/X 
 
	XM1 = 1.- X 
 
	DM1 = 1.- D 
 
	XODM1 = X/DM1 
 
	XM1OD = XM1/D 
 
	XMDM = XM1/DM1 
 
	X2 = X**2 
 
	D2 = D**2 
 
	XD2 = X2/D2 
 
	XDM1 = 1. - XD 
 
	XMDM1 = 1. - XMDM 
 
	XLND = LOG(D) 
 
	XLNX = LOG(X) 
 
	XLNXD = LOG(XD) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNDM1 = LOG(DM1) 
 
	XLNXDM1 = LOG(XDM1) 
 
	XLNDXM1 = LOG(DX - 1.) 
 
	XLNXMDM = LOG(XMDM) 
 
	XLNXMDM1 = LOG(XMDM1) 
 
	XLND2 = XLND**2 
 
	XLNX2 = XLNX**2 
 
	XLNXD2 = XLNXD**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNDM12 = XLNDM1**2 
 
	XLNXDM12 = XLNXDM1**2 
 
	XLNDXM12 = XLNDXM1**2 
 
	XLNXMDM12 = XLNXMDM1**2 
 
	XLNXMDM2 = XLNXMDM**2 
 
	FF1 = FF1F(X,D) 
 
	HFF = HFFD(X,D) 
 
	HFFB = HFFDB(X,D)  
 
	CF2 = 1./2.*XLNXMDM12*((1. - XD*XMDM)/XM1 + X*XMDM1/XM1)  
     > + 1./2.*XLNX2*(2.*X*XMDM1/XM1 + XD/DM1) - FF1*(XLNXMDM1*XLNXMDM 
     > + XLNXM1*XLNX) - 2.*XLNXMDM1*XMDM1*(1. - XM1 - DM1 - X*D/4.) 
     > /(D*XM1) - 1./2.*XLNX*(4.*XD/DM1 + 3.*X*XMDM1/XM1) + 3.*XMDM 
     > + (4./3. - PI2/3.)*FF1 
 
	CFBO = 5./3.*FF1 + XMDM + XLNXMDM1*(1. - XD*XMDM)/XM1  
     > + XLNX*(X*XMDM1/XM1 + XD/DM1) 
 
	CFCFMCG = 2.*XMDM + 4./3.*FF1 + HFF + HFFB 
 
	PNSMD = X*(CF**2*CF2 - CF*BO/2.*CFBO - CF*(CF - CG/2.)*CFCFMCG) 

	RETURN 
 
C	------------------------ 
 
	ENTRY PNSPD (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	D = YY 
 
	XD = X/D 
 
	DX = D/X 
 
	XM1 = 1.- X 
 
	DM1 = 1.- D 
 
	XODM1 = X/DM1 
 
	XM1OD = XM1/D 
 
	XMDM = XM1/DM1 
 
	X2 = X**2 
 
	D2 = D**2 
 
	XD2 = X2/D2 
 
	XDM1 = 1. - XD 
 
	XMDM1 = 1. - XMDM 
 
	XLND = LOG(D) 
 
	XLNX = LOG(X) 
 
	XLNXD = LOG(XD) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNDM1 = LOG(DM1) 
 
	XLNXDM1 = LOG(XDM1) 
 
	XLNXMDM = LOG(XMDM) 
 
	XLNXMDM1 = LOG(XMDM1) 
 
	XLND2 = XLND**2 
 
	XLNX2 = XLNX**2 
 
	XLNXD2 = XLNXD**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNDM12 = XLNDM1**2 
 
	XLNXDM12 = XLNXDM1**2 
 
	XLNXMDM12 = XLNXMDM1**2 
 
	XLNXMDM2 = XLNXMDM**2 
 
	FF1 = FF1F(X,D) 
 
	HFF = HFFD(X,D) 
 
	HFFB = HFFDB(X,D)  
 
	CF2 = 1./2.*XLNXMDM12*((1. - XD*XMDM)/XM1 + X*XMDM1/XM1)  
     > + 1./2.*XLNX2*(2.*X*XMDM1/XM1 + XD/DM1) - FF1*(XLNXMDM1*XLNXMDM 
     > + XLNXM1*XLNX) - 2.*XLNXMDM1*XMDM1*(1. - XM1 - DM1 - X*D/4.) 
     > /(D*XM1) - 1./2.*XLNX*(4.*XD/DM1 + 3.*X*XMDM1/XM1) + 3.*XMDM 
     > + (4./3. - PI2/3.)*FF1 
 
	CFBO = 5./3.*FF1 + XMDM + XLNXMDM1*(1. - XD*XMDM)/XM1  
     > + XLNX*(X*XMDM1/XM1 + XD/DM1) 
 
	CFCFMCG = 2.*XMDM + 4./3.*FF1 + HFF - HFFB 
 
	PNSPD = X*(CF**2*CF2 - CF*BO/2.*CFBO - CF*(CF - CG/2.)*CFCFMCG) 

	RETURN 
 
C	------------------------ 
 
	ENTRY SFF2DA (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	D = YY 
 
	XD = X/D 
 
	XM1 = 1.- X 
 
	DM1 = 1.- D 
 
	XMDM = XM1/DM1 
 
	X2 = X**2 
 
	D2 = D**2 
 
	XD2 = X2/D2 
 
	XDM1 = 1. - XD 
 
	XMDM1 = 1. - XMDM 
 
	XLNX = LOG(X) 
 
	XLNXMDM1 = LOG(XMDM1) 
 
	CFNFTF = X*(XDM1*XLNXMDM1*(XLNXMDM1 - 1.)  
     > + XD*XLNX*(XLNX - 1.)/DM1) 
 
	SFF2DA = -2.*CF*TFNF*CFNFTF 
 
	RETURN 
 
C	------------------------- 
 
	ENTRY FG2DA (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	D = YY 
 
	XD = X/D 
 
	DX = D/X 
 
	XM1 = 1.- X 
 
	DM1 = 1.- D 
 
	DMXM1 = DM1/XM1 - 1. 
 
	XMDM = XM1/DM1 
 
	X2 = X**2 
 
	D2 = D**2 
 
	XD2 = X2/D2 
 
	XDM1 = 1. - XD 
 
	XMDM1 = 1. - XMDM 
 
	XLNDMXM1 = LOG(DMXM1) 
 
	XLND = LOG(D) 
 
	XLNX = LOG(X) 
 
	XLNXMX = LOG(XM1/X) 
 
	XLNXD = LOG(XD) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNDM1 = LOG(DM1) 
 
	XLNDXM1 = LOG(1.- DX) 
 
	XLNXMDM = LOG(XMDM) 
 
	XLNXMDM1 = LOG(XMDM1) 
 
	XLND2 = XLND**2 
 
	XLNX2 = XLNX**2 
 
	XLNXD2 = XLNXD**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNDM12 = XLNDM1**2 
 
	XLNDXM12 = XLNDXM1**2 
 
	XLNXMDM12 = XLNXMDM1**2 
 
	XLNXMDM2 = XLNXMDM**2 
 
	XLNDMXM12 = XLNDMXM1**2 
 
	XLNXMX2 = XLNXMX**2 
 
	FG1 = FG1FA(X,D) 
 
	HFG = HFGDA(X,D) 
 
	HFGB = HFGDBA(X,D) 
 
	CGNFTF = -(3./2.*FG1 + 4.*XDM1/DM1**2 + 3.)*XLNXMDM12  
     > - 1./2.*FG1*(XLNXMDM2 + XLNXM12) - XLNX2*XD*(1. + 3./DM1**2) 
     > + 2.*XLNXMDM1*(FG1 + 2. 
     > + 3.*XDM1/DM1**2) + 2.*XLNX*XD*(1. + 2./DM1**2)  
     > + 2.*(XLNXMDM + XLNXM1) 
     > *(1./2.*FG1 - 1./DM1) - 2.*FG1 + 6./DM1 - 1./2.*HFG  
     > - 1./2.*HFGB 
 
	CFNFTF = XD*(XLNXMDM12 - XLNX2) + 1./2.*FG1*(XLNDMXM12  
     > + XLNXMX2 - XLNX2) 
     > + (XLNX + XLNDXM1)*XMDM1*(7./D + 3.*D - 8.)/DM1  
     > + XLNDM1*(3.*XDM1 - 2.*XD*(1. + DM1)/DM1) + 4.*XLNXM1*XMDM 
     > *(1. + DM1)/DM1 - XLNX*(XD*(2. + 5./DM1) + 3.*X/DM1**2)  
     > + (5. - PI2/3.)*FG1 - 5./DM1 
 
	FG2DA = X*(CG*TFNF*CGNFTF + CF*TFNF*CFNFTF) 
 
	RETURN 
 
C	-------------------------- 
 
	ENTRY GF2DA (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	D = YY 
 
	XD = X/D 
 
	DX = D/X 
 
	XM1 = 1.- X 
 
	DM1 = 1.- D 
 
	DMXM1 = DM1/XM1 - 1. 
 
	XMDM = XM1/DM1 
 
	X2 = X**2 
 
	D2 = D**2 
 
	XD2 = X2/D2 
 
	XDM1 = 1. - XD 
 
	XMDM1 = 1. - XMDM 
 
	XLNDMXM1 = LOG(DMXM1) 
 
	XLND = LOG(D) 
 
	XLNX = LOG(X) 
 
	XLNXMX = LOG(XM1/X) 
 
	XLNXD = LOG(XD) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNDM1 = LOG(DM1) 
 
	XLNDXM1 = LOG(1.- DX) 
 
	XLNXMDM = LOG(XMDM) 
 
	XLNXMDM1 = LOG(XMDM1) 
 
	XLND2 = XLND**2 
 
	XLNX2 = XLNX**2 
 
	XLNXD2 = XLNXD**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNDM12 = XLNDM1**2 
 
	XLNDXM12 = XLNDXM1**2 
 
	XLNXMDM12 = XLNXMDM1**2 
 
	XLNXMDM2 = XLNXMDM**2 
 
	XLNDMXM12 = XLNDMXM1**2 
 
	XLNXMX2 = XLNXMX**2 
 
	GF1 = GF1FA(X,D) 
 
	HGF = HGFDA(X,D) 
 
	HGFB = HGFDBA(X,D)	 
 
	CF2  = 1./2.*XLNXMDM12*(GF1 - X*XD) - 1./2.*GF1*(XLNXMDM2  
     > + XLNXM12) + X*XD*XLNX2/2. - 1./2.*(XLNX + XLNDXM1)*XMDM1 
     > *(X - D + 3. - 5.*XD) + 1./2.*XLNDM1*(XMDM1*XDM1*(4. - D)  
     > + X*(2. - XD)) - 1./2.*XLNX*(5.*XD*XMDM1*(1. - 4./5.*D)  
     > + X*(2. + 1./DM1)) - XLNXM1*(XMDM1*(2.  
     > + X - D) + D) - 3./2.*GF1 - 3./2.*X - 3./2.*XMDM1 
 
	CFNFTF = GF1*(XLNXMDM + XLNXM1 - 2./3.) + 2.*X + 2.*XMDM1 
 
	CGCF = 1./2.*XLNXMDM12*(3.*GF1*DM1 - X*XD*(1. + 3.*DM1))  
     > + 1./2.*XLNX2*(XD*X*(4. - D)/DM1 - GF1) + 1./2.*GF1 
     > *(XLNDMXM12 + XLNXMX2) + XLNXMDM1*XMDM1*(4. - 5.*X  
     > + 3.*DM1*XMDM1/D) + 1./6.*(XLNXMDM + XLNXM1)*(5.*X  
     > + XMDM1*(5. + X)) - XLNX*(4.*XD*XMDM1 
     > + X*(XMDM1 - XD)) + (28./9. - PI2/3.)*GF1   
     > + 8./3.*(XMDM1 + X) + 1./2.*HGFB - 1./2.*HGF 
 
	GF2DA = CF**2*CF2 - 2./3.*CF*TFNF*CFNFTF + CG*CF*CGCF 
 
	RETURN 
 
C	--------------------------- 
 
	ENTRY GG2DA (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	D = YY 
 
	XD = X/D 
 
	DX = D/X 
 
	XM1 = 1.- X 
 
	DM1 = 1.- D 
 
	DMXM1 = DM1/XM1 - 1. 
 
	XMDM = XM1/DM1 
 
	X2 = X**2 
 
	D2 = D**2 
 
	XD2 = X2/D2 
 
	XODM1 = X/DM1 
 
	XDM1 = 1. - XD 
 
	XMDM1 = 1. - XMDM 
 
	XLNDMXM1 = LOG(DMXM1) 
 
	XLND = LOG(D) 
 
	XLNX = LOG(X) 
 
	XLNXMX = LOG(XM1/X) 
 
	XLNXD = LOG(XD) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNDM1 = LOG(DM1) 
 
	XLNDXM1 = LOG(1.- DX) 
 
	XLNXMDM = LOG(XMDM) 
 
	XLNXMDM1 = LOG(XMDM1) 
 
	XLND2 = XLND**2 
 
	XLNX2 = XLNX**2 
 
	XLNXD2 = XLNXD**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNDM12 = XLNDM1**2 
 
	XLNDXM12 = XLNDXM1**2 
 
	XLNXMDM12 = XLNXMDM1**2 
 
	XLNXMDM2 = XLNXMDM**2 
 
	XLNDMXM12 = XLNDMXM1**2 
 
	XLNXMX2 = XLNXMX**2 
 
	GG1 = GG1FA(X,D) 
 
	GG1A = GGA1F(X,D) 
 
	HGG = HGGDA(X,D) 
 
	HGGB = HGGDBA(X,D) 
 
	CGNFTF = XLNXMDM1*XDM1**2*D - XODM1**2*XLNX/D - 8./3.*GG1  
     > + X*XMDM1*(XMDM1 + X)/XM1 - 3./2.*GG1A + 2.*XODM1  
     > - D*(1. - X*D)/DM1**2 
 
	CFNFTF = XLNXMDM12*D*XDM1**2 - XODM1**2*XLNX2/D  
     > - 2.*XLNXMDM1*(XMDM1 - D)*XDM1 - 2.*XLNX*(XD*XMDM1  
     > + X*(2. - X)/DM1)/DM1 - 8.*XODM1*XMDM 
     > + 4.*D*(1. - X2)/DM1**2 - 2.*(D/DM1)**2*XM1 
 
	CG2 = 1./2.*XLNXMDM12*(1. + 9.*X - 4.*D - 2.*(4. - 3.*X2)/DM1 
     > - 8.*XMDM1**2/D + XMDM*(7. - 6.*X)/DM1) + XLNX2*XODM1**2* 
     > (4.*DM1/D + 2.*D + DM1**2/XM1) + 1./2.*GG1*(XLNDMXM12  
     > - XLNXMDM2 - 2.*XLNX*XLNXM1) - XLNXMDM1*XDM1*(2.*D/DM1  
     > + 31./6.*X + 5./6.*D) + XLNX*XODM1*(6.  
     > + 2.*DM1 - 31./6.*XD)/DM1 + (67./18. - PI2/3.)*GG1  
     > + 5./6.*(D/DM1)**2*XM1 + 11./3.*XODM1*XMDM  
     > - 11./6.*D*(1. - X2)/DM1**2 + 1./2.*HGGB - 1./2.*HGG 
 
	GG2DA = 2./3.*CG*TFNF*CGNFTF + CF*TFNF*CFNFTF + CG**2*CG2 
 
	RETURN 
 
C	-------------------------------- 
 
	ENTRY SFF2D (N,XX,YY) 
 
	TFNF = TF*N 
 
	X = XX 
 
	D = YY 
 
	XD = X/D 
 
	XM1 = 1.- X 
 
	DM1 = 1.- D 
 
	XMDM = XM1/DM1 
 
	X2 = X**2 
 
	XDM1 = 1. - XD 
 
	XMDM1 = 1. - XMDM 
 
	TOD = 1. - 1./D 
 
	XLND = -LOG(1./D) 
 
	XLNX = -LOG(1./X) 
 
	XLNXMDM1 = -LOG(1./XMDM1) 
 
	XLNXMDM12 = XLNXMDM1**2 
 
	XLNX2 = XLNX**2 
 
c	print *, XDM1 
 
	SPENXDM1 = SPENCE(XDM1)  
 
	SPEN1OMD = SPENCE(TOD) 
 
	CFNFTF = XDM1*(1. - 4.*XD/D)*XLNXMDM12   
     > + 8.*XD*XDM1*XLNXMDM1*(2.*XLND - XLNX)/D  
     > - 16.*XD*XDM1*(SPENXDM1 - SPEN1OMD)/D   
     > + XLNX2*XD*(3.- 4.*XD + 1./DM1)/D  
     > - XLNXMDM1*XDM1*(3. - 8.*XD)  
     > + XLNX*XD*(5. - 8.*XD)/DM1 - 3.*XMDM  
     > + 8.*(1. - X2)/(D*DM1)  
     > - 16.*XD*XMDM/D 
 
	SFF2D = -2.*CF*TFNF*CFNFTF*X 
  
	RETURN 
 
C	-------------------------------- 
 
	ENTRY FG2D (N,XX,YY) 
 
	TFNF = TF*N 
 
	X = XX 
 
	D = YY 
 
	XD = X/D 
 
	DX = D/X 
 
	XM1 = 1.- X 
 
	DM1 = 1.- D 
 
	DMXM1 = DM1/XM1 - 1. 
 
	XMDM = XM1/DM1 
 
	X2 = X**2 
 
	D2 = D**2 
 
	XD2 = X2/D2 
 
	XODM1 = X/DM1 
 
	XDM1 = 1. - XD 
 
	XMDM1 = 1. - XMDM 
 
	XLNDMXM1 = LOG(DMXM1) 
 
	XLND = LOG(D) 
 
	XLNX = LOG(X) 
 
	XLNXMX = LOG(XM1/X) 
 
	XLNXD = LOG(XD) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNDM1 = LOG(DM1) 
 
	XLNDXM1 = LOG(1.- DX) 
 
	XLNXMDM = LOG(XMDM) 
 
	XLNXMDM1 = LOG(XMDM1) 
 
	XLND2 = XLND**2 
 
	XLNX2 = XLNX**2 
 
	XLNXD2 = XLNXD**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNDM12 = XLNDM1**2 
 
	XLNDXM12 = XLNDXM1**2 
 
	XLNXMDM12 = XLNXMDM1**2 
 
	XLNXMDM2 = XLNXMDM**2 
 
	XLNDMXM12 = XLNDMXM1**2 
 
	XLNXMX2 = XLNXMX**2 
 
	FG1 = FG1F(X,D) 
 
	FG1A = FG1FA(X,D) 
 
	FG1C = FGC1F(X,D) 
 
	HFG = HFGD(X,D) 
 
	HFGB = HFGDB(X,D) 
 
	CGNFTF = -2.*XLNXMDM12*(FG1A + 3./2. + XDM1* 
     > (1. + 6.*XD - 4.*XD/D) 
     > /DM1**2) - 1./2.*FG1*(XLNXMDM12 + XLNXMDM2  
     > + XLNXM12 - XLNX2)  
     > + 2.*XLNX2*(XDM1*(1. + 6.*XD - 4.*XD/D)  
     > - 3./2.)/DM1**2  
     > + 2.*XLNXMDM1*(6.*FG1 - 7.*FG1C + 2.  
     > + XDM1*(3. - 2.*D + 20.*X 
     > - 8.*XD)/DM1**2) + (XLNXMDM + XLNXMX)* 
     > (FG1 - 2./DM1) + 2.*XLNX* 
     > (1. + D - XDM1*(1. + 4.*X + 8.*XD)  
     > + 2.*DM1*XD)/DM1**2 
     > - 2.*(FG1 + FG1C) - 32.*XMDM*XDM1/D  
     > + 2.*(3. + 8.*XM1/D)/DM1 
     > + 1./2.*HFGB - 1./2.*HFG 
 
	CFNFTF = - XD*XLNXMDM12 + 1./2.*FG1*(XLNDMXM12  
     > + XLNXMX2)  
     > - XLNX2*XDM1/DM1**2 +(XLNX + XLNDXM1)*XMDM1 
     > *(4.*X - D - 1./D)/DM1  
     > - XLNDM1*(XDM1 + 2.*X/DM1) + XLNX*(XD  
     > + 3.*X*XMDM1/DM1  
     > + (X/DM1)**2) + 4.*XLNXM1*XMDM*(2.*X - D)/DM1  
     > - PI2/3.*FG1  
     > - 20.*XMDM*X/DM1 + 8.*XMDM/DM1 + 6.*D*XMDM/DM1  
     > + 5./DM1 
 
	FG2D = X*(CG*TFNF*CGNFTF + CF*TFNF*CFNFTF) 

	RETURN 
 
C	------------------------------------- 
 
	ENTRY GF2D (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	D = YY 
 
	XD = X/D 
 
	DX = D/X 
 
	XM1 = 1.- X 
 
	DM1 = 1.- D 
 
	DMXM1 = DM1/XM1 - 1. 
 
	XMDM = XM1/DM1 
 
	X2 = X**2 
 
	D2 = D**2 
 
	XD2 = X2/D2 
 
	XODM1 = X/DM1 
 
	XDM1 = 1. - XD 
 
	XMDM1 = 1. - XMDM 
 
	XLNDMXM1 = LOG(DMXM1) 
 
	XLND = LOG(D) 
 
	XLNX = LOG(X) 
 
	XLNXMX = LOG(XM1/X) 
 
	XLNXD = LOG(XD) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNDM1 = LOG(DM1) 
 
	XLNDXM1 = LOG(1.- DX) 
 
	XLNXMDM = LOG(XMDM) 
 
	XLNXMDM1 = LOG(XMDM1) 
 
	XLND2 = XLND**2 
 
	XLNX2 = XLNX**2 
 
	XLNXD2 = XLNXD**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNDM12 = XLNDM1**2 
 
	XLNDXM12 = XLNDXM1**2 
 
	XLNXMDM12 = XLNXMDM1**2 
 
	XLNXMDM2 = XLNXMDM**2 
 
	XLNDMXM12 = XLNDMXM1**2 
 
	XLNXMX2 = XLNXMX**2 
 
	GF1 = GF1F(X,D) 
 
	GF1A = GF1FA(X,D) 
 
	HGF = HGFD(X,D) 
 
	HGFB = HGFDB(X,D) 
 
	CF2 = 1./2.*XLNXMDM12*(GF1 - 2. + XD*X) - 1./2.*X*XD*XLNX2 
     > - 1./2.*GF1*(XLNXMDM2 + XLNXM12) - 1./2.*XLNXMDM1*XDM1*(3.*D  
     > + X*D*(3. - 1./D)/DM1) - 2.*XLNXM1*(X*XMDM1  
     > + 3./2.*(1. + XMDM*XM1)) + 1./2.*XLNDM1*(5.*X*XMDM1  
     > + 3.*XMDM*(1. + DM1)) + XLNX*X*(1.+ X*(1.+ 1./2.*1./D)/DM1) 
     > - 3./2.*XMDM1*(1. + 3.*X) - 3./2.*X*XMDM 
 
	CFNFTF = GF1*(XLNXMDM + XLNXM1) + 10./3.*GF1 - 2.*GF1A + 2.*X   
     > + 2.*XMDM1 
 
	CGCF = 1./2.*XLNXMDM12*(XDM1**2*(4. - 3.*D + 8.*XD) - 2. + X*XD) 
     > + 1./2.*GF1*(XLNDMXM12 + XLNXMX2) + 1./2.*XLNX2*XD2*(3./DM1 + DM1 
     > + 8.*XDM1 - GF1/XD2) + (XLNXMDM + XLNXMX)* 
     > (17./6.*GF1 - XMDM*(2. - D))  
     > + XLNXMDM1*(XDM1*(2.*X + 3.*D - 5.*XD + 8.*X*DM1*XDM1/D) 
     > + XD*XMDM - XD2) + XLNX*(11./6. - 11./3.*X*XD  
     > + (11. - 18.*X - 11.*X2)/(6.*DM1) - XD*XDM1*(4. + 8./3.*D  
     > - 8.*XD)/DM1) + 4.*XD*XMDM*(2. + X - 2.*XD) + 20./3.*X  
     > + 16./9. - XMDM*(7. - 118./9.*XM1) - PI2/3.*GF1 + 11./3.*XMDM1   
     > - 1./2.*HGF - 1./2.*HGFB 
 
	GF2D = CF**2*CF2 - 2./3.*CF*TFNF*CFNFTF + CG*CF*CGCF  

	RETURN 
 
C	---------------------------- 
 
	ENTRY GG2D (N,XX,YY) 
 
	TFNF = TF*N 
 
	BO = 4.*TFNF/3. - 11.*CG/3. 
 
	X = XX 
 
	D = YY 
 
	XD = X/D 
 
	DX = D/X 
 
	XM1 = 1.- X 
 
	DM1 = 1.- D 
 
	DMXM1 = DM1/XM1 - 1. 
 
	XMDM = XM1/DM1 
 
	X2 = X**2 
 
	D2 = D**2 
 
	XD2 = X2/D2 
 
	XODM1 = X/DM1 
 
	XDM1 = 1. - XD 
 
	XMDM1 = 1. - XMDM 
 
	XLNDMXM1 = LOG(DMXM1) 
 
	XLND = LOG(D) 
 
	XLNX = LOG(X) 
 
	XLNXMX = LOG(XM1/X) 
 
	XLNXD = LOG(XD) 
 
	XLNXM1 = LOG(XM1) 
 
	XLNDM1 = LOG(DM1) 
 
	XLNDXM1 = LOG(1.- DX) 
 
	XLNXMDM = LOG(XMDM) 
 
	XLNXMDM1 = LOG(XMDM1) 
 
	XLND2 = XLND**2 
 
	XLNX2 = XLNX**2 
 
	XLNXD2 = XLNXD**2 
 
	XLNXM12 = XLNXM1**2 
 
	XLNDM12 = XLNDM1**2 
 
	XLNDXM12 = XLNDXM1**2 
 
	XLNXMDM12 = XLNXMDM1**2 
 
	XLNXMDM2 = XLNXMDM**2 
 
	XLNDMXM12 = XLNDMXM1**2 
 
	XLNXMX2 = XLNXMX**2 
 
	GG1 = GG1F(X,D) 
 
	GG1C = GGC1F(X,D) 
 
	GG1A = GGA1F(X,D) 
 
	FG1 = FG1FA(X,D) 
 
	HGG = HGGD(X,D) 
 
	HGGB = HGGDB(X,D) 
 
	CGNFTF = XLNXMDM1*XDM1**2*D - XD*X*XLNX/DM1**2 - 5./3.*GG1  
     > - 13./3.*GG1C - 9./2.*GG1A - (D/DM1)**2*XM1 
 
	CFNFTF = XLNXMDM12*XDM1**2*D -XD*X*XLNX2/DM1**2  
     > + 2.*XLNX*(XMDM*X/DM1 - 2.* XD*X/DM1**2)  
     > + 2.*XLNXMDM1*((D - 2.*X)*XDM1 - DM1*XMDM1**2) 
     > - 20./3.*GG1C - 12.*GG1A - 4.*FG1 + 8.*XMDM1/DM1  
     > - 4.*XM1*(D/DM1)**2 
 
	CA2 = XLNXMDM12*(1. - 1./XM1 - 2.*D + 5.*X - 2.*XD2*((2. + D) 
     > *XDM1 + 1. + D)) + GG1*(XLNXMDM1*XLNDMXM1 - XLNXM1*XLNX)  
     > + XLNX2*(GG1 - 1./XM1 - (1. - 3.*X)/DM1**2  
     > + 2.*XD2*((3. - 1./DM1)*XDM1 + 1.)/DM1) 
     > + XLNXMDM1*(4.*XD2*(2.*XD*DM1 - XMDM - 2.)  
     > + XD*(4. + 55./6.*X) 
     > - 7./3.*X - 5./6.*D + 2.*D*XMDM) + XLNX*(2.*XD*XDM1 
     > *(2. - 4.*XDM1- D2) 
     > + XD2*DM1*(4.*X + 5./6.*D2/DM1 - 31./6.*D))/DM1**2 
     > + (67./18. - PI2/3.)*GG1 + 67./9.*GG1C + 11./2.*GG1A  
     > + 8.*XD*(1.- X2 - XM1*XD)/DM1**2 
     >  - 4.*XMDM*(2. - D*XM1)/DM1 + 5./6.*(D/DM1)**2*XM1 
     > - 1./2.*HGG - 1./2.*HGGB 
 
	GG2D = 2./3.*CG*TFNF*CGNFTF + CF*TFNF*CFNFTF + CG**2*CA2 
  
	RETURN 
 
C	--------------------  
 
	END 

C
C     Function to determine the value of the DGLAP kernels at y = Delta/x
C

       FUNCTION PNSMCD (N,X)

       IMPLICIT DOUBLE PRECISION (A-H, O-Z)

       PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)

       EXTERNAL SPENCE

       XM1 = 1. - X

       SPENX = SPENCE(X)

       XLNX = LOG(X)

       XLNX2 = XLNX**2

       XLNXM1 = LOG(XM1)

       XN = 1D0*N

       PNSMCD = X*(- 4./9.*SPENX/XM1 + (4./3. + 10./9.*X/XM1)*XLNX2
     >- 20./9.*XLNX*XLNXM1/XM1 + 34./9.*XLNX/XM1 - 16./27.*PI2/XM1
     >+ 28. + 134./9.*X/XM1 - XN*(32./27. + 20./27.*X/XM1 
     >+ 4./9.*XLNX/XM1))

       RETURN
C     ------------------------------------------------------

       ENTRY PNSPCD (N,X)

        XM1 = 1. - X

       SPENX = SPENCE(X)

       XLNX = LOG(X)

       XLNX2 = XLNX**2

       XLNXM1 = LOG(XM1)

       XN = 1D0*N

       PNSPCD = X*(- 4./9.*SPENX/XM1 + (8./9. + 10./9.*X/XM1)*XLNX2
     >- 20./9.*XLNX*XLNXM1/XM1 + 34./9.*XLNX/XM1 - 16./27.*PI2/XM1
     >+ 28. + 134./9.*X/XM1 - XN*(32./27. + 20./27.*X/XM1 
     >+ 4./9.*XLNX/XM1))

       RETURN
C     ------------------------------------------------------

       ENTRY GG2CDA (N,X)

       XM1 = 1. - X

       XM12 = XM1**2

       SPENX = SPENCE(X)

       XLNX = LOG(X)

       XLNX2 = XLNX**2

       XLNXM1 = LOG(XM1)

       XN = 1D0*N
       
       GG2CDA = 9.*PI2 - 9.*PI2/XM1 + 3.*PI2*X*(1. + X)/XM1
     >- 126. + 126./XM1 - 85./2.*X*(1. + X)/XM1 + 27.*(1. -
     >1./XM1)*SPENX + 9.*SPENX*X*(1. + X)/XM1 - 18.*X*(1. + 
     >1./XM1)*XLNX*XLNXM1 + X*XLNX2*(27./2. + 18./XM1 
     >+ 45./2./XM12) + X*XLNX*(18. + 15./2./XM1)/XM1 + XN*X*(
     >- 2./3.*XLNX2/XM12 - (4./3. + 7./3./XM1)*XLNX/XM1 - 9./2.
     >- 4./XM1)

       RETURN
C     ---------------------------------------------------------

       ENTRY GG2CD (N,X)

       XM1 = 1. - X

       XM12 = XM1**2

       SPENX = SPENCE(X)

       XLNX = LOG(X)

       XLNX2 = XLNX**2

       XLNXM1 = LOG(XM1)

       XN = 1D0*N

       GG2CD = - 9.*PI2 + 3.*PI2/XM1 - 3.*PI2*X*(1. + X)/XM1
     >+ 66. - 53./2./XM1 + 45./2.*(1. - X**3)/XM12 + 27.*SPENX
     >- 9.*SPENX/XM1 + 9.*SPENX*X*(1. + X)/XM1 + 27./2.*X*XLNX2
     >/XM12 + 18.*XLNX - 51./2.*XLNX/XM1 + 15./2.*XLNX/XM12 
     > - 9./2.*X*XLNX2
     > + XN*(- 29./18.*XM1 - 1./6.*(1. + X)/XM1 - 5./XM1 - 1./3.*
     >X*XLNX*(7. + 4.*X + 2.*XLNX)/XM12)

       RETURN

C     ---------------------------------------------------------

       ENTRY PNSML(N,z)

       tem1 = -2./27.*(-390.*z**3+24.*z**2*log(1.-z)+6.*SPENCE(z)
     > +6.*log(1.-z)
     > +24.*log(1.-z)**2*z**2-26.*z**2*N+16.*z**3*N-30.*log(1.-z)*z
     > +591.*z**2
     > -12.*z**2*N*log(z)+6.*z**3*N*log(z)+48.*log(z)*z**3*log(1.-z)
     > -96.*log(z)*z**2*log(1.-z)-12.*log(1.-z)**2*z**3
     > +48.*log(z)**2*z**2
     > -63.*log(z)*z**3-24.*log(z)**2*z**3+8.*z**3*PI2-16.*z**2*PI2
     > +150.*log(z)*z**2+12.*z**3*SPENCE(z)-6.*SPENCE(z)*z
     > -18.*z**2*SPENCE(z)-12.*log(1.-z)**2*z)/(z*(-1.+z))

       PNSML = tem1

       RETURN

C     ---------------------------------------------------------

       ENTRY PNSPL(N,z)

       tem1 = -2./27.*(-390.*z**3-12.*z**2*log(1.-z)-6.*SPENCE(z)
     > -6.*log(1.-z)
     > +24.*log(1.-z)**2*z**2-26.*z**2*N+16.*z**3*N+6.*log(1.-z)*z
     > +591.*z**2
     > +12.*log(1.-z)*z**3-12.*z**2*N*log(z)+6.*z**3*N*log(z)
     > +48.*log(z)*z**3*log(1.-z)-96.*log(z)*z**2*log(1.-z)
     > -12.*log(1.-z)**2*z**3
     > +48.*log(z)**2*z**2-63.*log(z)*z**3-24.*log(z)**2*z**3
     > +8.*z**3*PI2
     > -16.*z**2*PI2+150.*log(z)*z**2+6.*SPENCE(z)*z-6.*z**2*SPENCE(z)
     > -12.*log(1.-z)**2*z)/(z*(-1.+z))

       PNSPL = tem1

       RETURN

C     ---------------------------------------------------------

       ENTRY SFF2L(N,z)

       IF (z.EQ.0D0) THEN

          tem1 = 0.0

          ELSE

       tem1 = 4./3.*N*(-2./3.*z*(-10.+19.*z)*LOG(z)+z*(5.-8.*z)
     > -((1.-z)*(-3.+8.*z)+2./3.*(1-z)*(9.-19.*z))*LOG(1.-z)
     > +2*LOG(z)*z-(1.-z-4.*(1-z)*z)*LOG(1.-z)**2+8.*(1.-z)*z
     > *SPENCE(1.-z) -580./9.*(1.-z)*z+3.*z+8.*(1-z)*z*LOG(1.-z)
     > *LOG(z))

       ENDIF

       SFF2L = tem1

       RETURN
C     ---------------------------------------------------------

       ENTRY SFF2LA(N,z)

       IF (z.EQ.0D0) THEN

          tem1 = -4./3.*N

          ELSE

       tem1 = -4./3.*N*(-LOG(z)*z+z*(1-LOG(z))-(1.-z)*LOG(1.-z)
     >        *(1.-LOG(1.-z)))

       ENDIF

       SFF2LA = tem1

       RETURN

C     ---------------------------------------------------------

       ENTRY FG2L(N,z)

       tem1 = -1./18.*(3165.*z**3-1692.*z**4+663.*z**2*log(1.-z)
     > -27.*log(1.-z)
     > +327.*log(1.-z)**2*z**2-27.*z-294.*log(1.-z)*z-1446.*z**2
     > -411.*log(1.-z)*z**3
     > +288.*log(z)*z**4*log(1.-z)-528.*log(z)*z**3*log(1.-z)
     > +240.*log(z)*z**2*log(1.-z)
     > -426.*log(1.-z)**2*z**3+48.*z**4*PI2+180.*log(1.-z)**2*z**4
     > -120.*log(z)**2*z**2
     > -405.*log(z)*z**3+264.*log(z)**2*z**3-144.*log(z)**2*z**4
     > -88.*z**3*PI2
     > +40.*z**2*PI2+15.*log(z)*z**2+54.*log(1.-z)*z**4
     > +1296.*z**3*SPENCE(z)
     > -648.*z**4*SPENCE(z)+54.*SPENCE(z)*z+378.*z**4*log(z)
     > -702.*z**2*SPENCE(z)-81.*log(1.-z)**2*z)*N/(z*(-1.+z))

       FG2L = tem1

       RETURN
C     ---------------------------------------------------------

       ENTRY FG2LA(N,z)

       tem1 = -1./18.*(159.*z**3+27.*log(1.-z)+27.*z+54.*z*SPENCE(z)
     > +48.*z**2*log(z)*log(1.-z)+177.*log(1.-z)**2*z**2
     > -24.*log(z)**2*z**2
     > -162.*z**2*SPENCE(z)-81.*log(1.-z)**2*z-48.*log(1.-z)*z
     > -15.*log(1.-z)*z**2
     > +81.*z**2*log(z)-186.*z**2+8.*z**2*PI2-96.*z**3*log(1.-z)**2
     > +108.*z**3*SPENCE(z)-48.*z**3*log(1.-z)*log(z)
     > +24.*z**3*log(z)**2
     > +21.*z**3*log(1.-z)-8.*z**3*PI2-93.*log(z)*z**3)*N/((-1.+z)*z)

       FG2LA = tem1

       RETURN

C     ---------------------------------------------------------

       ENTRY GF2L(N,z)


       tem1 = -1264./9.*z**3+4.*log(1.-z)-10./3.*z-4.*SPENCE(z)
     > +160./27.*N*z**3
     > -8.*z*SPENCE(z)+152./3.*log(1.-z)*z**3+152./3.*log(z)*z**3
     > +4./9.*N*z**2*log(z)
     > +4./9.*N*log(1.-z)*z**2-4./3.*N*z-164./27.*N*z**2
     > -8./9.*N*log(1.-z)
     > -40.*z**2*log(z)*log(1.-z)-64./9.*log(1.-z)**2*z**2
     > +20.*log(z)**2*z**2
     > +48.*z**2*SPENCE(z)+12.*log(1.-z)**2*z+592./9.*log(1.-z)*z
     > -1042./9.*log(1.-z)*z**2-146./3.*z**2*log(z)+170.*z**2
     > -20./3.*PI2*z**2
     > +16./3.*PI2*z**3+32.*log(z)*z**3*log(1.-z)+4.*z*log(z)
     > -34./9.*log(1.-z)**2-32.*z**3*SPENCE(z)-16.*log(z)**2*z**3

       GF2L = tem1

       RETURN

C     ---------------------------------------------------------

       ENTRY GF2LA(N,z)


       tem1 = -64./9.*log(1.-z)-130./9.*z+4.*SPENCE(z)-8.*z*SPENCE(z)
     > -4./9.*N*z**2*log(z)-4./9.*N*log(1.-z)*z**2-4./9.*N*z
     > -4./27.*N*z**2
     > -8.*z**2*log(z)*log(1.-z)-44./9.*log(1.-z)**2*z**2
     > +4.*log(z)**2*z**2
     > +8.*z**2*SPENCE(z)+12.*log(1.-z)**2*z-16./9.*log(1.-z)*z
     > +34./3.*log(1.-z)*z**2
     > -34./3.*z**2*log(z)+278./9.*z**2-4./3.*z**2*PI2+4.*z*log(z)
     > -6.*log(1.-z)**2

       GF2LA = tem1

       RETURN
       
C     ---------------------------------------------------------

       ENTRY GG2L(N,z)

       tem1 = 2./3.*N*(20.*z**3+z**2*LOG(1.-z)**2+10.*LOG(1.-z)*z**2
     > -3.*LOG(z)*z**2-32.*z**2+11.*z-2.*LOG(1.-z)**2*z-18.*LOG(1.-z)*z
     > +8.*LOG(1.-z)+LOG(1.-z)**2)

       tem2 = 1./3.*N*(99.*z**2-3.*z+3.*LOG(1.-z)*z**3-9.*LOG(1.-z)*
     > z**2
     > +9.*LOG(1.-z)*z-3.*LOG(1.-z)-160.*z**3+69.*z**4)/(-1.+z)

       tem3 = 1./4.*(108.*PI2*z**2+336.*z+18.-3411.*z**2+72.*PI2*z**4
     > -168.*PI2*z**3-2112.*z**4+5035.*z**3)/(-1.+z)
     > -9./2.*z*LOG(z)*(z+2.-28.*LOG(z)*z**2+12.*LOG(z)*z**3+18.
     > *LOG(z)*z
     > -14.*z**3+14.*z**2)/(-1.+z)
     > + 3./2.*LOG(1.-z)*((-168.*z**3*LOG(z)+72.*z**4*LOG(z)+108.
     > *LOG(z)*z**2
     > +354.*z**2-180.*z+20.+78.*z**4-281.*z**3)/(-1.+z)+9/((-1.+z)*z))
     > -18.*(z**2-2.*z+1.)*LOG(1.-z)**2  
     > -9.*SPENCE(z)*(-1./((-1+z)*z) + (-30.*z**3+12.*z**4+22.*z**2
     > -3.*z+1.)/(-1.+z))

       GG2L = tem1 + tem2 + tem3

       RETURN
C     ---------------------------------------------------------

       ENTRY GG2LA(N,z)

       tem1 = 1./12.*(144.*z**3*N-2961.*z**3-162.*log(1.-z)-54.*z
     > -108.*SPENCE(z)
     > +1713.*z**4+108.*z*SPENCE(z)-8.*N*log(1.-z)**2*z
     > +24.*N*log(1.-z)**2*z**2
     > -28.*N*log(1.-z)*z+68.*N*log(1.-z)*z**2-36.*N*z**2
     > -648.*log(1.-z)**2*z**2
     > +324.*z**2*SPENCE(z)+216.*log(1.-z)**2*z+900.*log(1.-z)*z
     > -864.*log(1.-z)*z**2
     > -108.*z**2*log(z)+900.*z**2+648.*z**3*log(1.-z)**2
     > -864.*z**3*SPENCE(z)
     > +648.*z**3*log(1.-z)*log(z)-324.*z**3*log(z)**2
     > -432.*z**4*log(1.-z)*log(z)
     > +12.*z**4*N*log(1.-z)+8.*N*log(1.-z)**2*z**4
     > -24.*N*log(1.-z)**2*z**3-88.*N*z**4
     > -8.*N*z**3*log(z)-52.*N*z**3*log(1.-z)+8.*N*z**4*log(z)
     > -72.*z**4*PI2
     > +432.*z**4*SPENCE(z)-324.*z**4*log(z)-216.*z**4*log(1.-z)**2
     > +234.*z**4*log(1.-z)
     > +216.*z**4*log(z)**2-108.*z**3*log(1.-z)+108.*z**3*PI2
     > +270.*log(z)*z**3)/((-1.+z)*z)

       GG2LA = tem1

       RETURN

C     ---------------------------------------------------------

       ENTRY SFF2BD(N,z1)

       TEM1 =  4./3.*N*(5.-13.*z1+8.*z1**2-log(z1)**2*z1
     >         -3.*log(z1)*z1)/(z1*(-1+z1))

       SFF2BD = TEM1 

       RETURN
C     ---------------------------------------------------------

       ENTRY FG2BD(N,z1)

       TEM1 = -1./3.*(24.-6.*log(z1-1.)*log(z1)-48.*z1**3+6.*log(z1)**2
     > *z1**2-PI2-6.*log(z1-1.)*log(z1)*z1**2+12.*log(z1-1.)*log(z1)
     > *z1-78.*z1+102.*z1**2+3.*log(z1-1.)**2-6.*log(z1-1.)**2*z1
     > +2.*PI2*z1+3.*log(z1-1.)**2*z1**2-PI2*z1**2-6.*SPENCE(1.-z1)
     > +12.*SPENCE(1.-z1)*z1-6.*log(z1)*z1-6.*SPENCE(z1)*z1**2
     > +18.*log(z1)*z1**2-6.*log(z1-1.)+6.*log(z1-1.)*z1)
     > /(z1**2*(z1-1.)**2)

       TEM2 = -1./3.*(-42.+2.*PI2-39.*z1**2+2.*PI2*z1**2-4.*PI2*z1
     > +81.*z1+3.*log(z1)*z1**2-3.*log(z1-1.)**2*z1**2
     > +6.*log(z1-1.)**2*z1-6.*log(z1-1.)*z1-3.*log(z1-1.)**2
     > +6.*log(z1-1.))/(z1**2*(z1-1.)**2)

       FG2BD = CG*TF*N*TEM1 + N*TF*CF*TEM2

       RETURN
C     ---------------------------------------------------------

       ENTRY GF2BD(N,z1)

       TEM1 = -1./2.*(-3.-3.*log(z1)+log(z1-1.)**2+3.*log(z1-1.)
     > -9.*z1*log(z1-1.)-4.*z1**2*log(z1-1.)*log(z1)
     > -2.*log(z1-1.)*log(z1)+2.*z1**2*log(z1-1.)**2
     > +2.*z1**2*log(z1)**2-3.*z1*log(z1-1.)**2
     > +6.*z1*log(z1-1.)*log(z1)+12.*log(z1)*z1-2.*log(z1)**2*z1
     > +6*z1**2*log(z1-1.)-6.*z1**2*log(z1)+3.*z1)/(z1*(z1-1.))

       TEM2 = 1./18.*(44.-66.*z1**2*log(z1)-3.*PI2+66.*z1**2*log(z1-1.)
     > +18.*log(z1)**2*z1+9.*log(z1-1.)**2-33.*log(z1)
     > +36.*SPENCE(1.-z1)*z1**2+153.*log(z1)*z1-54.*SPENCE(1.-z1)*z1
     > -6.*PI2*z1**2+18.*z1**2*log(z1-1.)**2+9.*PI2*z1
     > -27.*z1*log(z1-1.)**2-99.*z1*log(z1-1.)+33.*log(z1-1.)
     > +18.*SPENCE(1.-z1)-186.*z1+142.*z1**2)/(z1*(z1-1.))

       TEM3 = 2./9.*(-6.*z1*log(z1-1.)+6.*log(z1)*z1-20.*z1
     > +3.*log(z1-1.)-3.*log(z1)+10.)/z1

       GF2BD = CF**2*TEM1 + CF*CG*TEM2 + CF*TF*N*TEM3

       RETURN
C     ---------------------------------------------------------

       ENTRY PNSPBD (N,z1)

       TEM1 = 1./2.*(-2.+4.*log(z1)*z1+log(z1)**2*z1-2.*log(z1)**2
     > +4.*SPENCE(1.-z1)*z1+2.*z1+2.*log(z1)*z1*log(-1.+z1))
     > /(z1*(-1.+z1))

       TEM2 = 2./9.*(3.-8.*z1+3.*log(z1)*z1)/(z1*(-1.+z1))

       TEM3 = -1./18.*(51.-9.*log(z1)**2+3*PI2*z1+33.*log(z1)*z1
     > +18.*SPENCE(z1)*z1-118.*z1+9.*log(z1)**2*z1)/(z1*(-1.+z1))

       PNSPBD = CF**2*TEM1 + CF*N*TF*TEM2 + CF*CG*TEM3

       RETURN
C     ---------------------------------------------------------

       ENTRY PNSMBD (N,z1)

       TEM1 = 1./2.*(-2.+4.*log(z1)*z1+log(z1)**2*z1-2.*log(z1)**2
     > +4.*SPENCE(1.-z1)*z1+2.*z1+2.*log(z1)*z1*log(-1.+z1))
     > /(z1*(-1.+z1)) -2.*log(z1)**2/z1

       TEM2 = 2./9.*(3.-8.*z1+3.*log(z1)*z1)/(z1*(-1.+z1))

       TEM3 = -1./18.*(51.-9.*log(z1)**2+3*PI2*z1+33.*log(z1)*z1
     > +18.*SPENCE(z1)*z1-118.*z1+9.*log(z1)**2*z1)/(z1*(-1.+z1))
     > + log(z1)**2/z1

       PNSPBD = CF**2*TEM1 + CF*N*TF*TEM2 + CF*CG*TEM3

       RETURN
C     ---------------------------------------------------------

       ENTRY GG2BD(N,z1)

       TEM1 = -1./3.*(4.+3.*log(z1)**2*z1**2-6.*log(z1)*z1
     > -6.*log(z1)*z1**2-24.*z1+24.*z1**2-4.*z1**3)
     > /(z1**2*(z1-1.)**2)

       TEM2 = 1./9.*(-60.*z1+15.-46.*z1**3+91.*z1**2+6.*log(z1)*z1**2)
     > /(z1**2*(z1-1.)**2)

       TEM3 = -1./18.*(18.*SPENCE(1.-z1)-3.*PI2+72.*log(z1)*z1**2*
     > log(z1-1.)-36.*log(z1)*z1**3*log(z1-1.)
     > -54.*log(z1)*z1*log(z1-1.)+45+18*log(z1)*log(z1-1.)+6.*PI2*z1**3
     > +9.*PI2*z1-54.*SPENCE(1.-z1)*z1+36.*log(z1)*z1+9.*z1*log(z1)**2
     > -21.*log(z1)*z1**2+18.*z1**3*log(z1)**2-54.*log(z1)**2*z1**2
     > +72.*SPENCE(1.-z1)*z1**2-36.*SPENCE(1.-z1)*z1**3-12.*PI2*z1**2
     > -132.*z1+211.*z1**2-124.*z1**3)/(z1**2*(z1-1.)**2)

       GG2BD = CF*N*TF*TEM1 + CG*TF*N*TEM2 + CG**2*TEM3

       RETURN
C     ---------------------------------------------------------

        ENTRY SFF2BDA(N,z1)

        SFF2BDA = -2.*CF*N*TF*log(z1)*(1.+log(z1))/(-1.+z1)

         RETURN
C     ---------------------------------------------------------

        ENTRY FG2BDA(N,z1)

        TEM1 = -1./3.*(30.-9.*log(z1)*z1**2+45.*z1**2-75.*z1+4.*PI2*z1
     >  -2.*PI2*z1**2-2.*PI2-12.*log(z1-1.)*z1**2+18.*log(z1-1.)*z1
     >  -6.*log(z1-1.)-6.*log(z1-1.)**2*z1+3.*log(z1-1.)**2*z1**2
     >  +3.*log(z1-1.)**2)/(z1**2*(z1-1.)**2)

        TEM2 = -1./3.*(-12.+6.*log(z1)**2*z1**2+6.*log(z1)*z1**2
     > -30.*z1**2+42.*z1+6.*SPENCE(1.-z1)-12.*z1*log(z1-1.)*log(z1)
     > +6.*z1**2*log(z1-1.)*log(z1)-3.*log(z1-1.)**2+6.*log(z1-1.)
     > +PI2+6.*SPENCE(1.-z1)*z1**2-12.*SPENCE(1.-z1)*z1+6.*log(z1)*z1
     > +PI2*z1**2-2.*PI2*z1-3.*log(z1-1.)**2*z1**2-18.*log(z1-1.)*z1
     > +12.*log(z1-1.)*z1**2+6.*log(z1-1.)*log(z1)+6.*log(z1-1.)**2*z1)
     > /(z1**2*(z1-1.)**2)

        FG2BDA = CF*N*TF*TEM1 + CG*N*TF*TEM2

         RETURN
C     ---------------------------------------------------------

       ENTRY GF2BDA(N,z1)

       TEM1 =-1./2.*(-6.-2.*z1*log(-1.+z1)*log(z1)+6.*z1+3.*log(z1)
     > -4.*log(z1)*z1+z1*log(-1.+z1)**2+2.*log(-1.+z1)*log(z1)
     > +z1*log(-1.+z1)-log(-1.+z1)-log(-1.+z1)**2)/((-1.+z1)*z1)
          
       TEM2 = -2./9.*(3.*log(-1.+z1)-3.*log(z1)+4.)/z1

       TEM3 = -1./18.*(104.-33.*log(z1)+9.*log(-1.+z1)**2
     > +18.*SPENCE(1.-z1)-3.*PI2-104.*z1+15.*log(-1.+z1)
     > +33.*log(z1)*z1+3.*PI2*z1-9.*z1*log(-1.+z1)**2
     > -18.*SPENCE(1.-z1)*z1-18.*log(z1)**2*z1-15.*z1*log(-1.+z1))
     > /((-1.+z1)*z1)

       GF2BDA = CF**2*TEM1 + CF*N*TF*TEM2 + CF*CG*TEM3

       RETURN
C     ---------------------------------------------------------

       ENTRY GG2BDA(N,z1)

       TEM1 = -(2.+4.*z1**2+2.*log(z1)*z1-6.*z1+log(z1)**2*z1**2
     > -4.*log(z1)*z1**2)/(z1**2*(z1-1.)**2)

       TEM2 = 1./9.*(54.*z1-19.-35.*z1**2+6.*log(z1)*z1**2)
     > /(z1**2*(z1-1.)**2)

       TEM3 = -1./18.*(-85.+3.*PI2-18.*SPENCE(1.-z1)-36.*log(z1)*z1
     > -167.*z1**2+6.*PI2*z1**2-9.*PI2*z1-36.*z1**2*SPENCE(1.-z1)
     > +54.*z1*SPENCE(1.-z1)-18.*log(z1)**2*z1**2+51.*log(z1)*z1**2
     > -9.*log(z1)**2*z1+252.*z1-18.*log(z1)*log(z1-1.)
     > +54.*log(z1)*log(z1-1.)*z1-36.*log(z1)*log(z1-1.)*z1**2)
     > /(z1**2*(z1-1.)**2)

       GG2BDA = CF*N*TF*TEM1 + CG*TF*N*TEM2 + CG**2*TEM3

        RETURN
C     ---------------------------------------------------------

C
       END
C
       SUBROUTINE KERNEL 
     >(FF1, GG1, PNSPBLK, PNSPBLC1, PNSPBLCX, PNSMBLK, PNSMBLC1,  
     > PNSMBLCX, SFF2BLK, GG2BLK, GG2BLC1,GG2BLCX, FG2BLK, GF2BLK, 
     > PNSP, PNSM, SFF2, FG2, GF2, GG2, FF1BLX,GG1BLX,FF1BL1,
     > GG1BL1,PNSMBLKD,PNSPBLKD,SFF2BLKD,GG2BLKD,FG2BLKD,GF2BLKD,
     > FG2BLK1,GF2BLK1,SFF2BLK1,PNSPC,PNSMC,GG2C, PNSPDGBL,PNSMDGBL,
     > GG2DGBL, NFL,IRT) 
 
C     New version with 'regularized' kernel functions which are smooth, hence 
C     are more suitable for interpolation. 
C 
C     Subroutine to calculate the values of the 1st and 2nd order evolution 
C                kernel function at a given value of X. 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      PARAMETER (D0 = 0.0, D1 = 1.0) 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6)	 
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)	 
 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX) 
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX, DEL, IV	 
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx 
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX)
 
      EXTERNAL FF1F, FF1BLF, GG1F, GG1BLF, GG1BLFA, GG1FA, PNSPBL, 
     > SFF2BLA, PNSPD, PNSMD, SFF2D, SFF2DA, FG2D, FG2DA, GF2D, GF2DA,  
     > GG2D, GG2DA, PNSMBL, SFF2BL, FG1BLF, GF1BLF, FG1BLFA, GF1BLFA,
     > PNSPCD, PNSMCD, GG2CDA, GG2CD,GG2L,SFF2L,SFF2LA,SFF2BD,FG2BD,
     > GF2BD,GG2BD,PNSPBD,PNSMBD,SFF2BDA,FG2BDA,GF2BDA,GG2BDA,GG2BL,
     > FG2BL,GG2BLA,FG2BLA,GF2BL,GF2BLA,PNSPL,PNSML,FG2L,FG2LA,GF2L,
     > GF2LA,GG2LA 
 
      DIMENSION FF1(0:MXX,0:MXX), GG1(0:MXX,0:MXX),    
     > FG2(0:MXX,0:MXX),GF2(0:MXX,0:MXX),GG2(0:MXX,0:MXX), 
     > PNSP(0:MXX,0:MXX), PNSM(0:MXX,0:MXX), SFF2(0:MXX,0:MXX), 
     > SFF2BLK(0:MXX,0:MXX), FG2BLK(0:MXX,0:MXX), GF2BLK(0:MXX,0:MXX),  
     > GG2BLK(0:MXX,0:MXX), PNSPBLK(0:MXX,0:MXX), PNSMBLK(0:MXX,0:MXX), 
     > GG2BLC1(0:MXX,0:MXX), PNSPBLC1(0:MXX,0:MXX),FF1BLX(0:MXX,0:MXX),  
     > GG2BLCX(0:MXX,0:MXX), PNSPBLCX(0:MXX,0:MXX),GG1BLX(0:MXX,0:MXX),  
     > PNSMBLC1(0:MXX,0:MXX), PNSMBLCX(0:MXX,0:MXX),FF1BL1(0:MXX,0:MXX), 
     > GG1BL1(0:MXX,0:MXX), SFF2BLKD(0:MXX,0:MXX), FG2BLKD(0:MXX,0:MXX), 
     > GF2BLKD(0:MXX,0:MXX), GG2BLKD(0:MXX,0:MXX),PNSPBLKD(0:MXX,0:MXX), 
     > PNSMBLKD(0:MXX,0:MXX),SFF2BLK1(0:MXX,0:MXX),FG2BLK1(0:MXX,0:MXX), 
     > GF2BLK1(0:MXX,0:MXX),PNSPC(0:MXX,0:MXX),PNSMC(0:MXX,0:MXX),
     > GG2C(0:MXX,0:MXX), PNSPDGBL(0:MXX,0:MXX), PNSMDGBL(0:MXX,0:MXX),
     > GG2DGBL(0:MXX,0:MXX)
      DIMENSION T(10),T1(10),T2(10),XX2(10)


C
      DATA IWRN / 0 /
      DATA CF, CG, TF, TINY / 1.333333333333, 3.0, 0.5, 1E-10/	
C
      IRT = 0
C
C     scaling for subtraction integrals
C

      IF (DEL.GE.0.001 .AND. DEL.LT.0.03) THEN
	ET2 = 1000.
	ET1 = 100.
        ET3 = 1000.
        ELSEIF (DEL.GE.0.03 .AND. DEL.LT.0.1) THEN
        ET2 = 100.
	ET1 = 50.
        ET3 = 100.  
        ELSEIF (DEL.GE.0.1) THEN
        ET2 = 100.
	ET1 = 50.
        ET3 = 100.  
        ELSEIF (DEL.LT.0.001) THEN
           ET2 = 1000.
           ET1 = 100.
           ET3 = 1000.
        ENDIF

      DIM = 1D0/ZV(IV)
      DIM1 = 1D0/DEL


C
C	Set up the kernel arrays for 2nd. order kernels plus subtraction and 
C	counterterm kernels. First unpolarized then polarized case. 
C

	IF (IKNL .EQ. 2) THEN

	DO 100 IX = 2, NX-1

        IF (IX.EQ.IV) GOTO 100

        X = XV(IX)
	XD = X*DIM1
        DX = 1./XD
	XMD = 1. - XD
        X1 = XV(IX)
	XD1 = X1*DIM1
        DX1 = 1./XD1
	XMD1 = 1. - XD1
	IN = IX - IV

	IF (IX.LE.IV) THEN

	IQ = 1

	IQ1 = NX

	ELSE

	IQ = IV+1

        IQ1 = NX

	ENDIF 

        K1 = 1

	DO 101 IY = IQ, IQ1 

           Y = XV(IY)
           Y1 = XV(IY)
           
        IF (IQ.EQ.IV+1.AND.IY.LT.IX) THEN

           GOTO 103

        ELSE

	IF (IY .EQ. IX .AND. IY.GT.IV) THEN

           SM = 0.0000001
           X11 = 0.9999999

           FG2(IN,IN) = 0.0

           GF2(IN,IN) = 0.0

           SFF2(IN,IN) = 0.0

           PNSM(IN,IN) = - (PNSMD(NFL,X1,DEL)/X1)/ET2

           PNSP(IN,IN) = - (PNSPD(NFL,X1,DEL)/X1)/ET2

           GG2(IN,IN) = - GG2D(NFL,X1,DEL)/ET1

	GOTO 103

	ELSE

	XY = X/Y
	YD = Y*DIM1
	DY = 1./YD
	YMD = 1. - YD
	XOY = XMD/YMD
	DYM1 = 1. - DY 
	XYM1 = 1. - XY

        XY1 = X1/XV(IY)
	YD1 = XV(IY)*DIM1
	DY1 = DEL/XV(IY)
	YMD1 = 1. - YD1
	XOY1 = XMD1/YMD1
	DYM11 = 1. - DY1 
	XYM11 = 1. - XY1
	
C
C	BL Region for the integral from x to 1! This integral is split in two regions  
C       x..DEL and DEL..1. In the first region there will be a regularization of the  
C       kernel, the second region is singularity free! IY + 1 in order to extrapolate  
C       down to DEL at position IY for the integral from DEL..1! 
C
	
	IF (IX.LE.IV.AND.IY.GT.IV) THEN	
           
	PNSMBLK (IX,IY) = PNSMBL(NFL,XD,YD)

	PNSPBLK (IX,IY) = PNSPBL(NFL,XD,YD)

	GG2BLK (IX,IY) = GG2BL(NFL,XD,YD)

	FG2BLK1 (IX,IY) = FG2BL(NFL,XD,YD)

	GF2BLK1 (IX,IY) = GF2BL(NFL,XD,YD)

	SFF2BLK1 (IX,IY) = SFF2BL(NFL,XD,YD)

        PNSMBLKD (IX,IY+1) = PNSMBL(NFL,XMD,YD)

        PNSPBLKD (IX,IY+1) = PNSPBL(NFL,XMD,YD)

        GG2BLKD (IX,IY+1) = GG2BL(NFL,XMD,YD)

        FG2BLKD (IX,IY+1) = FG2BL(NFL,XMD,YD)

        GF2BLKD (IX,IY+1) = GF2BL(NFL,XMD,YD)

        SFF2BLKD (IX,IY+1) = SFF2BL(NFL,XMD,YD)
 

	ELSEIF (IX.LT.IV.AND.IY.GT.IX.AND.IY.LT.IV) THEN

        
           IF (ABS(YD1-XMD1).LE.TINY) THEN
              
        GG2BLCX (IX,IY) = 0.0

	PNSPBLCX (IX,IY) = 0.0

	PNSMBLCX (IX,IY) = 0.0
	
        ELSE

	GG2BLCX (IX,IY) = (GG2BL(NFL,XD1,YD1)
     >                    - GG2BL(NFL,YMD1,XMD1))/ET3

	PNSPBLCX (IX,IY) = (PNSPBL(NFL,XD1,YD1)
     >                     - PNSPBL(NFL,YMD1,XMD1))/ET2

	PNSMBLCX (IX,IY) = (PNSMBL(NFL,XD1,YD1)
     >                     - PNSMBL(NFL,YMD1,XMD1))/ET2

        ENDIF

	PNSMBLK (IX,IY) = PNSMBL(NFL,XD,YD)

	PNSPBLK (IX,IY) = PNSPBL(NFL,XD,YD)

	GG2BLK (IX,IY) = GG2BL(NFL,XD,YD)

	FG2BLK1 (IX,IY) = FG2BL(NFL,XD,YD)

	GF2BLK1 (IX,IY) = GF2BL(NFL,XD,YD)

	SFF2BLK1 (IX,IY) = SFF2BL(NFL,XD,YD)

C
C	BL region for the integration from 0 to x! 
C 

	ELSEIF (IX.GT.1.AND.IX.LE.IV.AND.IY.LT.IX.AND.IY.GT.1) THEN

        
        IF (ABS(YMD1-XD1).LE.TINY) THEN

        PNSMBLC1 (IX,IY) = 0.0

	PNSPBLC1 (IX,IY) = 0.0

	GG2BLC1 (IX,IY) = 0.0
 
           ELSE

	PNSMBLC1 (IX,IY) = (PNSMBL(NFL,XMD1,YMD1)
     >                    -PNSMBL(NFL,YD1,XD1))/ET2

	PNSPBLC1 (IX,IY) = (PNSPBL(NFL,XMD1,YMD1)
     >                    -PNSPBL(NFL,YD1,XD1))/ET2

	GG2BLC1 (IX,IY) = (GG2BL(NFL,XMD1,YMD1)
     >                    -GG2BL(NFL,YD1,XD1))/ET3

        ENDIF

	PNSMBLK (IX,IY) = PNSMBL(NFL,XMD,YMD)

	PNSPBLK (IX,IY) = PNSPBL(NFL,XMD,YMD)

	GG2BLK (IX,IY) = GG2BL(NFL,XMD,YMD)

	FG2BLK (IX,IY) = -FG2BL(NFL,XMD,YMD)

	GF2BLK (IX,IY) = -GF2BL(NFL,XMD,YMD)

	SFF2BLK (IX,IY) = SFF2BL(NFL,XMD,YMD)

C
C	DGLAP Region! x dependent kernels only! The factors of 1/y stem from the  
C       integral measure in the DGLAP region of dy/y, where the 1/y is now shifted into  
C       the kernels. An array like PNSM(IN1,IN) represents the counterterms for the  
C       regularized kernels like P^QQ(y,del) - P^QQ(y,y*del/x) - (x/del)*V^qq(y*x/del,x/del)
C       which are x dependent 
C       and have to be integrated from del..1! Note that the extra 1/y in f.ex.  
C       PNSM(IN1,IN) is due to the extra factor of x in the definition of the QQ and QG 
C       kernels. P^QQ(x/y,del/y) = x/y*p^QQ(x/y,del/y) with p^QQ from the hep/ph paper 
C       mentioned above. 
C

	ELSEIF (IX.GT.IV .AND. IY.EQ.NX) THEN

	IN1 = IY - IV
        Y = XV(IY)

	PNSM(IN,IN1) = PNSMD(NFL,XY,DY)/Y

        PNSP(IN,IN1) = PNSPD(NFL,XY,DY)/Y

	SFF2(IN,IN1) = SFF2D(NFL,XY,DY)/Y

	GG2(IN,IN1) = GG2D(NFL,XY,DY)/Y

	FG2(IN,IN1) = FG2D(NFL,XY,DY)/Y

	GF2(IN,IN1) = GF2D(NFL,XY,DY)/Y

	PNSM(IN1,IN) = 0.0

	PNSP(IN1,IN) = 0.0

	GG2(IN1,IN) = 0.0

	ELSEIF (IX .GT. IV .AND. IY .GT. IX .AND. IY.LT.NX) THEN

	IN1 = IY - IV
	Y = XV(IY)
	YP = Y*DEL/X

        Y1 = XV(IY)
	YP1 = Y1*DEL/X1

	PNSM(IN,IN1) = PNSMD(NFL,XY,DY)/Y

	PNSP(IN,IN1) = PNSPD(NFL,XY,DY)/Y

	SFF2(IN,IN1) = SFF2D(NFL,XY,DY)/Y

	GG2(IN,IN1) = GG2D(NFL,XY,DY)/Y

	FG2(IN,IN1) = FG2D(NFL,XY,DY)/Y

	GF2(IN,IN1) = GF2D(NFL,XY,DY)/Y

C     counterterms 

        PNSM(IN1,IN) = - (PNSMD(NFL,Y1,YP1)/Y1)/ET2

        PNSP(IN1,IN) = - (PNSPD(NFL,Y1,YP1)/Y1)/ET2

        GG2(IN1,IN) = - GG2D(NFL,Y1,YP1)/ET1

	ENDIF

	ENDIF

        ENDIF

C     counterterms in DGLAP region

 103    IF (Y1.GT.DX1 .AND. Y1.LT.1.0) THEN

            IF (K1.EQ.1) THEN

        PNSMC(IX,IY-1) = (PNSMCD(NFL,DX1)/DX1)/ET2

        PNSPC(IX,IY-1) = (PNSPCD(NFL,DX1)/DX1)/ET2

        GG2C(IX,IY-1) = GG2CD(NFL,DX1)/ET1

        K1 = K1 + 1

        GOTO 122

        ELSE

        GOTO 122

        ENDIF

 122    PNSPC(IX,IY) = (PNSPD(NFL,Y1,DX1)/Y1)/ET2

        PNSMC(IX,IY) = (PNSMD(NFL,Y1,DX1)/Y1)/ET2

        GG2C(IX,IY) = GG2D(NFL,Y1,DX1)/ET1

        ENDIF

 101	CONTINUE
        
C     counterterms in DGLAP region

        IF (IX.GE.IV+1) then

        DO 202 IT = 2,NX-1

           IF(XV(IT).GT.DX) THEN

              ITEMP = IT

              GOTO 203

              ENDIF

           DYZOX = XV(IT)*XD

           PNSPDGBL(IX,IT) = XD*PNSPBL(NFL,DYZOX,XD)/ET2

           PNSMDGBL(IX,IT) = XD*PNSMBL(NFL,DYZOX,XD)/ET2

           GG2DGBL(IX,IT) = XD*GG2BL(NFL,DYZOX,XD)/ET3


 202       CONTINUE
 
C     
C  Value of Counterterm V(y*x/del,x/del)*x/del at y=0
C           

 203     PNSPDGBL(IX,1) = 0D0  

         PNSMDGBL(IX,1) = 0D0

         GG2DGBL(IX,1) = 0D0

C     
C  Value of Counterterm V(y*x/del,x/del)*x/del at y=del/x
C           

         XTE = 0.99999999D0

         PNSPDGBL(IX,ITEMP) = XD*PNSPBL(NFL,XTE,XD)/ET2

         PNSMDGBL(IX,ITEMP) = XD*PNSMBL(NFL,XTE,XD)/ET2

         GG2DGBL(IX,ITEMP) = XD*GG2BL(NFL,XTE,XD)/ET3
        
         ENDIF

        IF (IX.LE.IV) THEN 

	XX = 0.9999999
        SM = 0.000000001

C
C     Counterterms for integral x..d at y=del!
C

	GG2BLCX (IX,IV) = GG2L(NFL,XD1)/ET3

	PNSPBLCX (IX,IV) = PNSPL(NFL,XD1)/ET2

	PNSMBLCX (IX,IV) = PNSML(NFL,XD1)/ET2

C
C     Terms for integral del..1 at y=del!
C

        GG2BLKD (IX,IV+1) =  GG2L(NFL,XMD)

        FG2BLKD (IX,IV+1) =  FG2L(NFL,XMD)

        GF2BLKD (IX,IV+1) =  GF2L(NFL,XMD)

        SFF2BLKD (IX,IV+1) =  SFF2L(NFL,XMD)

        PNSMBLKD (IX,IV+1) = PNSML(NFL,XMD)

        PNSPBLKD (IX,IV+1) = PNSPL(NFL,XMD)

C
C     Terms for integral x..d at y=del!
C

        GG2BLK(IX,IV) = GG2L(NFL,XD)

        PNSPBLK(IX,IV) = PNSPL(NFL,XD)

        PNSMBLK(IX,IV) = PNSML(NFL,XD)

        FG2BLK1(IX,IV) = FG2L(NFL,XD)

        GF2BLK1(IX,IV) = GF2L(NFL,XD)

        SFF2BLK1(IX,IV) = SFF2L(NFL,XD)

C
C     Terms for the integral from x..d at y=x!
C

        FG2BLK1(IX,IX) = 0.0

        GF2BLK1(IX,IX) = 0.0

        SFF2BLK1(IX,IX) = 0.0 

C
C     Terms for integral from 0..x at y=x!
C

        FG2BLK (IX,IX) = 0.0

	GF2BLK (IX,IX) = 0.0

	SFF2BLK (IX,IX) = 0.0

C
C     terms for integral 0..x at y =0!
C

        GG2BLK(IX,1) = GG2L(NFL,XMD)

        PNSPBLK(IX,1) = PNSPL(NFL,XMD)

        PNSMBLK(IX,1) = PNSML(NFL,XMD)

        FG2BLK(IX,1) = FG2L(NFL,XMD)

        GF2BLK(IX,1) = GF2L(NFL,XMD)

        SFF2BLK(IX,1) = SFF2L(NFL,XMD)

C
C     counterterms for integral 0..x at y=0!
C

        GG2BLC1 (IX,1) = GG2L(NFL,XMD1)/ET3

        PNSMBLC1 (IX,1) = PNSML(NFL,XMD1)/ET2

        PNSPBLC1 (IX,1) = PNSPL(NFL,XMD1)/ET2

        ENDIF

 100    CONTINUE
	
C
C	Polarized case now! Same as above just different kernels! 
C
	ELSEIF (IKNL.EQ.-2) THEN

           DO 500 IX = 2, NX-1

        IF (IX.EQ.IV) GOTO 500

        X = XV(IX)
	XD = X*DIM1
        DX = 1./XD
	XMD = 1. - XD
        X1 = XV(IX)
	XD1 = X1*DIM1
        DX1 = 1./XD1
	XMD1 = 1. - XD1
	IN = IX - IV

	IF (IX.LE.IV) THEN

	IQ = 1

	IQ1 = NX

	ELSE

	IQ = IV+1

        IQ1 = NX

	ENDIF 

        K1 = 1

	DO 501 IY = IQ, IQ1 

           Y = XV(IY)
           Y1 = XV(IY)

          IF (IQ.EQ.IV+1.AND.IY.LT.IX) THEN

           GOTO 503

           ELSE

	IF (IY .EQ. IX .AND. IY.GT.IV) THEN

           SM = 0.0000001
           X11 = 0.9999999

           FG2(IN,IN) = 0.0

           GF2(IN,IN) = 0.0

           SFF2(IN,IN) = 0.0

           PNSM(IN,IN) = - (PNSPD(NFL,X1,DEL)/X1)/ET2

           PNSP(IN,IN) = - (PNSMD(NFL,X1,DEL)/X1)/ET2

           GG2(IN,IN) = - GG2DA(NFL,X1,DEL)/ET1

	GOTO 503

	ELSE

	XY = X/XV(IY)
	YD = XV(IY)*DIM1
	DY = 1./YD
	YMD = 1. - YD
	XOY = XMD/YMD
	DYM1 = 1. - DY 
	XYM1 = 1. - XY

        XY1 = X1/XV(IY)
	YD1 = XV(IY)*DIM1
	DY1 = DEL/XV(IY)
	YMD1 = 1. - YD1
	XOY1 = XMD1/YMD1
	DYM11 = 1. - DY1 
	XYM11 = 1. - XY1
	
C
C	BL Region for the integral from x to 1! This integral is split in two regions  
C       x..DEL and DEL..1. In the first region there will be a regularization of the  
C       kernel, the second region is singularity free! IY + 1 in order to extrapolate  
C       down to DEL at position IY for the integral from DEL..1! 
C
	
	IF (IX.LE.IV.AND.IY.GT.IV) THEN	

          
	PNSMBLK (IX,IY) = PNSPBL(NFL,XD,YD)

	PNSPBLK (IX,IY) = PNSMBL(NFL,XD,YD)

	GG2BLK (IX,IY) = GG2BLA(NFL,XD,YD)

	FG2BLK1 (IX,IY) = FG2BLA(NFL,XD,YD)

	GF2BLK1 (IX,IY) = GF2BLA(NFL,XD,YD)

	SFF2BLK1 (IX,IY) = SFF2BLA(NFL,XD,YD)

        PNSMBLKD (IX,IY+1) = PNSPBL(NFL,XMD,YD)

        PNSPBLKD (IX,IY+1) = PNSMBL(NFL,XMD,YD)

        GG2BLKD (IX,IY+1) = GG2BLA(NFL,XMD,YD)

        FG2BLKD (IX,IY+1) = FG2BLA(NFL,XMD,YD)

        GF2BLKD (IX,IY+1) = GF2BLA(NFL,XMD,YD)

        SFF2BLKD (IX,IY+1) = SFF2BLA(NFL,XMD,YD)

	ELSEIF (IX.LT.IV.AND.IY.GT.IX.AND.IY.LT.IV) THEN

        
           IF (ABS(YD1-XMD1).LE.TINY) THEN
              
        GG2BLCX (IX,IY) = 0.0

	PNSPBLCX (IX,IY) = 0.0

	PNSMBLCX (IX,IY) = 0.0
	
        ELSE

	GG2BLCX (IX,IY) = (GG2BLA(NFL,XD1,YD1)
     >                    - GG2BLA(NFL,YMD1,XMD1))/ET3

	PNSPBLCX (IX,IY) = (PNSMBL(NFL,XD1,YD1)
     >                     - PNSMBL(NFL,YMD1,XMD1))/ET2

	PNSMBLCX (IX,IY) = (PNSPBL(NFL,XD1,YD1)
     >                     - PNSPBL(NFL,YMD1,XMD1))/ET2

        ENDIF

	PNSMBLK (IX,IY) = PNSPBL(NFL,XD,YD)

	PNSPBLK (IX,IY) = PNSMBL(NFL,XD,YD)

	GG2BLK (IX,IY) = GG2BLA(NFL,XD,YD)

	FG2BLK1 (IX,IY) = FG2BLA(NFL,XD,YD)

	GF2BLK1 (IX,IY) = GF2BLA(NFL,XD,YD)

	SFF2BLK1 (IX,IY) = SFF2BLA(NFL,XD,YD)

C
C	BL region for the integration from 0 to x! 
C 

	ELSEIF (IX.GT.1.AND.IX.LE.IV.AND.IY.LT.IX.AND.IY.GT.1) THEN


        IF (ABS(YMD1-XD1).LE.TINY) THEN

        PNSMBLC1 (IX,IY) = 0.0

	PNSPBLC1 (IX,IY) = 0.0

	GG2BLC1 (IX,IY) = 0.0
 
           ELSE

	PNSMBLC1 (IX,IY) = (PNSPBL(NFL,XMD1,YMD1)
     >                    -PNSPBL(NFL,YD1,XD1))/ET2

	PNSPBLC1 (IX,IY) = (PNSMBL(NFL,XMD1,YMD1)
     >                    -PNSMBL(NFL,YD1,XD1))/ET2

	GG2BLC1 (IX,IY) = (GG2BLA(NFL,XMD1,YMD1)
     >                    -GG2BLA(NFL,YD1,XD1))/ET3

        ENDIF

	PNSMBLK (IX,IY) = PNSPBL(NFL,XMD,YMD)

	PNSPBLK (IX,IY) = PNSMBL(NFL,XMD,YMD)

	GG2BLK (IX,IY) = GG2BLA(NFL,XMD,YMD)

	FG2BLK (IX,IY) = -FG2BLA(NFL,XMD,YMD)

	GF2BLK (IX,IY) = -GF2BLA(NFL,XMD,YMD)

	SFF2BLK (IX,IY) = SFF2BLA(NFL,XMD,YMD)

C
C	DGLAP Region! x dependent kernels only! The factors of 1/y stem from the  
C       integral measure in the DGLAP region of dy/y, where the 1/y is now shifted into  
C       the kernels. An array like PNSM(IN1,IN) represents the counterterms for the  
C       regularized kernels like P^QQ(y,del) - P^QQ(y,y*del/x) - (x/del)*V^qq(y*x/del,x/del)
C       which are x dependent 
C       and have to be integrated from del..1! Note that the extra 1/y in f.ex.  
C       PNSM(IN1,IN) is due to the extra factor of x in the definition of the QQ and QG 
C       kernels. P^QQ(x/y,del/y) = x/y*p^QQ(x/y,del/y) with p^QQ from the hep/ph paper 
C       mentioned above. 
C

	ELSEIF (IX.GT.IV .AND. IY.EQ.NX) THEN

	IN1 = IY - IV
        Y = XV(IY)

	PNSM(IN,IN1) = PNSPD(NFL,XY,DY)/Y

        PNSP(IN,IN1) = PNSMD(NFL,XY,DY)/Y

	SFF2(IN,IN1) = SFF2DA(NFL,XY,DY)/Y

	GG2(IN,IN1) = GG2DA(NFL,XY,DY)/Y

	FG2(IN,IN1) = FG2DA(NFL,XY,DY)/Y

	GF2(IN,IN1) = GF2DA(NFL,XY,DY)/Y

	PNSM(IN1,IN) = 0.0

	PNSP(IN1,IN) = 0.0

	GG2(IN1,IN) = 0.0

	ELSEIF (IX .GT. IV .AND. IY .GT. IX .AND. IY.LT.NX) THEN

	IN1 = IY - IV
	Y = XV(IY)
	YP = Y*DEL/X

        Y1 = XV(IY)
	YP1 = Y1*DEL/X1

	PNSM(IN,IN1) = PNSPD(NFL,XY,DY)/Y

	PNSP(IN,IN1) = PNSMD(NFL,XY,DY)/Y

	SFF2(IN,IN1) = SFF2DA(NFL,XY,DY)/Y

	GG2(IN,IN1) = GG2DA(NFL,XY,DY)/Y

	FG2(IN,IN1) = FG2DA(NFL,XY,DY)/Y

	GF2(IN,IN1) = GF2DA(NFL,XY,DY)/Y

        PNSM(IN1,IN) = - (PNSPD(NFL,Y1,YP1)/Y1)/ET2

        PNSP(IN1,IN) = - (PNSMD(NFL,Y1,YP1)/Y1)/ET2

        GG2(IN1,IN) = - GG2DA(NFL,Y1,YP1)/ET1

	ENDIF

	ENDIF

        ENDIF

 503    IF (Y1.GT.DX1 .AND. Y1.LT.1.0) THEN

            IF (K1.EQ.1) THEN

        PNSMC(IX,IY-1) = (PNSPCD(NFL,DX1)/DX1)/ET2

        PNSPC(IX,IY-1) = (PNSMCD(NFL,DX1)/DX1)/ET2

        GG2C(IX,IY-1) = GG2CDA(NFL,DX1)/ET1

        K1 = K1 + 1

        GOTO 522

        ELSE

        GOTO 522

        ENDIF

 522    PNSPC(IX,IY) = (PNSMD(NFL,Y1,DX1)/Y1)/ET2

        PNSMC(IX,IY) = (PNSPD(NFL,Y1,DX1)/Y1)/ET2

        GG2C(IX,IY) = GG2DA(NFL,Y1,DX1)/ET1

        ENDIF

 501	CONTINUE

        IF (IX.GE.IV+1) then

        DO 602 IT = 2,NX-1

           IF(XV(IT).GT.DX) THEN

              ITEMP = IT

              GOTO 603

              ENDIF

           DYZOX = XV(IT)*XD

           PNSPDGBL(IX,IT) = XD*PNSMBL(NFL,DYZOX,XD)/ET2

           PNSMDGBL(IX,IT) = XD*PNSPBL(NFL,DYZOX,XD)/ET2

           GG2DGBL(IX,IT) = XD*GG2BLA(NFL,DYZOX,XD)/ET3

 602       CONTINUE
 
C     
C  Value of Counterterm V(y*x/del,x/del)*x/del at y=0
C           

 603     PNSPDGBL(IX,0) = 0D0 

         PNSMDGBL(IX,0) = 0D0

         GG2DGBL(IX,0) = 0D0

C     
C  Value of Counterterm V(y*x/del,x/del)*x/del at y=del/x
C           

         XTE = 0.99999999D0

         PNSPDGBL(IX,ITEMP) = XD*PNSMBL(NFL,XTE,XD)/ET2

         PNSMDGBL(IX,ITEMP) = XD*PNSPBL(NFL,XTE,XD)/ET2

         GG2DGBL(IX,ITEMP) =  XD*GG2BLA(NFL,XTE,XD)/ET3

         ENDIF

        
	IF (IX.LE.IV) THEN 

	XX = 0.9999999
        SM = 0.000000001

C
C     Counterterms for integral x..d at y=del!
C

	GG2BLCX (IX,IV) = GG2LA(NFL,XD1)/ET3

	PNSPBLCX (IX,IV) = PNSML(NFL,XD1)/ET2

	PNSMBLCX (IX,IV) = PNSPL(NFL,XD1)/ET2

C
C     Terms for integral del..1 at y=del!
C

        GG2BLKD (IX,IV+1) =  GG2LA(NFL,XMD)

        FG2BLKD (IX,IV+1) =  FG2LA(NFL,XMD)

        GF2BLKD (IX,IV+1) =  GF2LA(NFL,XMD)

        SFF2BLKD (IX,IV+1) =  SFF2LA(NFL,XMD)

        PNSMBLKD (IX,IV+1) = PNSPL(NFL,XMD)

        PNSPBLKD (IX,IV+1) = PNSML(NFL,XMD)

C
C     Terms for integral x..d at y=del!
C

        GG2BLK(IX,IV) = GG2LA(NFL,XD)

        PNSPBLK(IX,IV) = PNSML(NFL,XD)

        PNSMBLK(IX,IV) = PNSPL(NFL,XD)

        FG2BLK1(IX,IV) = FG2LA(NFL,XD)

        GF2BLK1(IX,IV) = GF2LA(NFL,XD)

        SFF2BLK1(IX,IV) = SFF2LA(NFL,XD)

C
C     Terms for the integral from x..d at y=x!
C

        FG2BLK1(IX,IX) = 0.0

        GF2BLK1(IX,IX) = 0.0

        SFF2BLK1(IX,IX) = 0.0 

C
C     Terms for integral from 0..x at y=x!
C

        FG2BLK (IX,IX) = 0.0

	GF2BLK (IX,IX) = 0.0

	SFF2BLK (IX,IX) = 0.0

C
C     terms for integral 0..x at y =0!
C

        GG2BLK(IX,1) = GG2LA(NFL,XMD)

        PNSPBLK(IX,1) = PNSML(NFL,XMD)

        PNSMBLK(IX,1) = PNSPL(NFL,XMD)

        FG2BLK(IX,1) = FG2LA(NFL,XMD)

        GF2BLK(IX,1) = GF2LA(NFL,XMD)

        SFF2BLK(IX,1) = SFF2LA(NFL,XMD)

C
C     counterterms for integral 0..x at y=0!
C

        GG2BLC1 (IX,1) = GG2LA(NFL,XMD1)/ET3

        PNSMBLC1 (IX,1) = PNSPL(NFL,XMD1)/ET2

        PNSPBLC1 (IX,1) = PNSML(NFL,XMD1)/ET2

        ENDIF
              
 500	CONTINUE

	ENDIF

        KT = 3

        DO 222 IX = 2,IV-1

           IF (IX.GT.3) THEN

        DO 900 I=KT,1,-1

           T(I) =  GG2BLC1 (IX,IX-I)
           T1(I) =  PNSMBLC1 (IX,IX-I)
           T2(I) =  PNSPBLC1 (IX,IX-I)
           XX2(I) = XV(IX-I)

 900       CONTINUE

           XX = XV(IX)

           CALL POLINT(XX2,T,KT,XX,TEM,ERR)
           GG2BLC1 (IX,IX) = TEM

           CALL POLINT(XX2,T1,KT,XX,TEM,ERR)
           PNSMBLC1 (IX,IX) = TEM

           CALL POLINT(XX2,T2,KT,XX,TEM,ERR)
           PNSPBLC1 (IX,IX) = TEM

           ENDIF

           IF (IX.LT.IV-2) THEN

           DO 901 I=KT,1,-1

           T(I) =  GG2BLCX (IX,IX+I)
           T1(I) =  PNSMBLCX (IX,IX+I)
           T2(I) =  PNSPBLCX (IX,IX+I)
           XX2(I) = XV(IX+I)
         
 901       CONTINUE

           XX = XV(IX)

           CALL POLINT(XX2,T,KT,XX,TEM,ERR)
           GG2BLCX (IX,IX) = TEM

           CALL POLINT(XX2,T1,KT,XX,TEM,ERR)
           PNSMBLCX (IX,IX) = TEM

           CALL POLINT(XX2,T2,KT,XX,TEM,ERR)
           PNSPBLCX (IX,IX) = TEM

           ENDIF

 222       CONTINUE
           


        RETURN

C                        ****************************
      END
C 	****************************************** 
C
C	Function to allow the use of Gauss integration routine for the counterterms 
C	by using interpolation routine splint (polynomial fit in bin) or fintrp 
C       (more elaborate but slower)! 
C

	FUNCTION PFF1 (IX,IX1,XX) 
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
       
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1) 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2) 
      PARAMETER (MX = 3) 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX, DEL, IV 
      COMMON / KRNL00 / DZ, XL(MX), NNX 
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX) 
      COMMON / KRN1ST / FF1(0:MXX,0:MXX), GG1(0:MXX,0:MXX),    
     > FG2(0:MXX,0:MXX),GF2(0:MXX,0:MXX),GG2(0:MXX,0:MXX), 
     > PNSP(0:MXX,0:MXX), PNSM(0:MXX,0:MXX), SFF2(0:MXX,0:MXX), 
     > SFF2BLK(0:MXX,0:MXX), FG2BLK(0:MXX,0:MXX), GF2BLK(0:MXX,0:MXX),  
     > GG2BLK(0:MXX,0:MXX), PNSPBLK(0:MXX,0:MXX), PNSMBLK(0:MXX,0:MXX), 
     > GG2BLC1(0:MXX,0:MXX), PNSPBLC1(0:MXX,0:MXX),FF1BLX(0:MXX,0:MXX),  
     > GG2BLCX(0:MXX,0:MXX), PNSPBLCX(0:MXX,0:MXX),GG1BLX(0:MXX,0:MXX), 
     > PNSMBLC1(0:MXX,0:MXX), PNSMBLCX(0:MXX,0:MXX),FF1BL1(0:MXX,0:MXX), 
     > GG1BL1(0:MXX,0:MXX),SFF2BLKD(0:MXX,0:MXX), FG2BLKD(0:MXX,0:MXX), 
     > GF2BLKD(0:MXX,0:MXX), GG2BLKD(0:MXX,0:MXX),PNSPBLKD(0:MXX,0:MXX), 
     > PNSMBLKD(0:MXX,0:MXX),FG2BLK1(0:MXX,0:MXX),GF2BLK1(0:MXX,0:MXX),
     > SFF2BLK1(0:MXX,0:MXX),PNSPC(0:MXX,0:MXX),PNSMC(0:MXX,0:MXX),
     > GG2C(0:MXX,0:MXX), PNSPDGBL(0:MXX,0:MXX), PNSMDGBL(0:MXX,0:MXX),
     > GG2DGBL(0:MXX,0:MXX)
      COMMON / DERIVAR / DPP(3,MXX),DPM(3,MXX),DGG(3,MXX),DPP1(3,MXX)
     > ,DPM1(3,MXX),
     > DGG1(3,MXX),DPPB(3,MXX),DPMB(3,MXX),DPPB1(3,MXX),DPMB1(3,MXX)
     > ,DGGB1(3,MXX),DGGB(3,MXX),DDGBLGG(3,MXX),DDGBLPM(3,MXX)
     > ,DDGBLPP(3,MXX)



      RETURN 
C                        ---------------------------- 
      ENTRY FNSP (IX,IX1,XX) 
 
      X = XX 
 
      IF (X .GT. D1) THEN 
 
        FNSP = 0 
        RETURN 
 
      ElseIF (X .GE. DEL) THEN  
        
         CALL SPLINT(DPP1,IX1,X,TEM)
           
      EndIf 

           A = PNSP(IX1-IV,IX)
            B = PNSP(IX1-IV+1,IX)
            IT = IV

            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN

           IT1 = IX1-IT
           
         TEM = FINTRP (DEL,IV,IT1,IX,PNSP,XV,-DZ,2,NX,X,ERR,IRT) 

	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

           IT1 = IX1 - IT

         TEM = FINTRP (DEL,IV,IT1,IX,PNSP,XV,-DZ,2,NX,X,ERR,IRT)  

	ENDIF

        ENDIF


      FNSP = TEM 

      RETURN 
C                        ---------------------------- 
      ENTRY FNSM (IX,IX1,XX) 
 
      X = XX 
 
      IF (X .GT. D1) THEN 
 
        FNSM = 0 
        RETURN 
 
      ElseIF (X .GE. DEL) THEN  

         CALL SPLINT(DPM1,IX1,X,TEM)

      ENDIF

           A = PNSM(IX1-IV,IX)
            B = PNSM(IX1-IV+1,IX)
            IT = IV

            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN
           
         IT1 = IX1-IT 
         TEM = FINTRP (DEL,IV,IT1,IX,PNSM,XV,-DZ,2,NX,X,ERR,IRT)

	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

         IT1 = IX1 - IT
         TEM = FINTRP (DEL,IV,IT1,IX,PNSM,XV,-DZ,2,NX,X,ERR,IRT)  

	ENDIF

        ENDIF
    
      FNSM = TEM       

      RETURN 

C                        ---------------------------- 
  
      ENTRY RGG2 (IX,IX1,XX) 
 
      X = XX 
 
      IF (X .GT. D1) THEN 
 
        RGG2 = 0 
        RETURN 
 
      ElseIF (X .GE. DEL) THEN 

         CALL SPLINT(DGG1,IX1,X,TEM)
    
      EndIf 
   
            A = GG2(IX1-IV,IX)
            B = GG2(IX1-IV+1,IX)
            IT = IV

            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN
         
         IT1 = IX1 - IT
         TEM = FINTRP (DEL,IV,IT1,IX,GG2,XV,-DZ,2,NX,X,ERR,IRT)

	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

           IT1 = IX1-IT
         TEM = FINTRP (DEL,IV,IT1,IX,GG2,XV,-DZ,2,NX,X,ERR,IRT)  

	ENDIF

        ENDIF

         RGG2 = TEM 
         
      RETURN 
 
C                        ---------------------------- 
      ENTRY FNSPC (IX,IX1,XX) 
 
      X = XX 

      DX = DEL/XV(IX)
 
      IF (X .GT. D1) THEN 
 
        FNSPC = 0 
        RETURN 
 
      ElseIF (X .GE. DEL) THEN 
 
         CALL SPLINT(DPP,IX1,X,TEM)         
           
      EndIf 

            A = PNSPC(IX,IX1)
            B = PNSPC(IX,IX1+1)
     
            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN
           
         TEM = FINTRP (DX,IV,IX1,IX,PNSPC,XV,-DZ,1,NX,X,ERR,IRT) 

	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

         TEM = FINTRP (DX,IV,IX1,IX,PNSPC,XV,-DZ,1,NX,X,ERR,IRT)  

	ENDIF

        ENDIF


      FNSPC = TEM 

      RETURN 
C                        ---------------------------- 
      ENTRY FNSMC (IX,IX1,XX) 
 
      X = XX 

      DX = DEL/XV(IX)
 
      IF (X .GT. D1) THEN 
 
        FNSMC = 0 
        RETURN 
 
      ElseIF (X .GE. DEL) THEN 
  
         CALL SPLINT(DPM,IX1,X,TEM)

      ENDIF

            A = PNSMC(IX,IX1)
            B = PNSMC(IX,IX1+1)
     
            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN
           
         TEM = FINTRP (DX,IV,IX1,IX,PNSMC,XV,-DZ,1,NX,X,ERR,IRT) 

	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

         TEM = FINTRP (DX,IV,IX1,IX,PNSMC,XV,-DZ,1,NX,X,ERR,IRT)  

	ENDIF

        ENDIF
          
      FNSMC = TEM       

      RETURN 

C                        ---------------------------- 
  
      ENTRY RGG2C (IX,IX1,XX) 
 
      X = XX 

      R = 1D0

      DX = DEL/XV(IX)
 
      IF (X .GT. D1) THEN 
 
        RGG2C = 0 
        RETURN 
 
      ElseIF (X .GE. DEL) THEN  
 
         CALL SPLINT(DGG,IX1,X,TEM)
    
      EndIf 

            A = GG2C(IX,IX1)
            B = GG2C(IX,IX1+1)
     
            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN
           
         TEM = FINTRP (DX,IV,IX1,IX,GG2C,XV,R,1,NX,X,ERR,IRT) 

	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

         TEM = FINTRP (DX,IV,IX1,IX,GG2C,XV,R,1,NX,X,ERR,IRT)  
           
	ENDIF

        ENDIF
        
         RGG2C = TEM 

      RETURN 
C                        ---------------------------- 

      ENTRY FNSMDGBL (IX,IX1,XX) 
 
      X = XX 

      R = 1D0

      DX = DEL/XV(IX)
 
      IF (X .GT. D1) THEN 
 
        FNSMDGBL = 0 
        RETURN 
 
      ElseIF (X .GE. DEL) THEN  
 
         CALL SPLINT(DDGBLPM,IX1,X,TEM)
    
      EndIf 

            A = PNSMDGBL(IX,IX1)
            B = PNSMDGBL(IX,IX1+1)
     
            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN
           
         TEM = FINTRP01 (DX,IX1,IX,PNSMDGBL,XV,R,NX,X,ERR,IRT) 

	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

         TEM = FINTRP01 (DX,IX1,IX,PNSMDGBL,XV,R,NX,X,ERR,IRT)  
           
	ENDIF

        ENDIF
        
         FNSMDGBL = TEM 

      RETURN 

C       **********************
 
      ENTRY FNSPDGBL (IX,IX1,XX) 
 
      X = XX 

      R = 1D0

      DX = DEL/XV(IX)
 
      IF (X .GT. D1) THEN 
 
        FNSPDGBL = 0 
        RETURN 
 
      ElseIF (X .GE. DEL) THEN  
 
         CALL SPLINT(DDGBLPP,IX1,X,TEM)
    
      EndIf 

            A = PNSPDGBL(IX,IX1)
            B = PNSPDGBL(IX,IX1+1)
     
            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN
           
         TEM = FINTRP01 (DX,IX1,IX,PNSPDGBL,XV,R,NX,X,ERR,IRT) 

	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

         TEM = FINTRP01 (DX,IX1,IX,PNSPDGBL,XV,R,NX,X,ERR,IRT)  
           
	ENDIF

        ENDIF

         FNSPDGBL = TEM 

      RETURN 

C       ********************** 
  
      ENTRY RGG2DGBL (IX,IX1,XX) 
 
      X = XX 

      R = 1D0

      DX = DEL/XV(IX)
 
      IF (X .GT. D1) THEN 
 
        RGG2DGBL = 0 
        RETURN 
 
      ElseIF (X .GE. DEL) THEN  
 
         CALL SPLINT(DDGBLGG,IX1,X,TEM)
    
      EndIf 

            A = GG2DGBL(IX,IX1)
            B = GG2DGBL(IX,IX1+1)
     
            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN
           
         TEM = FINTRP01 (DX,IX1,IX,GG2DGBL,XV,R,NX,X,ERR,IRT) 

	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

         TEM = FINTRP01 (DX,IX1,IX,GG2DGBL,XV,R,NX,X,ERR,IRT)  
           
	ENDIF

        ENDIF
        
         RGG2DGBL = TEM 

      RETURN 

C       ********************** 

      ENTRY FNSPB (IX,IX1,XX) 
 
      X = XX 
      ID = 1
 
      IF (X .LT. 0.0 .OR. X.GT.DEL) THEN 
 
        FNSPB = 0 
        RETURN 
 
      ElseIF (X .GE. 0.0) THEN  
         
         CALL SPLINT(DPPB,IX1,X,TEM)
     
      EndIf 

            A = PNSPBLC1(IX,IX1)
            B = PNSPBLC1(IX,IX1+1)

            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN
           
        TEM = FINTRP2 (ID,DEL,IV,IX1,IX,PNSPBLC1,XV,-DZ,DZ,NX,X,ERR,IRT)

	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

        TEM = FINTRP2 (ID,DEL,IV,IX1,IX,PNSPBLC1,XV,-DZ,DZ,NX,X,ERR,IRT)

	ENDIF

        ENDIF

         FNSPB = TEM 

      RETURN 
C                        ---------------------------- 
      ENTRY FNSMB (IX,IX1,XX) 
 
      X = XX 
      ID = 1 
 
      IF (X .LT. 0.0 .OR. X.GT.DEL) THEN 
 
        FNSMB = 0 
        RETURN 
 
      ElseIF (X .GE. 0.0) THEN 
      
         CALL SPLINT(DPMB,IX1,X,TEM)

         endif

            A = PNSMBLC1(IX,IX1)
            B = PNSMBLC1(IX,IX1+1)

            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN
           
        TEM = FINTRP2 (ID,DEL,IV,IX1,IX,PNSMBLC1,XV,-DZ,DZ,NX,X,ERR,IRT)

	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

        TEM = FINTRP2 (ID,DEL,IV,IX1,IX,PNSMBLC1,XV,-DZ,DZ,NX,X,ERR,IRT)

	ENDIF

        ENDIF

      FNSMB = TEM 

      RETURN 
C                        ---------------------------- 
 
      ENTRY RGG2B (IX,IX1,XX) 
 
      X = XX 
      ID = 1
 
      IF (X .LT. 0.0 .OR. X.GT.DEL) THEN 
 
        RGG2B = 0 
        RETURN 
 
      ElseIF (X .GE. 0.0) THEN 
 
         CALL SPLINT(DGGB,IX1,X,TEM)
    
      EndIf 

            A = GG2BLC1(IX,IX1)
            B = GG2BLC1(IX,IX1+1)

            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN
           
        TEM = FINTRP2 (ID,DEL,IV,IX1,IX,GG2BLC1,XV,-DZ,DZ,NX,X,ERR,IRT)

	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

        TEM = FINTRP2 (ID,DEL,IV,IX1,IX,GG2BLC1,XV,-DZ,DZ,NX,X,ERR,IRT)

	ENDIF

        ENDIF

        RGG2B = TEM

      RETURN 
 
C       ********************** 

      ENTRY FNSPB1 (IX,IX1,XX) 
 
      X = XX
      ID = 2
 
      IF (X .LT. 0.0 .OR. X.GT.DEL) THEN 
 
        FNSPB1 = 0 
        RETURN 
 
      ElseIF (X .GE. 0.0) THEN  
      
        CALL SPLINT(DPPB1,IX1,X,TEM)

      EndIf 

            A = PNSPBLCX(IX,IX1)
            B = PNSPBLCX(IX,IX1+1)

            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN
           
        TEM = FINTRP2 (ID,DEL,IV,IX1,IX,PNSPBLCX,XV,-DZ,DZ,NX,X,ERR,IRT)

	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

        TEM = FINTRP2 (ID,DEL,IV,IX1,IX,PNSPBLCX,XV,-DZ,DZ,NX,X,ERR,IRT)

	ENDIF

        ENDIF

        FNSPB1 = TEM 

      RETURN 
C                        ---------------------------- 
      ENTRY FNSMB1 (IX,IX1,XX) 
 
      X = XX
      ID = 2
 
      IF (X .LT. 0.0 .OR. X.GT.DEL) THEN 
 
        FNSMB1 = 0 
        RETURN 
 
      ElseIF (X .GE. 0.0) THEN  
      
         CALL SPLINT(DPMB1,IX1,X,TEM)

      EndIf 

            A = PNSMBLCX(IX,IX1)
            B = PNSMBLCX(IX,IX1+1)

            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN
           
        TEM = FINTRP2 (ID,DEL,IV,IX1,IX,PNSMBLCX,XV,-DZ,DZ,NX,X,ERR,IRT)
          
	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

        TEM = FINTRP2 (ID,DEL,IV,IX1,IX,PNSMBLCX,XV,-DZ,DZ,NX,X,ERR,IRT)
        
	ENDIF

        ENDIF

         FNSMB1 = TEM 

      RETURN 

C                        ---------------------------- 
   
      ENTRY RGG2B1 (IX,IX1,XX) 
 
      X = XX
      ID = 2
 
      IF (X .LT. 0.0 .OR. X.GT.DEL) THEN 
 
        RGG2B1 = 0 
        RETURN 
 
      ElseIF (X .GE. 0.0) THEN 
             
         CALL SPLINT(DGGB1,IX1,X,TEM)

      EndIf 

            A = GG2BLCX(IX,IX1)
            B = GG2BLCX(IX,IX1+1)

            IF (A.LT.B) THEN

	IF (TEM.LT.A.OR.TEM.GT.B) THEN

        TEM = FINTRP2 (ID,DEL,IV,IX1,IX,GG2BLCX,XV,-DZ,DZ,NX,X,ERR,IRT)
       
	ENDIF

	ELSEIF (A.GT.B) THEN

        IF (TEM.LT.B.OR.TEM.GT.A) THEN

        TEM = FINTRP2 (ID,DEL,IV,IX1,IX,GG2BLCX,XV,-DZ,DZ,NX,X,ERR,IRT)

	ENDIF

        ENDIF

         RGG2B1 = TEM 

      RETURN 

C       ********************** 
      END 
C	********************** 

	SUBROUTINE STUPKL (NFL) 
 
C                                 Set up the common block containing the arrays 
C                                        for the first and second order kernels 
C 
C                                 Also calculates the integrals and constants 
C                                   needed to complete the [p(x,del)]sub+ integrals 
 
C                         Real calculation is done in the routine KERNEL which 
C                           follows strictly the BFM notation. 
C 
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      LOGICAL LSTX 
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1) 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (MX = 3) 
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2) 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX, DEL,IV 
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX) 
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX) 
      COMMON / KRN1ST / FF1(0:MXX,0:MXX), GG1(0:MXX,0:MXX),    
     >FG2(0:MXX,0:MXX),GF2(0:MXX,0:MXX),GG2(0:MXX,0:MXX), 
     >PNSP(0:MXX,0:MXX), PNSM(0:MXX,0:MXX), SFF2(0:MXX,0:MXX), 
     >SFF2BLK(0:MXX,0:MXX), FG2BLK(0:MXX,0:MXX), GF2BLK(0:MXX,0:MXX),  
     >GG2BLK(0:MXX,0:MXX), PNSPBLK(0:MXX,0:MXX), PNSMBLK(0:MXX,0:MXX), 
     >GG2BLC1(0:MXX,0:MXX), PNSPBLC1(0:MXX,0:MXX),FF1BLX(0:MXX,0:MXX),  
     >GG2BLCX(0:MXX,0:MXX), PNSPBLCX(0:MXX,0:MXX),GG1BLX(0:MXX,0:MXX), 
     >PNSMBLC1(0:MXX,0:MXX), PNSMBLCX(0:MXX,0:MXX),FF1BL1(0:MXX,0:MXX), 
     >GG1BL1(0:MXX,0:MXX),SFF2BLKD(0:MXX,0:MXX), FG2BLKD(0:MXX,0:MXX), 
     >GF2BLKD(0:MXX,0:MXX), GG2BLKD(0:MXX,0:MXX),PNSPBLKD(0:MXX,0:MXX), 
     >PNSMBLKD(0:MXX,0:MXX),FG2BLK1(0:MXX,0:MXX),GF2BLK1(0:MXX,0:MXX),
     >SFF2BLK1(0:MXX,0:MXX),PNSPC(0:MXX,0:MXX),PNSMC(0:MXX,0:MXX),
     >GG2C(0:MXX,0:MXX), PNSPDGBL(0:MXX,0:MXX), PNSMDGBL(0:MXX,0:MXX),
     > GG2DGBL(0:MXX,0:MXX)
      COMMON / KRNL00 / DZ, XL(MX), NNX 
      COMMON / KRNL01 / AGG2(0:MXX), ANSP(0:MXX), ANSM(0:MXX) 
      COMMON / DERIVAR / DPP(3,MXX),DPM(3,MXX),DGG(3,MXX),
     > DPP1(3,MXX),DPM1(3,MXX),
     > DGG1(3,MXX),DPPB(3,MXX),DPMB(3,MXX),DPPB1(3,MXX),DPMB1(3,MXX)
     > ,DGGB1(3,MXX),DGGB(3,MXX),DDGBLGG(3,MXX),DDGBLPM(3,MXX)
     > ,DDGBLPP(3,MXX)

      EXTERNAL PFF1, RGG2, FNSP, FNSM,
     > RGG2B, RGG2B1, FNSPB, FNSPB1, FNSMB, FNSMB1,RGG2C, FNSPC, FNSMC,
     > RGG2DGBL,FNSMDGBL, FNSPDGBL
 
      DIMENSION X2(10),FX(10), TFM(MXX), TFP(MXX), TGG(MXX),T1FM(MXX), 
     >T1FP(MXX), T1GG(MXX), TFMD(MXX), TFPD(MXX), TGGD(MXX),T1FMD(MXX),
     >T1FPD(MXX), T1GGD(MXX),TDGBLGG(MXX),TDGBLPM(MXX),TDGBLPP(MXX),
     >FX1(10),FX2(10) 	 
 
      SAVE 
      DATA AERR, RERR / 0.0, 0.02 / 
 
 
C 
C                          First past these quantities to FUNCTION PFF1 ... etc 
C                   via / KRNL00 / for interpolation and extrapolation purposes 
      NNX = NX 
      DIM = 1D0/DEL 

	DZ = 1./(NX-1.) 
        IF (DEL.GE.0.001 .AND. DEL.LT.0.03) THEN
	ET2 = 1000.
	ET1 = 100.
        ET3 = 1000.
        ELSEIF (DEL.GE.0.03 .AND. DEL.LT.0.1) THEN
        ET2 = 100.
	ET1 = 50.
        ET3 = 100.  
        ELSEIF (DEL.GE.0.1) THEN
        ET2 = 100.
	ET1 = 50.
        ET3 = 100.  
        ELSEIF (DEL.LT.0.001) THEN
           ET2 = 1000.
           ET1 = 100.
           ET3 = 1000.
        ENDIF
 
C 
C      Calculation of Kernels in subroutine KERNEL! 
C  
       CALL KERNEL (FF1, GG1, PNSPBLK, PNSPBLC1, PNSPBLCX, PNSMBLK, 
     > PNSMBLC1, PNSMBLCX, SFF2BLK, GG2BLK, GG2BLC1, GG2BLCX, FG2BLK, 
     > GF2BLK, PNSP, PNSM, SFF2, FG2,GF2, GG2, FF1BLX,GG1BLX,FF1BL1,
     > GG1BL1,PNSMBLKD,PNSPBLKD,SFF2BLKD,GG2BLKD,FG2BLKD, GF2BLKD,
     > FG2BLK1,GF2BLK1,SFF2BLK1,PNSPC,PNSMC,GG2C,PNSPDGBL,PNSMDGBL,
     > GG2DGBL,NFL, IRT) 
 
C 
C     Int of PNSP, and PNSM etc. to determine counterterms! 
C                                        

C 
C     Set counterterms for NX equal 0! And set integration precision
C     for GAUSS integration
C 
      IF (DEL.GE.0.001 .AND. DEL.LT.0.03) THEN
	AER = 0.000001
        RER = 0.00002
        AER1 = 0.0000001
        RER1 = 0.000002
        AER2 = 0.000001
        RER2 = 0.00001
        AER3 = 0.000001
        RER3 = 0.00001
        ELSEIF(DEL.Ge.0.03 .AND. DEL.LT.0.1) THEN
          AER = 0.000001
          RER = 0.00002 
          AER1 = 0.0000001
          RER1 = 0.000002 
          AER2 = 0.000001
          RER2 = 0.00001
          AER3 = 0.000001
          RER3 = 0.00001
         ELSEIF(DEL.Ge.0.1) THEN
          AER = 0.000001
          RER = 0.00002 
          AER1 = 0.0000001
          RER1 = 0.000002 
          AER2 = 0.000001
          RER2 = 0.00001
          AER3 = 0.000001
          RER3 = 0.00001
        ELSE
           AER = 0.000001
           AER1 = 0.0000001
           RER = 0.00002
           RER1 = 0.000002
           AER2 = 0.0000001
           RER2 = 0.000001
           AER3 = 0.000001
           RER3 = 0.00001
        ENDIF

	AGG2(NX) = 0D0
	
	ANSM(NX) = 0D0
        
	ANSP(NX) = 0D0
	
C
C	Initialize a couple of variables used in the integration to store intermediate 
C       results!  
C 
 	
	TEMFM2 = 0 
	TEMFM22 = 0 
	TEMFP2 = 0 
	TEMFP22 = 0 
	TEMG2 = 0 
	TEMG22 = 0 
 
	TEMF = 0 
	TEMG = 0 
	TEMFP = 0 
	TEMGP = 0 

C
C     Start integration loops!
C

      
      DO 21 I1 = IV+1,NX-1 

         IN = I1 - IV 

         DX = DEL/XV(I1)
C
C     Does the integration of V(y*x/del,x/del)*x/del from 0..del/x for fixed x!
C


         DO 833 I3 = IV,NX-1 

            IF (XV(I3).GE.DX) THEN

               KK1 = I3
               GOTO 899

            ENDIF

 833     CONTINUE   

C
C     Set up array for interpolation!
C


 899    CALL SPLINED2(XV,PNSPDGBL,KK1,1,I1,0,DX,DDGBLPP)
        CALL SPLINED2(XV,PNSMDGBL,KK1,1,I1,0,DX,DDGBLPM)
        CALL SPLINED2(XV,GG2DGBL,KK1,1,I1,0,DX,DDGBLGG)

         DO 911 I2 = 1,KK1-1

            IF (I2.EQ.KK1-1) THEN 

C
C     Integration from d/x..x_i2
C

      TEMFP2 = GAUSS(I1,I2,FNSPDGBL,XV(I2),DX,AER,RER,ER21,IRT)*ET2 
     >            + TEMFP2
      TEMFM2 = GAUSS(I1,I2,FNSMDGBL,XV(I2),DX,AER,RER,ER211,IRT)*ET2
     >            + TEMFM2 
      TEMG2 = GAUSS(I1,I2,RGG2DGBL,XV(I2),DX,AER1,RER1,ER2,IRT)*ET3
     >            + TEMG2

         ELSE


      TEMFP2=GAUSS(I1,I2,FNSPDGBL,XV(I2),XV(I2+1),AER,RER,ER21,IRT)*ET2
     >  + TEMFP2 
 
      TEMFM2=GAUSS(I1,I2,FNSMDGBL,XV(I2),XV(I2+1),AER,RER,ER211,IRT)*ET2
     >  + TEMFM2 
 
      TEMG2=GAUSS(I1,I2,RGG2DGBL,XV(I2),XV(I2+1),AER1,RER1,ER2,IRT)*ET3
     >  + TEMG2

        ENDIF

 911    CONTINUE

C
C     Save result of del/x..1 integration.
C

        TDGBLPP(I1) = TEMFP2
        TDGBLPM(I1) = TEMFM2
        TDGBLGG(I1) = TEMG2

C
C     reset variables!
C

        TEMFP2 = 0.0
        TEMFM2 = 0.0
        TEMG2  = 0.0


C
C     Does the integration of P(y,del/x) from del/x..1 for fixed x!
C

         IF (DX.GT.XV(NX-6)) THEN 

        TEMFM2 = 0  
	TEMFP2 = 0  
	TEMG2 = 0 

           NS = I1
        
        GOTO 201

        ELSE 

         DO 333 I3 = IV+1,NX-2 

            IF (XV(I3).GE.DX) THEN

               KK1 = I3
               GOTO 999

            ENDIF

 333     CONTINUE   

C
C     Set up array for interpolation!
C


 999    CALL SPLINED(XV,PNSPC,NX-2,KK1,I1,1,DX,DPP)
        CALL SPLINED(XV,PNSMC,NX-2,KK1,I1,1,DX,DPM)
        CALL SPLINED(XV,GG2C,NX-2,KK1,I1,1,DX,DGG)

         DO 211 I2 = KK1,NX-2

            IF (I2.EQ.KK1) THEN 

C
C     Integration from d/x..x_i2
C

        TEMFP2 = GAUSS(I1,I2-1,FNSPC,DX,XV(I2),AER,RER,ER21,IRT)*ET2 
 
	TEMFM2 = GAUSS(I1,I2-1,FNSMC,DX,XV(I2),AER,RER,ER211,IRT)*ET2
 
	TEMG2 = GAUSS(I1,I2-1,RGG2C,DX,XV(I2),AER1,RER1,ER2,IRT)*ET1


         ELSE

        
        TEMFP2=GAUSS(I1,I2,FNSPC,XV(I2),XV(I2+1),AER,RER,ER21,IRT)*ET2
     >  + TEMFP2 
 
	TEMFM2=GAUSS(I1,I2,FNSMC,XV(I2),XV(I2+1),AER,RER,ER211,IRT)*ET2
     >  + TEMFM2 
 
	TEMG2=GAUSS(I1,I2,RGG2C,XV(I2),XV(I2+1),AER1,RER1,ER2,IRT)*ET1
     >  + TEMG2

        ENDIF

 211    CONTINUE

        ENDIF

C
C     Save result of del/x..1 integration.
C

 201    T1FPD(I1) = TEMFP2
        T1FMD(I1) = TEMFM2
        T1GGD(I1) = TEMG2

C
C     reset variables!
C

        TEMFP2 = 0.0
        TEMFM2 = 0.0
        TEMG2  = 0.0

 
C 
C	Does the integration of - P(y,y*del/x) from x to 1 for fixed x! 
C 
        IF (I1.GE.NX-6) THEN 
 
	TEMFP22 = 0 
	TEMFM22 = 0 
	TEMG22 = 0 

	GOTO 200 
 
	ELSE 

        CALL SPLINED(XV,PNSP,NX-1,I1,IN,IV,DEL,DPP1)
        CALL SPLINED(XV,PNSM,NX-1,I1,IN,IV,DEL,DPM1)
        CALL SPLINED(XV,GG2,NX-1,I1,IN,IV,DEL,DGG1)

	DO 20 I2 = I1, NX-2 

	IT = I2 - IV 

	TEMFP22=GAUSS(IN,I2,FNSP,XV(I2),XV(I2+1),AER,RER,ER21,IRT)*ET2
     >  +TEMFP22 

	TEMFM22=GAUSS(IN,I2,FNSM,XV(I2),XV(I2+1),AER,RER,ER211,IRT)*ET2
     >  +TEMFM22 

	TEMG22=GAUSS(IN,I2,RGG2,XV(I2),XV(I2+1),AER1,RER1,ER2,IRT)*ET1
     >  +TEMG22 

   20 CONTINUE 
 
	ENDIF 
	
 
C
C save the result of x..1!
C

 200    TFPD(I1) = TEMFP22 
	TFMD(I1) = TEMFM22 
	TGGD(I1) = TEMG22 

C 
C Reset the variables for the inner loop! 
C 
 
	 
	TEMFP22 = 0 
	TEMFM22 = 0 
	TEMG22 = 0 
 
 21 	CONTINUE 

C 
C	extrapolate the last 7 bins with 10th. order polynomial. 
C

        KXL = 7

 	DO 44 IB = NX-KXL, NX-1 
 	 
 	DO 45 IA = 1,KXL 
 	X2(IA) = XV(IB-IA) 
  45	CONTINUE 
  	XX = XV(IB) 
 
  	DO 46 IC = 1,KXL 
  	FX(IC) = TFPD(IB-IC) 
  46	CONTINUE 
 
  	CALL RATINT(X2,FX,KXL,XX,TEM,ERR) 
        TFPD(IB) = TEM 
  	 
  	DO 47 IC = 1,KXL 
  	FX(IC) = TFMD(IB-IC) 
  47	CONTINUE 
 
  	CALL RATINT(X2,FX,KXL,XX,TEM,ERR) 
        TFMD(IB) = TEM
 
  	DO 48 IC = 1,KXL 
  	FX(IC) = TGGD(IB-IC)

  48	CONTINUE 
 
  	CALL RATINT(X2,FX,KXL,XX,TEM,ERR) 
        TGGD(IB) = TEM 

  44	CONTINUE 

C
C     Extrapolation of the first bins using a 10th order rational function
C     
        KL = 10

 	DO 441 IB = NS, IV+1,-1 
 	 
 	DO 451 IA = 1,KL 

        X2(IA) = 1D0*IA
 451    CONTINUE 
 
        XX = 0D0
  	DO 461 IC = 1,KL 
  	FX(IC) = T1FPD(IB+IC) 
 461    CONTINUE 
 
  	CALL RATINT(X2,FX,KL,XX,TEM,ERR) 
        T1FPD(IB) = TEM 

  	DO 471 IC = 1,KL 
  	FX(IC) = T1FMD(IB+IC) 
 471    CONTINUE 
 
  	CALL RATINT(X2,FX,KL,XX,TEM,ERR) 
        T1FMD(IB) = TEM

  	DO 481 IC = 1,KL 
  	FX(IC) = T1GGD(IB+IC)
 481    CONTINUE 
 
  	CALL RATINT(X2,FX,KL,XX,TEM,ERR) 
        T1GGD(IB) = TEM 

 441    CONTINUE 


C
C     Assemble final result of the counterterm integrals in the DGLAP region. 
C
	DO 49 I = IV+1,NX-1

        DXX = (XV(NX)-XV(NX-1))/2.
           
        ANSP(I) = TFPD(I) + T1FPD(I) + TDGBLPP(I)
     >  + (PNSPC(I,NX-1) + PNSP(NX-1-IV,I-IV))*DXX 
        ANSM(I) = TFMD(I) + T1FMD(I) + TDGBLPM(I)
     >  + (PNSMC(I,NX-1) + PNSM(NX-1-IV,I-IV))*DXX
        AGG2(I) = TGGD(I) + T1GGD(I) + TDGBLGG(I)
     >  + (GG2C(I,NX-1) + GG2(NX-1-IV,I-IV))*DXX

 49	Continue

        DO 50 IP = NX-2,NX-1

        DO 51 IU = 1,KL
         X2(IU) = XV(IP-IU)
         FX(IU) = ANSP(IP-IU)
         FX1(IU) = ANSM(IP-IU)
         FX2(IU) = AGG2(IP-IU)
 51     CONTINUE

        XX = XV(IP)

        CALL RATINT(X2,FX,KL,XX,TEM,ERR)

        ANSP(IP) = TEM

        CALL RATINT(X2,FX1,KL,XX,TEM,ERR)

        ANSM(IP) = TEM

        CALL RATINT(X2,FX2,KL,XX,TEM,ERR)

        AGG2(IP) = TEM

 50     CONTINUE


C 
C	reset temp. variables to 0 
C 
	TEMF1 = 0 
	TEMF11 = 0 
	TEMG1 = 0 
	TEMG11 = 0 
	TEMFM2 = 0 
	TEMFM22 = 0 
	TEMFP2 = 0 
	TEMFP22 = 0 
	TEMG2 = 0 
	TEMG22 = 0 

C 
C	Counterterms for ERBL region. First the integration between 0..x. 
C 
 
       
 
	DO 15 IX = 2,IV-1 

           IF (IX.GT.6) THEN

        CALL SPLINED1(XV,PNSPBLC1,IX,1,IX,0,DEL,DPPB)
        CALL SPLINED1(XV,PNSMBLC1,IX,1,IX,0,DEL,DPMB)
        CALL SPLINED1(XV,GG2BLC1,IX,1,IX,0,DEL,DGGB)

        ENDIF

	DO 16 IY = 1, IX-1 
 
	IF (IX.GT.6) THEN 

	TEMFP22=GAUSS(IX,IY,FNSPB,XV(IY),XV(IY+1),AER3,RER3,ER21,IRT)
     >  *ET2+TEMFP22 

	TEMFM22=GAUSS(IX,IY,FNSMB,XV(IY),XV(IY+1),AER3,RER3,ER211,IRT)
     >  *ET2+TEMFM22 

	TEMG22 = GAUSS(IX,IY,RGG2B,XV(IY),XV(IY+1),AER2,RER2,ER2,IRT)
     >  *ET3+TEMG22 

	ELSEIF (IX.LE.6) THEN 
 
	TEMFP22 = 0 
	TEMFM22 = 0 
	TEMG22 = 0 

 
	ENDIF 
 
 16	CONTINUE 

        TFP(IX) = TEMFP22 
	TFM(IX) = TEMFM22 
	TGG(IX) = TEMG22 

C 
C	now integral from x..del! 
C	
 	
	TEMFP2 = 0.0
	TEMFM2 = 0.0
	TEMG2 =  0.0

        IF (IX.LT.IV-5) THEN


        CALL SPLINED1(XV,PNSPBLCX,IV,IX,IX,0,DEL,DPPB1)
        CALL SPLINED1(XV,PNSMBLCX,IV,IX,IX,0,DEL,DPMB1)
        CALL SPLINED1(XV,GG2BLCX,IV,IX,IX,0,DEL,DGGB1)

        ENDIF

	DO 17 IY1 = IX,IV 
 
	IF (IX.LT.IV-5 .AND. IY1.LT.IV) THEN 

      TEMFP2=GAUSS(IX,IY1,FNSPB1,XV(IY1),XV(IY1+1),AER3,RER3,ER21,IRT)
     >       *ET2+TEMFP2 

      TEMFM2=GAUSS(IX,IY1,FNSMB1,XV(IY1),XV(IY1+1),AER3,RER3,ER211,IRT)
     >       *ET2+TEMFM2 

      TEMG2=GAUSS(IX,IY1,RGG2B1,XV(IY1),XV(IY1+1),AER2,RER2,ER2,IRT)
     >       *ET3+TEMG2  

	ELSEIF (IX.GE.IV-5) THEN 
 
	TEMFP2 = 0 
	TEMFM2 = 0 
	TEMG2 = 0 
	
	ENDIF 
 
 17	CONTINUE 


	T1FM(IX) = TEMFM2
	T1FP(IX) = TEMFP2
	T1GG(IX) = TEMG2

	TEMF1 = 0 
	TEMF11 = 0 
	TEMG1 = 0 
	TEMG11 = 0 
	TEMFM2 = 0 
	TEMFM22 = 0 
	TEMFP2 = 0 
	TEMFP22 = 0 
	TEMG2 = 0 
	TEMG22 = 0 

 15	CONTINUE

C 
C	extrapolate with 10th. order polynomial the first 5 bins. 
C 
 
        KT = 10
        KT1 = 9
	DO 1 IE = 6,2,-1 
	 
	DO 2 ID = 1,KT 
	X2(ID) = XV(IE+ID) 
 2	CONTINUE 
 	XX = XV(IE) 
 	 
 	DO 3 ID = 1,KT1
 	FX(ID) = TFP(IE+ID) 
 3	CONTINUE 
 
 	CALL RATINT(X2,FX,KT1,XX,TEM,ERR) 

 	TFP(IE) = TEM 
 	 
 	DO 4 ID = 1,KT1
 	FX(ID) = TFM(IE+ID) 
 4	CONTINUE 
 
 	CALL RATINT(X2,FX,KT1,XX,TEM,ERR) 
 	 
 	TFM(IE) = TEM

 	DO 5 ID = 1,KT
 	FX(ID) = TGG(IE+ID) 
 5	CONTINUE 
 
 	CALL RATINT(X2,FX,KT,XX,TEM,ERR) 
 	 
 	TGG(IE) = TEM 

 1	CONTINUE 

C 
C	extrapolate with 10th. order polynomial the last 6 bins. 
C 
        KT = 10

	DO 71 IE = IV-5,IV-1 
	 
	DO 72 ID = 1,KT
	X2(ID) = XV(IE-ID) 
 72	CONTINUE 
 	XX = XV(IE) 
 	 
 	DO 73 ID = 1,KT
 	FX(ID) = T1FP(IE-ID) 
 73	CONTINUE 
 
 	CALL RATINT(X2,FX,KT,XX,TEM,ERR) 
 	 
 	T1FP(IE) = TEM 

 	DO 74 ID = 1,KT
 	FX(ID) = T1FM(IE-ID) 
 74	CONTINUE 
 
 	CALL RATINT(X2,FX,KT,XX,TEM,ERR) 
 	 
 	T1FM(IE) = TEM 

 	DO 75 ID = 1,KT
 	FX(ID) = T1GG(IE-ID) 
 75	CONTINUE 
 
 	CALL RATINT(X2,FX,KT,XX,TEM,ERR) 
 	 
 	T1GG(IE) = TEM 
c        print *,"GG",IE,TEM,ERR

 71	CONTINUE 

        DO 18 J = 2,IV-1

	ANSM(J) = (TFM(J) + T1FM(J))*DIM 
	ANSP(J) = (TFP(J) + T1FP(J))*DIM 
	AGG2(J) = (TGG(J) + T1GG(J))*DIM 

 18	CONTINUE

        JK = 1
           KT = 5

      DO 811 IE = JK,1,-1
	 
	DO 821 ID = 1,KT
	X2(ID) = XV(IV-ID) 
 821  CONTINUE 
 	XX = XV(IV) 
 	 
 	DO 831 ID = 1,KT
 	FX(ID) = ANSP(IV-ID) 
 831	CONTINUE 
 
 	CALL POLINT(X2,FX,KT,XX,TEM,ERR) 
 	 
 	ANSP(IV) = TEM 

 	DO 841 ID = 1,KT
 	FX(ID) = ANSM(IV-ID) 
 841	CONTINUE 
 
 	CALL POLINT(X2,FX,KT,XX,TEM,ERR) 
 	 
 	ANSM(IV) = TEM 

 	DO 851 ID = 1,KT
 	FX(ID) = AGG2(IV-ID) 
 851  CONTINUE 
 
 	CALL POLINT(X2,FX,KT,XX,TEM,ERR) 
 	 
 	AGG2(IV) = TEM 

 811  CONTINUE  

      ANSM(1) = ANSM(IV)
      ANSP(1) = ANSM(IV)
      AGG2(1) = AGG2(IV)


      RETURN 
C                        **************************** 
      END 
 
      SUBROUTINE HQRK (NX, TT, NQRK, Y, F) 
C 
C       Subroutine to compute the (heavy) quark distribution from (given) 
C       gluon distribution, as the result of a change in the renormalization 
C       scheme (from MS-bar to BPH) as the threshold for quark flavor Nqrk 
C       is crossed. 
C 
C       Nx is the number of mesh-points, Tt is the Log Q variable. 
C       Y is the input g-distribution function defined on the mesh points 
C       F is the outpur Qrk-distribution function defined on the same pts. 
C 
C       If the threshold is chosen at the 'natural boundary' Mu = Mass(qrk), 
C       then there is no renormalization and this routine returns zero. 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
C 
      DIMENSION Y(NX), F(NX) 
C      DIMENSION W0(MXX), W1(MXX), W2(MXX), WH(MXX), WM(MXX) 
C 
C                                    Returns zero, assuming 'natural boundary'. 
      IF (NX .GT. 1) GOTO 11 
 
C      Q = EXP(TT) 
C      AL = ALPI(Q) 
C      AMS = AMASS(NQRK) 
C      AMH = AMHATF(NQRK) 
C      FAC = 2.* LOG (AMS / AMH) 
C      FAB = AL / 4. 
C 
C      CALL INTEGR (NX, 0, Y, W0, IR1) 
C      CALL INTEGR (NX, 1, Y, W1, IR2) 
C      CALL INTEGR (NX, 2, Y, W2, IR3) 
C 
   11 CONTINUE 
      DO 230 IZ = 1, NX 
C                                                       Returns zero answer. 
        IF (NX .GT. 1) THEN 
        F(IZ) = 0 
        GOTO 230 
        EndIf 
C 
C      F(IZ) = FAB * ( FAC * W0(IZ) 
C     >            + 2.* (FAC -1.) * (-W1(IZ) + W2(IZ) ) ) 
C      F(Iz) = fab * ( Fac * (W0(Iz) - Y(Iz) / 2.) 
C     >            + 2.* (fac -1.) * (-W1(Iz) + W2(Iz) + Y(Iz) / 12.)) 
C 
  230 CONTINUE 
C 
      RETURN 
C                        **************************** 
      END 
 
C 
C $Id: dilog64.F,v 1.1.1.1 1996/04/01 15:02:05 mclareni Exp $ 
C 
C $Log: dilog64.F,v $ 
C Revision 1.1.1.1  1996/04/01 15:02:05  mclareni 
C Mathlib gen 
C 
C 
      FUNCTION SPENCE(X) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION C(0:19) 
 
      PARAMETER (Z1 = 1, HF = Z1/2) 
      PARAMETER (PI = 3.14159 26535 89793 24D0) 
      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12) 
 
      DATA C( 0) / 0.42996 69356 08136 97D0/ 
      DATA C( 1) / 0.40975 98753 30771 05D0/ 
      DATA C( 2) /-0.01858 84366 50145 92D0/ 
      DATA C( 3) / 0.00145 75108 40622 68D0/ 
      DATA C( 4) /-0.00014 30418 44423 40D0/ 
      DATA C( 5) / 0.00001 58841 55418 80D0/ 
      DATA C( 6) /-0.00000 19078 49593 87D0/ 
      DATA C( 7) / 0.00000 02419 51808 54D0/ 
      DATA C( 8) /-0.00000 00319 33412 74D0/ 
      DATA C( 9) / 0.00000 00043 45450 63D0/ 
      DATA C(10) /-0.00000 00006 05784 80D0/ 
      DATA C(11) / 0.00000 00000 86120 98D0/ 
      DATA C(12) /-0.00000 00000 12443 32D0/ 
      DATA C(13) / 0.00000 00000 01822 56D0/ 
      DATA C(14) /-0.00000 00000 00270 07D0/ 
      DATA C(15) / 0.00000 00000 00040 42D0/ 
      DATA C(16) /-0.00000 00000 00006 10D0/ 
      DATA C(17) / 0.00000 00000 00000 93D0/ 
      DATA C(18) /-0.00000 00000 00000 14D0/ 
      DATA C(19) /+0.00000 00000 00000 02D0/ 
 
      IF(X .EQ. 1) THEN 
       H=PI6 
      ELSEIF(X .EQ. -1) THEN 
       H=-PI12 
      ELSE 
       T=-X 
       IF(T .LE. -2) THEN 
        Y=-1/(1+T) 
        S=1 
        A=-PI3+HF*(LOG(-T)**2-LOG(1+1/T)**2) 
       ELSEIF(T .LT. -1) THEN 
        Y=-1-T 
        S=-1 
        A=LOG(-T) 
        A=-PI6+A*(A+LOG(1+1/T)) 
       ELSE IF(T .LE. -HF) THEN 
        Y=-(1+T)/T 
        S=1 
        A=LOG(-T) 
        A=-PI6+A*(-HF*A+LOG(1+T)) 
       ELSE IF(T .LT. 0) THEN 
        Y=-T/(1+T) 
        S=-1 
        A=HF*LOG(1+T)**2 
       ELSE IF(T .LE. 1) THEN 
        Y=T 
        S=1 
        A=0 
       ELSE 
        Y=1/T 
        S=-1 
        A=PI6+HF*LOG(T)**2 
       ENDIF 
       H=Y+Y-1 
       ALFA=H+H 
       B1=0 
       B2=0 
       DO 1 I = 19,0,-1 
 
       B0=C(I)+ALFA*B1-B2 
       B2=B1 
    1  B1=B0 
       H=-(S*(B0-H*B2)+A) 
      ENDIF 
      SPENCE=H 
      RETURN 
      END 

C 
       SUBROUTINE QARRAY (NINI) 
 
 
C ==================================================================== 
C      SUBROUTINE QARRAY (NINI) 
C 
C       Given Qini, Qmax, and Nt, this routine go through the various 
C       flavor thresholds; determines the step-sizes in-between each pair 
C       of thresholds, computes a new NT, if necessary; and returns: 
C        the variable NINI, NFMX as arguments;  
C        the (NT+1)-dim Q- and T=lnlnQ/Lamda arrays in the common block QARAY1; 
C        the Neff-dependent variable (see below) in QARAY2; and  
C        KF in the common block EvlPac 
C 
C       Note in particular that the value of NT may increase by 1 or 2 in a 
C       call to this routine because of the logistics of setting up the steps 
C       for the flavor thresholds. 
C 
C      For given Neff, the variables TLN, DTN, NTL and NTN has the following 
C      meaning: 
 
C     Flavor thresholds 
C     Nini        Nf                                          Nf+1  ....NFmx 
C     |.. ........|<- DTN ->|<- DTN ->| .......................|..........| 
C                  <------------  NTL(Nf) steps  -------------> 
C     0         NTN(Nf)                                    NTN(Nf+1)     NT 
C               TLN(Nf)       (TLN=lnln(Q/Lamda))          TNL(Nf+1) 
C    TV(0).... TV(NTN(Nf)).....TV(NTN(Nf)+2).............TV(NTN(Nf+1))..TV(NT) 
C    QV(0).... QV(NTN(Nf)).....QV(NTN(Nf)+2).............QV(NTN(Nf+1))..QV(NT) 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (MXPN = MXF * 2 + 2) 
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN) 
C 
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG 
      COMMON / QARAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF) 
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx 
 
      NCNT = 0 
 
      IF (NT .GE. mxq) NT = mxq - 1 
 
C                                        Set up overall t=log(Q) parameters 
 
      S = LOG(QINI/AL) 
      TINI = LOG(S) 
      S = LOG(QMAX/AL) 
      TMAX = LOG(S) 
C                              Nt is the provisional # of mesh-points in t; 
C                                    Dt0 is the approximate increment in T, 
C                           Actual values of (Nt, Dt) are determined later. 
    1 DT0 = (TMAX - TINI) / float(NT) 
C                            Determine effective Nflv at Qini and at Qmax 
      NINI = NFL(QINI) 
      NFMX = NFL(QMAX) 
      Call ParQcd (2, 'ORDER', Ord, Ir) 
      Call ParQcd (2, 'ALAM', Al0, Ir) 
      Call ParQcd (2, 'NFL', Afl0, Ir) 
C                                      Set the total number of Quark flavors 
C                                      in the QCD package to NfMx 
      AFL = NfMx 
      Call ParQcd (1, 'NFL', AFL, Ir) 
C                                      and restore the QCD coupling 
      Iordr = Nint (Ord) 
      Ifl0  = Nint (Afl0) 
      Call SetLam (Ifl0, Al0, Iordr) 
C 
C      Q2 evolution is carried out in stages separated by flavor thresholds 
C 
      NG = NFMX - NINI + 1 
C                                                  ----------------------- 
C                                      Determine the threshold points and 
C                                     Set up detailed  t- mesh structure 
      QIN  = QINI 
      QOUT = QINI 
      S = LOG(QIN/AL) 
      TIN  = LOG(S) 
      TLN(1) = TIN 
      NTL(1)  = 0 
      QV(0) = QINI 
      TV(0) = Tin 
C 
      DO 20 NEFF = NINI, NFMX 
C 
        ICNT = NEFF - NINI + 1 
        IF (NEFF .LT. NFMX) THEN 
          THRN = AMHATF (NEFF + 1) 
          QOUN = MIN (QMAX, THRN) 
        Else 
          QOUN = QMAX 
        EndIf 
C 
        IF (QOUN-QOUT .LE. 0.0001) THEN 
          DT   = 0 
          NITR = 0 
        Else 
          QOUT = QOUN 
          S = LOG(QOUT/AL) 
          TOUT = LOG(S) 
          TEM = TOUT - TIN 
C                               Nitr = Number of iterations in this stage 
          NITR = INT (TEM / DT0) + 1 
          DT  = TEM / NITR 
        EndIf 
C                                       Save book-keeping data on array 
        DTN (ICNT) = DT 
        NTN (ICNT) = NITR 
        TLN (ICNT) = TIN 
        NTL (ICNT+1) = NTL(ICNT) + NITR 
C 
C                      QV is the physical Q-value for this lattice point. 
        IF (NITR .NE. 0) THEN 
        DO 205 I = 1, NITR 
           TV (NTL(ICNT)+I) = TIN + DT * I 
           S = EXP (TV(NTL(ICNT)+I)) 
           QV (NTL(ICNT)+I) = AL * EXP (S) 
  205   CONTINUE 
        EndIf 
C                                        Initialize for next iteration 
        QIN = QOUT 
        TIN = TOUT 
C 
   20 CONTINUE 
C            Redefine Nt, as the actual number of t-points; may have been 
C                             increased by 1 or 2 by the above algorithm 
C                       If NT > MXQ, reduce NT until it is within bounds 
      NCNT = NCNT + 1 
      NTP = NTL (NG + 1) 
      ND  = NTP - NT 
 
      IF (NTP .GE. MXQ) THEN 
         NT = MXQ - ND - NCNT 
         GOTO 1 
      EndIf 
      NT = NTP 
C 
      RETURN 
C                        **************************** 
      END 
 
 
      SUBROUTINE XARRAY 
C 
C       Given NX, this routine fills the x-arrays of the common blocks 
C       related to the x-variable for the use of various other routines. 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      LOGICAL LSTX 
C 
      PARAMETER (D0 = 0.0, D10=10.0) 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (MXPN = MXF * 2 + 2) 
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN) 
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2) 
C 
      Character Msg*80 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX) 
      COMMON / VARBAB / GB(NDG, NDH, MXX), H(NDH, MXX, M1:M2) 
      COMMON / HINTEC / GH(NDG, MXX) 
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX, DEL,IV 
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX)
      COMMON / DIFFD / DIFFDEL  
      
C 
      DIMENSION G1(NDG,NDH), G2(NDG,NDH), A(NDG) 
C 
      DATA F12, F22, F32 / 1D0, 1D0, 1D0 / 
      DATA (G1(I,NDH), G2(I,1), I=1,NDG) / 0.0,0.0,0.0,0.0,0.0,0.0 / 
      DATA PUNY / 1D-30 / 
C 
C                                                ---------------------------- 
C                                                   Change of variable X <--> Z 
C                        Range of Z is [0, 1] for X in [XM, 1] and I in [1, NX] 
C                                                 Use evenly spaced points in Z 
      IV = 1
      XV(0) = 0D0 
      DZ = 1D0 / (NX-1) 
      ZV(NX) = 1D0
      DO 10 I = 1, NX - 1 
         Z = DZ * (I-1) 
         ZV(I) = Z 

         IF(DEL.LT.0.1) then

            KX = NX/4

            ELSE

            KX = NX/2

         ENDIF

         IF (I.LE.KX) THEN

            X = XFRMZN(1,I,Z)

         DXTZ(I) = DXDZN(1,I,Z)

         IF (I.EQ.KX) DIFFDEL = DXDZN(2,I,Z)

         ELSE
       
            X = XFRMZN (2,I,Z) 
        
C                               DXDZ is the Jacobian for the change of variable 
         DXTZ(I) = DXDZN(2,I,Z)

         ENDIF

C
C	Determine the last counter position where x is smaller than del! 
C          
	IF (X.LE.DEL) THEN 
	 IV = I 
	ELSE 
	IV = IV 
	ENDIF 
C                                                      Fill x - arrays 
         XV (I)= X 

         XA(I, 1) = X 

         XA(I, 0) = LOG (X) 
         DO 20 L = L1, L2 
          IF (L .NE. 0 .AND. L .NE. 1)  XA(I, L) = X ** L 
   20    CONTINUE 
   10 CONTINUE 
 
C                                                         Fill last point, I=Nx 
         XV(NX) = 1D0 
         DXTZ(NX) = DXDZN(2,NX,1D0) 
         DO 21 L = L1, L2 
            XA (NX, L) = 1D0 
   21    CONTINUE 
         XA (NX, 0) = 0D0 
C                                                       Fill ELY = Log (1. - x) 
 
      DO 11 I = 1, NX-1 
         ELY(I) = LOG(1D0 - XV(I)) 
   11 CONTINUE 
C                     Log(1-x) is infinite at x=1, for the purpose of numerical 
C                       Calculations, the following extrapolated value is used 
C                       to avoid an artificial discontinuity. 
C 
       ELY(NX) = 3D0* ELY(NX-1) - 3D0* ELY(NX-2) + ELY(NX-3) 
C                                                  ---------------------------- 
C                                    Matrix elements for 2nd order calculations 
C 
C                                     1 . . . . . . . . . 
C                                     . 1 . . . . . . . . 
C                                     . . 1             . 
C                                     .     1   Z (X/Y) . 
C                                     .       1         . 
C                                     .         .       . 
C                                     .  X/Y      1     . 
C                                     .             1   . 
C                                     . . . . . . . . 1 . 
C                                     . . . . . . . . . 1 
      DO 17 IX = 1, NX 
      ZZ (IX, IX) = 1. 
      DO 17 IY = IX+1, NX 
         XY = XV(IX) / XV(IY) 
c         ZZ (IX, IY) = ZFRMX (XY) 
         ZZ (IX,IY) = 0D0
         ZZ (NX-IX+1, NX-IY+1) = XY 
   17 CONTINUE 
C                                                ------------------------------ 
C                                     Start of x - loop to compute ceefficients 
C                                   for integrals used in INTEG, AMOM, ... etc. 
      DO 30 I = 1, NX-1 
C                                                  "F" matrix {a(i)} --> {f(i)} 
      IF (I .NE. NX-1) THEN 
        F11 = 1D0/XV(I) 
        F21 = 1D0/XV(I+1) 
        F31 = 1D0/XV(I+2) 
        F13 = XV(I) 
        F23 = XV(I+1) 
        F33 = XV(I+2) 
C                                                        Determinant for matrix 
C 
        DET = F11*F22*F33 + F21*F32*F13 + F31*F12*F23 
     >      - F31*F22*F13 - F21*F12*F33 - F11*F32*F23 
        IF (ABS(DET) .LT. PUNY) THEN 
           Msg='Determinant close to zero; will be arbitrarily set to:' 
           CALL WARNR(IWRN, NWRT, Msg, 'DET', PUNY, D0, D0, 0) 
           DET = PUNY 
        EndIf 
C                                            Compute "G" matrix -- inverse of F 
C                                             G1 is only needed from I=2 and on 
        G2(1,2) = (F22*F33 - F23*F32) / DET 
        G2(1,3) = (F32*F13 - F33*F12) / DET 
        G2(1,4) = (F12*F23 - F13*F22) / DET 
C 
        G2(2,2) = (F23*F31 - F21*F33) / DET 
        G2(2,3) = (F33*F11 - F31*F13) / DET 
        G2(2,4) = (F13*F21 - F11*F23) / DET 
C 
        G2(3,2) = (F21*F32 - F22*F31) / DET 
        G2(3,3) = (F31*F12 - F32*F11) / DET 
        G2(3,4) = (F11*F22 - F12*F21) / DET 
C                                               Compute coefficients for HINTEG 
        B2 = LOG (XV(I+2)/XV(I)) 
        B3 = XV(I) * (B2 - 1.) + XV(I+2) 
        GH (1,I) = B2 * G2 (2,2) + B3 * G2 (3,2) 
        GH (2,I) = B2 * G2 (2,3) + B3 * G2 (3,3) 
        GH (3,I) = B2 * G2 (2,4) + B3 * G2 (3,4) 
      EndIf 
C                                         "G-bar" is the "average" of G1 and G2 
        DO 51 J = 1, NDH 
           DO 52 L = 1, NDG 
C                                                                First interval 
              IF     (I .EQ. 1) THEN 
                 GB(L,J,I) = G2(L,J) 
C                                                                 last interval 
              ElseIF (I .EQ. NX-1) THEN 
                 GB(L,J,I) = G1(L,J) 
C                                                        intermidiate intervals 
              Else 
                 GB(L,J,I) = (G1(L,J) + G2(L,J)) / 2D0 
              EndIf 
   52      CONTINUE 
   51   CONTINUE 
C                                                            Compute "A" matrix 
        DO 35 MM = M1, M2 
           DO 40 K = 1, NDG 
             KK = K + MM - 2 
             IF (KK .EQ. 0) THEN 
               A(K) = XA(I+1, 0) - XA(I, 0) 
             Else 
               A(K) = (XA(I+1, KK) - XA(I, KK)) / DBLE(KK) 
             EndIf 
   40      CONTINUE 
C                                                            Compute "H" matrix 
           DO 41 J = 1, NDH 
             TEM = 0 
             DO 43 L = 1, NDG 
               TEM = TEM + A(L) * GB(L,J,I) 
   43        CONTINUE 
             H(J,I,MM) = TEM 
   41      CONTINUE 
   35   CONTINUE 
C                                               ------------------------------ 
C                                                         Initialize G1 matrix 
      DO 42 J = 1, NDG 
        DO 44 L = 1, NDG 
           G1(L,J) = G2(L,J+1) 
   44 CONTINUE 
   42 CONTINUE 
C 
   30 CONTINUE 
C                                    End of x - loop to calculate coefficients 
C                                               ------------------------------ 
      LSTX = .TRUE. 
      RETURN 
C                        **************************** 
      END 
 
      SUBROUTINE NEWARRAY (NX, DEL, F, FF) 
 
C 
C     In this routine the input parton distribution f(y) will be approximated 
C     in the bin x_i, x_i+1 by making apolynomial approximation: 
C                    f(y)= a*y^2 + b*y +c 
C     where the coefficients a,b,c are determined by knowing f(y) in 
C     the points x_i-1, x_i, x_i+1. This relationship yields a 3x3 matrix 
C     which has then to be inverted to find a,b,c!  
C     The coeffiecints for the bins are 
C     stored in the array FF as FF(i) = a, FF(i+1) = b, FF(i+2) = c 
C     corresponding to the bin x_i, x_i+1.  
 
 
            IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C 
      CHARACTER MSG*80 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (MXPN = MXF * 2 + 2) 
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN) 
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2) 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX) 
      COMMON / VARBAB / GB(NDG, NDH, MXX), H(NDH, MXX, M1:M2) 
      DIMENSION F(NX), FF(3*MXX),XT(10),F1(10),F2(10),F3(10) 
C 
C     Initialize counter 
C  
            K = 1 
C 
C     Loop to calculate a,b,c except first bin! 
C 
      DO 10 I = NX-1, 1, -1 
C 
 
      IF (XA(I,1).GT.DEL .AND. XA(I-1,1).GT.DEL) THEN 
		 
      DET =1D0/((XA(I+1,1) - XA(I-1,1))*(XA(I+1,1) - XA(I,1))) 
      DET1 =1D0/((XA(I,1) - XA(I-1,1))*(XA(I+1,1) - XA(I,1))) 
      DET2 =1D0/((XA(I,1) - XA(I-1,1))*(XA(I+1,1) - XA(I-1,1))) 
C 
      FF(K) = F(I+1)*DET -F(I)*DET1 +F(I-1)*DET2 
C 
      FF(K+1) = - (XA(I-1,1) + XA(I,1))*F(I+1)*DET  
     >          + (XA(I+1,1) + XA(I-1,1))*F(I)*DET1  
     >          - (XA(I,1) + XA(I+1,1))*F(I-1)*DET2 
C 
      FF(K+2) = + XA(I,1)*XA(I-1,1)*F(I+1)*DET  
     >          - XA(I-1,1)*XA(I+1,1)*F(I)*DET1  
     >          + XA(I+1,1)*XA(I,1)*F(I-1)*DET2 
	 
 
        ELSEIF (XA(I,1).GT.DEL .AND. XA(I-1,1).LE.DEL) THEN


            DET =1D0/((XA(I+2,1) - XA(I,1))*(XA(I+2,1) - XA(I+1,1))) 
            DET1 =1D0/((XA(I+1,1) - XA(I,1))*(XA(I+2,1) - XA(I+1,1))) 
            DET2 =1D0/((XA(I+1,1) - XA(I,1))*(XA(I+2,1) - XA(I,1))) 
C 
      FF(K) = F(I+2)*DET -F(I+1)*DET1 +F(I)*DET2 
C 
      FF(K+1) = - (XA(I,1) + XA(I+1,1))*F(I+2)*DET  
     >          + (XA(I+2,1) + XA(I,1))*F(I+1)*DET1  
     >          - (XA(I+1,1) + XA(I+2,1))*F(I)*DET2 
C 
      FF(K+2) = + XA(I+1,1)*XA(I,1)*F(I+2)*DET  
     >          - XA(I,1)*XA(I+2,1)*F(I+1)*DET1  
     >          + XA(I+2,1)*XA(I+1,1)*F(I)*DET2 


        ELSEIF (XA(I,1).LE.DEL .AND. XA(I+1,1).GT.DEL) THEN

            DET =1D0/((XA(I+1,1) - XA(I-1,1))*(XA(I+1,1) - XA(I,1))) 
            DET1 =1D0/((XA(I,1) - XA(I-1,1))*(XA(I+1,1) - XA(I,1))) 
            DET2 =1D0/((XA(I,1) - XA(I-1,1))*(XA(I+1,1) - XA(I-1,1))) 
C 
      FF(K) = F(I+1)*DET -F(I)*DET1 +F(I-1)*DET2 
C 
      FF(K+1) = - (XA(I-1,1) + XA(I,1))*F(I+1)*DET  
     >          + (XA(I+1,1) + XA(I-1,1))*F(I)*DET1  
     >          - (XA(I,1) + XA(I+1,1))*F(I-1)*DET2 
C 
      FF(K+2) = + XA(I,1)*XA(I-1,1)*F(I+1)*DET  
     >          - XA(I-1,1)*XA(I+1,1)*F(I)*DET1  
     >          + XA(I+1,1)*XA(I,1)*F(I-1)*DET2 


      ELSEIF (XA(I,1).LT.DEL .AND. XA(I+1,1).LE.DEL .AND. I.GT.1) THEN

      DET =1D0/((XA(I+1,1) - XA(I-1,1))*(XA(I+1,1) - XA(I,1))) 
      DET1 =1D0/((XA(I,1) - XA(I-1,1))*(XA(I+1,1) - XA(I,1))) 
      DET2 =1D0/((XA(I,1) - XA(I-1,1))*(XA(I+1,1) - XA(I-1,1))) 
C 
      FF(K) = F(I+1)*DET -F(I)*DET1 +F(I-1)*DET2 
C 
      FF(K+1) = - (XA(I-1,1) + XA(I,1))*F(I+1)*DET  
     >          + (XA(I+1,1) + XA(I-1,1))*F(I)*DET1  
     >          - (XA(I,1) + XA(I+1,1))*F(I-1)*DET2 
C 
      FF(K+2) = + XA(I,1)*XA(I-1,1)*F(I+1)*DET  
     >          - XA(I-1,1)*XA(I+1,1)*F(I)*DET1  
     >          + XA(I+1,1)*XA(I,1)*F(I-1)*DET2 


	ELSEIF (I.EQ.1 .AND. XA(1,1).LT.DEL) THEN 
C	FIX last bin in BL Region 

           DET =1D0/((XA(I+2,1) - XA(I,1))*(XA(I+2,1) - XA(I+1,1))) 
            DET1 =1D0/((XA(I+1,1) - XA(I,1))*(XA(I+2,1) - XA(I+1,1))) 
            DET2 =1D0/((XA(I+1,1) - XA(I,1))*(XA(I+2,1) - XA(I,1))) 
C 
      FF(K) = F(I+2)*DET -F(I+1)*DET1 +F(I)*DET2 
C 
      FF(K+1) = - (XA(I,1) + XA(I+1,1))*F(I+2)*DET  
     >          + (XA(I+2,1) + XA(I,1))*F(I+1)*DET1  
     >          - (XA(I+1,1) + XA(I+2,1))*F(I)*DET2 
C 
      FF(K+2) = + XA(I+1,1)*XA(I,1)*F(I+2)*DET  
     >          - XA(I,1)*XA(I+2,1)*F(I+1)*DET1  
     >          + XA(I+2,1)*XA(I+1,1)*F(I)*DET2 


	ELSEIF (I.EQ.1 .AND. XA(1,1).GT.DEL) THEN 
C	Fix the last bin in DGLAP region if evolution only in DGLAP region! 
 
        DET =1D0/((XA(I+2,1) - XA(I,1))*(XA(I+2,1) - XA(I+1,1))) 
            DET1 =1D0/((XA(I+1,1) - XA(I,1))*(XA(I+2,1) - XA(I+1,1))) 
            DET2 =1D0/((XA(I+1,1) - XA(I,1))*(XA(I+2,1) - XA(I,1))) 
C 
        FF(K) = F(I+2)*DET -F(I+1)*DET1 +F(I)*DET2 
C 
        FF(K+1) = - (XA(I,1) + XA(I+1,1))*F(I+2)*DET  
     >          + (XA(I+2,1) + XA(I,1))*F(I+1)*DET1  
     >          - (XA(I+1,1) + XA(I+2,1))*F(I)*DET2 
C 
        FF(K+2) = + XA(I+1,1)*XA(I,1)*F(I+2)*DET  
     >          - XA(I,1)*XA(I+2,1)*F(I+1)*DET1  
     >          + XA(I+2,1)*XA(I+1,1)*F(I)*DET2 

 
 
	ENDIF 

       K = K + 3 
       
   10 CONTINUE 

         RETURN 

 
C         ***************************************************** 
 
         END 

      SUBROUTINE INTEGR (NX, M, F,   G, IR) 
C 
C     Computes (x/y) ** M * F(y)dy/y integrated over the range [x, 1]; 
C              Result is G. 
C     Integand function F must be defined on an array of size NX which 
C              covers the range [0, 1] of the variables x and y. 
C     The output function G is returned on the same array. 
C     The use of integration variable z defined by 
C 
C                  y = f(x) 
C 
C     where f(x) can be any monotonic function. 
C 
C     IR is an error return code: IR = 1   --  NX out of range 
C                                    = 2   --  M  out of range 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      CHARACTER MSG*80 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (MXPN = MXF * 2 + 2) 
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN) 
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2) 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX) 
      COMMON / VARBAB / GB(NDG, NDH, MXX), H(NDH, MXX, M1:M2) 
C 
      DIMENSION   F(NX), G(NX)
C 
      DATA IWRN1, IWRN2 / 0, 0 / 
C 
      IRR = 0 
C 
      IF (NX .LT. 1 .OR. XA(NX-1,1) .EQ. 0D0) THEN 
        MSG = 'NX out of range in INTEGR call' 
        CALL WARNI (IWRN1, NWRT, MSG, 'NX', NX, 0, MXX, 0) 
        IRR = 1 
      EndIf 
C 
      IF (M .LT. M1 .OR. M .GT. M2) THEN 
        MSG ='Exponent M out of range in INTEGR' 
        CALL WARNI (IWRN2, NWRT, MSG, 'M', M, M1, M2, 1) 
        IRR = 2 
      EndIf 
C                                                                         NX 
      G(NX) = 0D0 
C                                                                       NX - 1 
      TEM = H(1, NX-1, -M) * F(NX-2) + H(2, NX-1, -M) * F(NX-1) 
     >    + H(3, NX-1, -M) * F(NX) 
      IF (M .EQ. 0) THEN 
         G(NX-1) = TEM 
      Else 
 
         G(NX-1) = TEM * XA(NX-1, M) 
      EndIf 
C                                                                     NX-2 : 2 
      DO 10 I = NX-2, 2, -1 
         TEM = TEM + H(1,I,-M)*F(I-1) + H(2,I,-M)*F(I) 
     >             + H(3,I,-M)*F(I+1) + H(4,I,-M)*F(I+2) 
         IF (M .EQ. 0) THEN 
            G(I) = TEM 
         Else 
            G(I) = TEM * XA(I, M) 
         EndIf 
   10 CONTINUE 
C                                                                            1 
      TEM = TEM + H(2,1,-M)*F(1) + H(3,1,-M)*F(2) + H(4,1,-M)*F(3) 
      IF (M .EQ. 0) THEN 
         G(1) = TEM 
      Else 
         G(1) = TEM * XA(1, M) 
      EndIf 
 
      IR = IRR 
C 
      RETURN 
C                        **************************** 
      END 
 
      SUBROUTINE NINTEGR (NX, DEL, M, FF, F, G, IR) 
C 
C     Computes the integral of the more complicated pieces appearing 
C     in the kernels of the ND case as found in .... over the range  
C     [x, 1];The integrals are computed analytically within a bin and then  
C     summed to give the final answer. 
C     The input function F(y) has been replaced by FF an array containing  
C     the constants a,b,c from the quadratic approximation of F(y). 
C     The Result is G. 
C     Integand function F must be defined on an array of size NX which 
C     covers the range [0, 1] of the variables x and y. 
C     The output function G is returned on the same array. 
C 
C     IR is an error return code: IR = 1   --  NX out of range 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      CHARACTER MSG*80 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (MXPN = MXF * 2 + 2) 
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN) 
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)           
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX) 
      COMMON / VARBAB / GB(NDG, NDH, MXX), H(NDH, MXX, M1:M2) 
C 
      DIMENSION   FF(3*MXX), G(NX), F(NX), X2(10),FX(10),FX1(10),FX2(10) 
C 
      DATA IWRN1, IWRN2 / 0, 0 / 
C 
      IRR = 0 
C 
      IF (NX .LT. 1 .OR. XA(NX-1,1) .EQ. 0D0) THEN 
        MSG = 'NX out of range in INTEGR call' 
        CALL WARNI (IWRN1, NWRT, MSG, 'NX', NX, 0, MXX, 0) 
        IRR = 1 
      EndIf 
 
C 
C     Difference between x_min and delta is used to determine which analytic 
C     expression for the integrals to use! T < 10^-10 is treated as  
C     x_min = delta i.e T = 0 to avoid accuracy errors in the computation. 
C     T > 0.98 x_min will be treated by using the analytic expression for 
C     the integrals expanded to O(delta). This case also allows delta = 0 
C     i.e the diagonal case.  
C 
       
      IF (DEL .EQ. 0.0) then       
          DIM = 0.0           
      ELSE        
           DIM = 1D0/DEL 
      ENDIF 

      ITC = 0

           DO 234 KP = 1,NX
              IF (XA(KP,1).LE.DEL) THEN
                 ITC = ITC + 1
                 ELSE
                    ITC = ITC
                    ENDIF

 234             CONTINUE


                T = XA(ITC,1) - DEL 
                T1 = 0.98*XA(ITC,1)  
              
	
C                                               NX 
      G(NX) = 0 
C                                                             Nx-1 -> 1 
       TEM = 0 
       TEM1 = 0 
       TEM2 = 0  
         K = 1 
         KK = 1     
      
C 
C     a = FF(K) b = FF(K+1) c = FF(K+2) from quadratic approximation 
C     of F(X) knowing the function at xi , x(i-1), x(i+1) 
C     The following are the analytical results of the integrals from 
C     the kernels. 
C  
      IF (M .EQ. 1) THEN 
 
      IF (T .GT. T1) THEN 
 
       DO 90 J = NX-1, 1, -1 
 
      TEM = TEM + (FF(K+1)+FF(K)*DEL)*(XA(J,-1) - XA(J+1,-1))  
     >  + FF(K)*(XA(J+1,0) - XA(J,0)) 
     >  + (FF(K+2) + FF(K+1)*DEL)*(XA(J,-2) - XA(J+1,-2))/2. 
     >  + FF(K+2)*DEL*(XA(J,-3) - XA(J+1,-3))/3. 
   
         G(J) = (XA(J,1)-DEL)*XA(J,1)*TEM   
 
            K = K + 3 
 
   90 CONTINUE     
 
      ELSE 
    
      DO 10 I = NX-1, 1, -1 
       
      IF (XA(I,1) .LE. DEL) then 
 
	GOTO 666 
 
	ELSE 
 
       A2 = XA(I+1,1)*DIM - 1D0 
       A1 = XA(I,1)*DIM - 1D0 
       B1 = DEL*XA(I,-1) 
       B2 = DEL*XA(I+1,-1) 
        
      TEM = TEM + FF(K+2)*(B2 - B1)  
     >  + (FF(K+2) + FF(K+1)*DEL)*(XA(I,0) - XA(I+1,0)) 
     >  + (FF(K+2) + FF(K+1)*DEL + FF(K)*DEL**2) 
     >  *(LOG(A2) - LOG(A1)) 
   
         G(I) = (A1+1D0)*(A1)*TEM   

            K = K + 3  
	ENDIF 
 
   10 CONTINUE     
 
      ENDIF    
 
      ENDIF  
   

C
      IF (M .EQ. 2) THEN


      IF (T .GT. T1) THEN 
   
      DO 91 JI = NX-1, 1, -1

      DO 92 IJ = NX-1, JI, -1

      TEM = TEM + (FF(K+1)*XA(JI,1) - 2.*FF(K)*XA(JI,2) + 2.*FF(K)
     > *XA(JI,1)*DEL)*(XA(IJ+1,0) - XA(IJ,0)) + FF(K)*XA(JI,1)
     > *(XA(IJ+1,1) - XA(IJ,1))
     > + (FF(K)*XA(JI,3) - 2.*FF(K+1)*XA(JI,2) + FF(K+2)*XA(JI,1)
     > + 2.*FF(K+1)*XA(JI,1)*DEL - 4.*FF(K)*XA(JI,2)*DEL)
     > *(XA(IJ,-1) - XA(IJ+1,-1)) + (FF(K+1)*XA(JI,3)/2. - FF(K+2)
     > *XA(JI,2) + FF(K+2)*XA(JI,1)*DEL - 2.*FF(K+1)*XA(JI,2)*DEL
     > + FF(K)*XA(JI,3)*DEL)*(XA(IJ,-2) - XA(IJ+1,-2))
     > + (FF(K+2)*XA(JI,3)/3. - 4.*FF(K+2)*XA(JI,2)*DEL/3.
     > + 2.*FF(K+1)*XA(JI,3)*DEL/3.)*(XA(IJ,-3) - XA(IJ+1,-3))
     > + FF(K+2)*XA(JI,3)*DEL*(XA(IJ,-4) - XA(IJ+1,-4))/2. 

      K = K + 3

 92   CONTINUE             
      
      G(JI) = TEM
       TEM = 0
         K = 1


 91   CONTINUE        

      ELSE

      DO 12 J = NX-1, 1, -1

	IF (XA(J,1) .LE. DEL) THEN

	GOTO 666

	ELSE

       AA1 = XA(J,1)*DIM

      DO 23 I = NX-1, J, -1

         COMP = DEL/XA(J,1)
         COMP1 = DEL/XA(I,1)

         IF (COMP .LT. 0.0005 .OR. COMP1.LT.0.0005) THEN

            TEM = TEM + (FF(K+1)*XA(J,1) - 2.*FF(K)*XA(J,2) + 2.*FF(K)
     > *XA(J,1)*DEL)*(XA(I+1,0) - XA(I,0)) + FF(K)*XA(J,1)
     > *(XA(I+1,1) - XA(I,1))
     > + (FF(K)*XA(J,3) - 2.*FF(K+1)*XA(J,2) + FF(K+2)*XA(J,1)
     > + 2.*FF(K+1)*XA(J,1)*DEL - 4.*FF(K)*XA(J,2)*DEL)
     > *(XA(I,-1) - XA(I+1,-1)) + (FF(K+1)*XA(J,3)/2. - FF(K+2)
     > *XA(J,2) + FF(K+2)*XA(J,1)*DEL - 2.*FF(K+1)*XA(J,2)*DEL
     > + FF(K)*XA(J,3)*DEL)*(XA(I,-2) - XA(I+1,-2))
     > + (FF(K+2)*XA(J,3)/3. - 4.*FF(K+2)*XA(J,2)*DEL/3.
     > + 2.*FF(K+1)*XA(J,3)*DEL/3.)*(XA(I,-3) - XA(I+1,-3))
     > + FF(K+2)*XA(J,3)*DEL*(XA(I,-4) - XA(I+1,-4))/2.

            ELSE

       A1 = XA(I,1)*DIM - 1D0
       A2 = XA(I+1,1)*DIM -1D0

      TEM = TEM + AA1*(FF(K+2) + FF(K+1)
     >  *DEL + FF(K)*DEL**2)*(1D0 - AA1)**2
     >  *(A1**-1 - A2**-1)
     >  + FF(K)*XA(J,1)*(XA(I+1,1) - XA(I,1))
     >  + FF(K+2)*(DEL*XA(I,-1) - DEL*XA(I+1,-1))*AA1**3
     >  + AA1**2*(FF(K+1)*XA(J,1) + 2D0*FF(K+2)*AA1
     >   - 2D0*FF(K+2))*(XA(I+1,0) - XA(I,0))
     >  + AA1*(1D0 - AA1)*(FF(K+1)*DEL + 2D0*FF(K)*DEL**2
     >  + 2D0*FF(K+2)*AA1 + FF(K+1)*XA(J,1))
     >  *(LOG(A2) - LOG(A1))

      ENDIF

      K = K + 3

 23   CONTINUE             
      
      G(J) = TEM
       TEM = 0
         K = 1


	ENDIF

 12   CONTINUE             

      ENDIF

      ENDIF
C
      IF (M .EQ. 3) THEN
        
      IF (T .GT. T1) THEN

      DO 93 J = NX-1, 1, -1

      TEM = TEM + FF(K)*(XA(J,-1) - XA(J+1,-1)) + (FF(K+1)/2.
     > + FF(K)*DEL)*(XA(J,-2) - XA(J+1,-2)) + (FF(K+2)/3.
     > + 2.*FF(K+1)*DEL/3.)*(XA(J,-3) - XA(J+1,-3))
     > + FF(K+2)*DEL*(XA(J,-4) - XA(J+1,-4))/2.
  
      G(J) = (XA(J,1)-DEL)*XA(J,2)*TEM  

         K = K + 3

 93   CONTINUE

      ELSE
  
      DO 13 I = NX-1, 1, -1

	IF (XA(I,1) .LE. DEL) THEN

	GOTO 666

	ELSE

           COMP = DEL/XA(I,1)

           IF(COMP .LT. 0.0005) THEN

         A1 = XA(I,1)*DIM - 1D0
         A2 = XA(I+1,1)*DIM -1D0

         TEM = TEM + DEL**3*(FF(K)*(XA(I,-1) - XA(I+1,-1)) 
     > + FF(K+1)/2.*(XA(I,-2) - XA(I+1,-2)) 
     > + FF(K+2)/3.*(XA(I,-3) - XA(I+1,-3)))

         ELSE

       A1 = XA(I,1)*DIM - 1D0
       A2 = XA(I+1,1)*DIM -1D0
      
      TEM = TEM + (FF(K+2)*(DEL*XA(I,-1) 
     > - DEL*XA(I+1,-1))
     > + (FF(K+1)*DEL + FF(K)*DEL**2 + FF(K+2))
     >  *(1D0/A1 - 1D0/A2)
     >  + (2D0*FF(K+2) + FF(K+1)*DEL)*(XA(I+1,0) - XA(I,0)
     >  + LOG(A1) - LOG(A2)))

      ENDIF
  
      G(I) = (A1+1D0)**2*A1*TEM 
       
         K = K + 3

	ENDIF

 13   CONTINUE
     
      ENDIF   

      ENDIF
C
      IF (M .EQ. 4) THEN

      IF (T .GT. T1) THEN

      DO 95 JI = NX-1, 1, -1      

      DO 96 IJ = NX-1, JI, -1

      TEM = TEM + FF(K)*(XA(IJ+1,2) - XA(IJ,2))/2. + (FF(K+1)
     > - 2.*FF(K)*XA(JI,1) + 2.*FF(K)*DEL)
     > *(XA(IJ+1,1) - XA(IJ,1)) + (FF(K+1)*XA(JI,2) - 2.*FF(K+2)
     > *(XA(JI,1) + DEL) - 4.*DEL*FF(K+1)*XA(JI,1) + 2.*FF(K)
     > *XA(JI,2)*DEL)*(XA(IJ,-1) - XA(IJ+1,-1)) + (FF(K+2)
     > *XA(JI,2)/2. - 2.*FF(K+2)*XA(JI,1)*DEL
     > + FF(K+1)*XA(JI,2)*DEL)*(XA(IJ,-2) - XA(IJ+1,-2)) 
     > + 2.*FF(K+2)*XA(JI,2)*DEL* (XA(IJ,-3) - XA(IJ+1,-3))/3. 
     > + (FF(K+2) - 2.*FF(K+1)*XA(JI,1) + FF(K)*XA(JI,2) 
     > + 2.*FF(K+1)*DEL
     > - 4.*FF(K)*XA(JI,1)*DEL)*(XA(IJ+1,0) - XA(IJ,0))

      K = K + 3

 96   CONTINUE             
      
      G(JI) = TEM
       TEM = 0
         K = 1

 95   CONTINUE
   
      ELSE

      DO 15 J = NX-1, 1, -1

	IF (XA(J,1) .LE. DEL) THEN

	GOTO 666

	ELSE

        AA1 = XA(J,1)*DIM      

      DO 26 I = NX-1, J, -1

         COMP = DEL/XA(J,1)
         COMP1 = DEL/XA(I,1)

         IF(COMP .LT. 0.0005 .OR. COMP1.LT.0.0005) THEN

             A1 = XA(I,1)*DIM - 1D0
             A2 = XA(I+1,1)*DIM -1D0

         TEM = TEM + FF(K)*(XA(I+1,2) - XA(I,2))/2. + (FF(K+1)
     > - 2.*FF(K)*XA(J,1) + 2.*FF(K)*DEL)
     > *(XA(I+1,1) - XA(I,1)) + (FF(K+1)*XA(J,2) - 2.*FF(K+2)
     > *(XA(J,1) + DEL) - 4.*DEL*FF(K+1)*XA(J,1) + 2.*FF(K)
     > *XA(J,2)*DEL)*(XA(I,-1) - XA(I+1,-1)) + (FF(K+2)
     > *XA(J,2)/2. - 2.*FF(K+2)*XA(J,1)*DEL
     > + FF(K+1)*XA(J,2)*DEL)*(XA(I,-2) - XA(I+1,-2)) 
     > + 2.*FF(K+2)*XA(J,2)*DEL* (XA(I,-3) - XA(I+1,-3))/3. 
     > + (FF(K+2) - 2.*FF(K+1)*XA(J,1) + FF(K)*XA(J,2) 
     > + 2.*FF(K+1)*DEL
     > - 4.*FF(K)*XA(J,1)*DEL)*(XA(I+1,0) - XA(I,0))


         ELSE

       A1 = XA(I,1)*DIM - 1D0
       A2 = XA(I+1,1)*DIM -1D0

      TEM = TEM + (FF(K+1) + 2D0*(DEL - XA(J,1))*FF(K))
     >  *(XA(I+1,1) - XA(I,1))
     >  + (FF(K+2) + FF(K+1)*DEL + FF(K)*DEL**2)
     >  *(1D0 - AA1)**2*(A1**-1
     >  - A2**-1) + FF(K+2)*AA1**2
     >  *(XA(I+1,0) - XA(I,0))
     >  + (AA1 - 1D0)*(FF(K)*XA(J,1)*DEL - 2D0*FF(K+1)*DEL
     >  - FF(K+2)*AA1 - 3D0*FF(K)*DEL**2
     >  - FF(K+2))*(LOG(A2) - LOG(A1)) + FF(K)*(XA(I+1,2)
     >  - XA(I,2))/2D0

      ENDIF

      K = K + 3

 26   CONTINUE             
      
      G(J) = TEM
       TEM = 0
         K = 1

	ENDIF

 15   CONTINUE             

      ENDIF

      ENDIF
C
      IF (M .EQ. 5) THEN

      IF (T .GT. T1) THEN

      DO 97 II = NX-1, 1, -1

      TEM = TEM + FF(K)*(XA(II+1,1) - XA(II,1)) + (FF(K+1)
     > + 2.*FF(K)*DEL)*(XA(II+1,0) - XA(II,0)) + (FF(K+2)
     > + 2.*FF(K+1)*DEL)*(XA(II,-1) - XA(II+1,-1)) + FF(K+2)*DEL
     > *(XA(II,-2) - XA(II+1,-2))
  
      G(II) = (XA(II,1)-DEL)*TEM
      K = K + 3

 97   CONTINUE

      ELSE

      DO 16 I = NX-1, 1, -1

	IF (XA(I,1) .LE. DEL) THEN

	GOTO 666

	ELSE

           COMP = DEL/XA(I,1)

           IF (COMP .LT. 0.0005) THEN

              A1 = XA(I,1)*DIM - 1D0
              A2 = XA(I+1,1)*DIM -1D0

       TEM = TEM + DEL*(FF(K)*(XA(I+1,1) - XA(I,1)) + (FF(K+1)
     > + 2.*FF(K)*DEL)*(XA(I+1,0) - XA(I,0)) + (FF(K+2)
     > + 2.*FF(K+1)*DEL)*(XA(I,-1) - XA(I+1,-1)) + FF(K+2)*DEL
     > *(XA(I,-2) - XA(I+1,-2)))

       ELSE

       A1 = XA(I,1)*DIM - 1D0
       A2 = XA(I+1,1)*DIM -1D0
      
      TEM = TEM + FF(K)*(XA(I+1,1) - XA(I,1))*DEL 
     >  + (FF(K+2) + FF(K+1)*DEL + FF(K)*DEL**2)
     >  *(A1**-1 - A2**-1) + (FF(K+1)*DEL 
     >  + FF(K)*2D0*DEL**2)*(LOG(A2) - LOG(A1))

      ENDIF
  
      G(I) = A1*TEM
      
      K = K + 3

	ENDIF

 16   CONTINUE

      ENDIF   

      ENDIF
C 
      IF (M .EQ. 6) THEN 
 
       IF (T .GT. T1) THEN 
 
       DO 98 JI = NX-1, 1, -1 
 
       DO 99 IJ = NX-1, JI, -1 
 
      TEM = TEM + FF(K)*(XA(IJ+1,2) - XA(IJ,2))/2. + (FF(K+1) 
     > - 2.*FF(K)*XA(JI,1) + FF(K)*DEL)*(XA(IJ+1,1) - XA(IJ,1)) 
     > + (FF(K+1)*XA(JI,2) - 2.*FF(K+2)*(XA(JI,1) - DEL/2.) 
     > - 2.*FF(K+1)*XA(JI,1)*DEL + FF(K)*XA(JI,2)*DEL) 
     > *(XA(IJ,-1) - XA(IJ+1,-1)) + (FF(K+2)*XA(JI,2)/2.  
     > - FF(K+2)*XA(JI,1)*DEL + FF(K+1)*XA(JI,2)*DEL/2.) 
     > *(XA(IJ,-2) - XA(IJ+1,-2)) + FF(K+2)*XA(JI,2)*DEL 
     > * (XA(IJ,-3) - XA(IJ+1,-3))/3. + (FF(K+2) - 2.*FF(K+1) 
     > *XA(JI,1) + FF(K)*XA(JI,2) + FF(K+1)*DEL - 2.*FF(K)*DEL 
     > *XA(JI,1))*(XA(IJ+1,0) - XA(IJ,0)) 
 
        K = K + 3 
 
 99   CONTINUE              
       
      G(JI) = TEM 
       TEM = 0 
         K = 1 
 
 98   CONTINUE 
 
       ELSE 
 
       DO 18 J = NX-1, 1, -1 
 
	IF (XA(J,1) .LE. DEL) THEN 
 
	GOTO 666 
 
	ELSE 
 
         AA1 = XA(J,1)*DIM 
 
       DO 29 I = NX-1, J, -1 
 
       A1 = XA(I,1)*DIM - 1D0 
       A2 = XA(I+1,1)*DIM -1D0 
       B1 = DEL*XA(I,-1) 
       B2 = DEL*XA(I+1,-1) 
 
      TEM = TEM + FF(K)*(XA(I+1,2) - XA(I,2))/2D0   
     >  + (FF(K+1) + FF(K)*DEL*(1D0 - 2D0*AA1)) 
     >  *(XA(I+1,1) - XA(I,1)) + FF(K+2)*AA1**2 
     >  *(B2 - B1) + (2D0*FF(K+2)  
     >  *AA1 - FF(K+1)*AA1*XA(J,1)  
     >  - FF(K+2)*AA1**2) 
     >  *(XA(I+1,0) - XA(I,0)) 
     >  + (AA1-1D0)**2*(FF(K+2) + FF(K+1)*DEL 
     >  + FF(K)*DEL**2)*(LOG(A2) - LOG(A1)) 
 
        K = K + 3 
 
 29   CONTINUE              
       
      G(J) = TEM 
       TEM = 0 
         K = 1 
 
	ENDIF 
 
 18   CONTINUE              
 
      ENDIF 
 
      ENDIF 
 
C 
      IF (M .EQ. 7) THEN 
 
      IF (T .GT. T1) THEN 
 
      DO 48 IJ = NX-1, 1, -1 
 
      TEM = TEM + FF(K)*(XA(IJ+1,1) - XA(IJ,1)) + (FF(K+1) + FF(K) 
     > *DEL)*(XA(IJ+1,0) - XA(IJ,0)) + (FF(K+2) + FF(K+1)*DEL) 
     > *(XA(IJ,-1) - XA(IJ+1,-1))  
 
       G(IJ) = (XA(IJ,1)-DEL)*TEM 
           K = K + 3              
       
 48   CONTINUE 
 
      ELSE 
       
      DO 222 I = NX-1, 1, -1 
 
	IF (XA(I,1) .LE. DEL) THEN 
  
	GOTO 666 
 
	ELSE 
 
       A1 = XA(I,1)*DIM - 1D0 
       A2 = XA(I+1,1)*DIM - 1D0 
 
      TEM = TEM + DEL*FF(K)*(XA(I+1,1) - XA(I,1))  
     >   + FF(K+2)*(XA(I,0) - XA(I+1,0))  
     >   + (FF(K+2) + FF(K+1)*DEL + FF(K)*DEL**2) 
     >   *(LOG(A2) - LOG(A1)) 
 
       G(I) = A1*TEM 
        
          K = K + 3              
 
	ENDIF 
       
 222  CONTINUE 
 
      ENDIF 
 
      ENDIF 
C 
      IF (M .EQ. 8) THEN 
 
C Integration of (y*f(y)-del*f(x))/(y*(y-del)) for x >> del, (y*f(y)-x*f(x))/(y*(y-x) 
C will be integrated in HINTEG and is directly embedded in SNRHS, NPRHS, NMRHS for 
C the singlet and non-singlet kernels.  
       
      IF (T .GT. T1) THEN      
       
      DO 83 LL = NX-1, 1, -1 
 
      DO 84 II = NX-1, LL, -1 
 
      TEM = TEM + FF(K)*(XA(II+1,2) - XA(II,2))/2D0 + (FF(K+1) 
     > + FF(K)*DEL)*(XA(II+1,1) - XA(II,1)) + (FF(K+2) + FF(K+1) 
     > *DEL)*(XA(II+1,0) - XA(II,0)) + DEL*XA(LL,1) 
     > *(XA(II+1,-1) - XA(II,-1))*(FF(K+1) + FF(K)*XA(LL,1)) 
 
         K = K + 3 
 
 84   CONTINUE     
 
      G(LL) = - TEM - F(LL)*LOG(1D0 - DEL) 
       TEM = 0 
         K = 1 
 
 
 83   CONTINUE 
 
      ELSE 
       
      DO 202 I = NX-1,1, -1 
 
C Determine whether Integration is in the BL or DGLAP region 
 
	IF (XA(I,1) .LE. DEL) THEN  
 
	GOTO 666 
 
	ELSE 
 
C Integration in the DGLAP region and only of (y*f(y)-del*f(x))/(y*(y-del))_+ the other 
C part of the regularization of (x-del)*f(y)/((y-x)*(y-del))_+ is taken care of 
C by HINTEGN in the kernels in the subroutines SNRHS, NSMRHS and NSPRHS. 
 
      DO 112 J = NX-1,I, -1 
 
       A1 = XA(J,1)*DIM - 1D0 
       A2 = XA(J+1,1)*DIM - 1D0 
 
       TEM = TEM + (FF(K+1)+FF(K)*DEL)*(XA(J+1,1)-XA(J,1)) 
     >   + FF(K)*(XA(J+1,2) - XA(J,2))/2D0  
     >   + (FF(K+1)*DEL - FF(KK+1)*XA(I,1) + FF(K)*DEL**2 
     >   - FF(KK)*XA(I,2) + FF(K+2) - FF(KK+2)) 
     >   *(LOG(A2)-LOG(A1)) + (FF(KK)*XA(I,2) + FF(KK+1)*XA(I,1) 
     >   + FF(KK+2))*(XA(J+1,0) - XA(J,0))  
 
         K = K + 3 
 
 112  CONTINUE 
  
   
       G(I) = - TEM - F(I)*LOG(1D0-DEL) 
        TEM = 0 
          K = 1 
         KK = KK + 3  
 
	ENDIF  
          
 202  CONTINUE 
 
      ENDIF 
 
      ENDIF 
 
C 
 
      IF (M .EQ. 9) THEN 
 
      IF (T .GT. T1) THEN 
 
      DO 995 JI = NX-1, 1, -1       
 
      DO 996 IJ = NX-1, JI, -1 
 
      TEM = TEM - XA(JI,1)*(2.*FF(K)*(DEL - XA(JI,1)) +  
     > FF(K+1))* 
     > (XA(IJ+1,0) - XA(IJ,0)) - XA(JI,1)*FF(K)*(XA(IJ+1,1) 
     > - XA(IJ,1)) + XA(JI,1)*(FF(K+2)*(XA(JI,1) - DEL)  
     > + 3./2.* 
     > DEL*FF(K+1)*XA(JI,1))*(XA(IJ,-2) - XA(IJ+1,-2))  
     > + XA(JI,1)* 
     > (-FF(K+2) + 3.*FF(K)*XA(JI,1)*DEL - 2.*FF(K+1)*(DEL  
     > - XA(JI,1))) 
     > *(XA(IJ,-1) - XA(IJ+1,-1)) + DEL*FF(K+2)*XA(JI,2)* 
     > (XA(IJ,-3) - XA(IJ+1,-3))  
 
      K = K + 3 
 
 996   CONTINUE              
       
      G(JI) = TEM 
       TEM = 0 
         K = 1 
 
 995   CONTINUE 
    
      ELSE 
 
      DO 155 J = NX-1, 1, -1 
 
	IF (XA(J,1) .LE. DEL) THEN 
  
	GOTO 666 
 
	ELSE 

        AA1 = XA(J,1)*DIM       

      DO 265 I = NX-1, J, -1 

         COMP = DEL/XA(J,1)
         COMP1 = DEL/XA(I,1)

         IF (COMP.LT.0.0005.OR.COMP1.LT.0.0005) THEN

      TEM = TEM - XA(J,1)*(2.*FF(K)*(DEL - XA(J,1)) +  
     > FF(K+1))* 
     > (XA(I+1,0) - XA(I,0)) - XA(J,1)*FF(K)*(XA(I+1,1) 
     > - XA(I,1)) + XA(J,1)*(FF(K+2)*(XA(J,1) - DEL)  
     > + 3./2.* 
     > DEL*FF(K+1)*XA(J,1))*(XA(I,-2) - XA(I+1,-2))  
     > + XA(J,1)* 
     > (-FF(K+2) + 3.*FF(K)*XA(J,1)*DEL - 2.*FF(K+1)*(DEL  
     > - XA(J,1))) 
     > *(XA(I,-1) - XA(I+1,-1)) + DEL*FF(K+2)*XA(J,2)* 
     > (XA(I,-3) - XA(I+1,-3))

      K = K + 3

      ELSE

       A1 = XA(I,1)*DIM - 1D0 
       A2 = XA(I+1,1)*DIM - 1D0 
 
      TEM = TEM + (XA(I,0) - XA(I+1,0) + LOG(A2) - LOG(A1)) 
     > *XA(J,1)*FF(K+1)*AA1 + XA(J,1)*(2.*FF(K)*DEL*(AA1  
     > - 1.) - FF(K+1))*(LOG(A2) - LOG(A1)) - XA(J,1)*FF(K+2) 
     > *AA1*(XA(I,-1) - XA(I+1,-1)) 
     > + AA1*FF(K)*(XA(I,2)*A1**-1 - XA(I+1,2)*A2**-1)  
     > + AA1*((AA1 - 1.)*(FF(K)*DEL**2 + FF(K+1)*DEL + FF(K+2))  
     > - DEL**2*FF(K))*(A1**-1 - A2**-1) 
      
 
      K = K + 3 

      ENDIF

 265   CONTINUE              
       
      G(J) = TEM 
       TEM = 0 
         K = 1 
 
	ENDIF 
 
 155   CONTINUE              
 
      ENDIF 
 
      ENDIF 
C     
       
      IF (M .EQ. 10) THEN 
 
      IF (T .GT. T1) THEN 
 
      DO 994 JK = NX-1, 1, -1       
 
      DO 997 IK = NX-1, JK, -1 
 
      TEM = TEM + (XA(IK+1,0) - XA(IK,0))*(XA(JK,1)*FF(K) 
     > *(2.*DEL -  
     > XA(JK,1)) + FF(K+1)*(2.*XA(JK,1) - DEL)) - XA(JK,1)* 
     > (FF(K+2)*(XA(JK,1) - 2.*DEL) + FF(K+1)*XA(JK,1)*DEL) 
     > *(XA(IK,-2) - XA(IK+1,-2))/2. + (FF(K+1)*XA(JK,1)* 
     > (2.*DEL - XA(JK,1)) + FF(K+2)*(2.*XA(JK,1) - DEL) 
     > - DEL*XA(JK,2)*FF(K))*(XA(IK,-1) - XA(IK+1,-1))  
     > + FF(K)*(2.*XA(JK,1) - DEL)*(XA(IK+1,1) - XA(IK,1))  
     > - 1./3.*DEL*XA(JK,2)*FF(K+2)*(XA(IK,-3) - XA(IK+1,-3)) 
 
      K = K + 3 
 
 997   CONTINUE              
       
      G(JK) = TEM 
       TEM = 0 
         K = 1 
 
 994   CONTINUE 
    
      ELSE 
 
      DO 154 JJ = NX-1, 1, -1 
 
	IF (XA(JJ,1) .LE. DEL) THEN 
  
	GOTO 666 
 
	ELSE 
 
        AA1 = XA(JJ,1)*DIM      
 
      DO 264 IJ = NX-1, JJ, -1 

         COMP = DEL/XA(JJ,1)
         COMP1 = DEL/XA(IJ,1)

         IF (COMP.LT.0.0005.OR.COMP1.LT.0.0005) THEN


       TEM = TEM + (XA(IJ+1,0) - XA(IJ,0))*(XA(JJ,1)*FF(K) 
     > *(2.*DEL -  
     > XA(JJ,1)) + FF(K+1)*(2.*XA(JJ,1) - DEL)) - XA(JJ,1)* 
     > (FF(K+2)*(XA(JJ,1) - 2.*DEL) + FF(K+1)*XA(JJ,1)*DEL) 
     > *(XA(IJ,-2) - XA(IJ+1,-2))/2. + (FF(K+1)*XA(JJ,1)* 
     > (2.*DEL - XA(JJ,1)) + FF(K+2)*(2.*XA(JJ,1) - DEL) 
     > - DEL*XA(JJ,2)*FF(K))*(XA(IJ,-1) - XA(IJ+1,-1))  
     > + FF(K)*(2.*XA(JJ,1) - DEL)*(XA(IJ+1,1) - XA(IJ,1))  
     > - 1./3.*DEL*XA(JJ,2)*FF(K+2)*(XA(IJ,-3) - XA(IJ+1,-3))      

       K = K + 3

       ELSE

       A1 = XA(IJ,1)*DIM - 1D0 
       A2 = XA(IJ+1,1)*DIM - 1D0 
 
      TEM = TEM - (1. - AA1)**2*(FF(K)*DEL**2 + FF(K+1) 
     > *DEL + FF(K+2))* 
     > (LOG(A2) - LOG(A1) + XA(IJ,0) - XA(IJ+1,0)) +  
     > (XA(IJ+1,0) - XA(IJ,0))*(FF(K+1)*DEL*(2.* 
     > AA1 - 1.) - FF(K)*DEL**2*(1. - AA1)**2) + FF(K)*DEL* 
     > (2.*AA1 - 1.)*(XA(IJ+1,1) - XA(IJ,1)) + XA(JJ,1) 
     > *AA1*FF(K+2)*(XA(IJ,-1) - XA(IJ+1,-1))   
      
 
      K = K + 3 

      ENDIF

 264   CONTINUE              
       
      G(JJ) = TEM 
       TEM = 0 
         K = 1 
 
	ENDIF 
 
 154   CONTINUE              
 
      ENDIF 
 
      ENDIF 
       
C     
       
      IF (M .EQ. 11) THEN 
 
      IF (T .GT. T1) THEN 
 
      DO 993 JL = NX-1, 1, -1       
 
      DO 998 IL = NX-1, JL, -1 
 
      TEM = TEM + (XA(IL+1,0) - XA(IL,0))*(FF(K+1)*(3.* 
     > XA(JL,1) - DEL) - 2.*XA(JL,1)*FF(K)*(2.*XA(JL,1) 
     > -3.*DEL)) - (XA(IL,-1) - XA(IL+1,-1))*(4.* 
     > XA(JL,1)*FF(K+1)*(XA(JL,1) - 3./2.*DEL) - FF(K+2) 
     > *(3.*XA(JL,1) - DEL) + 6.*XA(JL,2)*DEL*FF(K)) +  
     > (XA(IL+1,1) - XA(IL,1))*FF(K)*(3.*XA(JL,1) - DEL) 
     > - (XA(IL,-2) - XA(IL+1,-2))*XA(JL,1)*(3.*XA(JL,1) 
     > *FF(K+1)*DEL + FF(K+2)*(2.*XA(JL,1) - 3.*DEL)) 
     > -2.*(XA(IL,-3) - XA(IL+1,-3))*XA(JL,2)*DEL*FF(K+2) 
      
 
      K = K + 3 
 
 998   CONTINUE              
       
      G(JL) = TEM 
       TEM = 0 
         K = 1 
 
 993   CONTINUE 
    
      ELSE 
 
      DO 153 J = NX-1, 1, -1 
 
	IF (XA(J,1) .LE. DEL) THEN 
  
	GOTO 666 
 
	ELSE 
 
        AA1 = XA(J,1)*DIM       
 
      DO 263 I = NX-1, J, -1 

         COMP = DEL/XA(J,1)
         COMP1 = DEL/XA(I,1)

         IF (COMP.LT.0.0005.OR.COMP1.LT.0.0005) THEN

       TEM = TEM + (XA(I+1,0) - XA(I,0))*(FF(K+1)*(3.* 
     > XA(J,1) - DEL) - 2.*XA(J,1)*FF(K)*(2.*XA(J,1) 
     > -3.*DEL)) - (XA(I,-1) - XA(I+1,-1))*(4.* 
     > XA(J,1)*FF(K+1)*(XA(J,1) - 3./2.*DEL) - FF(K+2) 
     > *(3.*XA(J,1) - DEL) + 6.*XA(J,2)*DEL*FF(K)) +  
     > (XA(I+1,1) - XA(I,1))*FF(K)*(3.*XA(J,1) - DEL) 
     > - (XA(I,-2) - XA(I+1,-2))*XA(J,1)*(3.*XA(J,1) 
     > *FF(K+1)*DEL + FF(K+2)*(2.*XA(J,1) - 3.*DEL)) 
     > -2.*(XA(I,-3) - XA(I+1,-3))*XA(J,2)*DEL*FF(K+2) 
     
      K = K + 3 

      ELSE

       A1 = XA(I,1)*DIM - 1D0 
       A2 = XA(I+1,1)*DIM - 1D0 
 
      TEM = TEM + 2.*(LOG(A2) - LOG(A1))*(1. - AA1) 
     > *((FF(K)*DEL**2 + FF(K+1)*DEL)*(AA1 - 1.) + 
     > FF(K)*XA(J,1)*DEL + 1./2.*FF(K+1)*DEL) 
     > + 2.*(XA(I+1,0) - XA(I,0))*AA1*XA(J,1)*FF(K+1) 
     > + (XA(I+1,1) - XA(I,1))*FF(K)*DEL*(3.*AA1 
     > - 1.) + 2.*XA(J,1)*AA1*FF(K+2)*(XA(I,-1) -  
     > XA(I+1,-1)) + (FF(K)*DEL**2 + FF(K+1)*DEL +  
     > FF(K+2))*(1. - AA1)*(1. - 2.*AA1)*(A2**-1 -  
     > A1**-1) 
      
 
      K = K + 3 

      ENDIF

 263   CONTINUE              
       
      G(J) = TEM 
       TEM = 0 
         K = 1 
 
	ENDIF 
 
 153   CONTINUE              
 
      ENDIF 
 
      ENDIF 
C 
 
	IF (M .EQ. 12) THEN 
 
          K = 3*NX-5
	 
          G(1) = 0.0

          G(ITC) = 0.0

          TEM = 0.0

	DO 371 IJ = 2, NX-1
 
	   IF (XA(IJ,1) .GE. DEL) THEN  
 
	     GOTO 666 
 
	ELSE 
  
	A1 = 1D0 - XA(IJ-1,1)*DIM 
        A2 = 1D0 - XA(IJ,1)*DIM 
             
      TEM = TEM + FF(K)*(XA(IJ-1,2) - XA(IJ,2))/2D0  
     >  + (FF(K)*DEL+FF(K+1))*(XA(IJ-1,1)-XA(IJ,1))  
     >  + (FF(K+2) + FF(K+1)*DEL + FF(K)*DEL**2) 
     >  *(LOG(A1) - LOG(A2)) 

      
	G(IJ) = TEM 

	K = K - 3 

	ENDIF 
 
 371	CONTINUE 

	ENDIF 

C 
 
      IF (M .EQ. 13) THEN 
 
          K = 3*NX-5 
  
	G(1) = 0.0

        G(ITC) =  0.0
        
        TEM = 0.0

	DO 171 IJ = 2, NX-1 
 
	IF (XA(IJ,1) .GE. DEL) THEN  
 
	  GOTO 666 
 
	ELSE 
 
	A1 = 1D0 - XA(IJ-1,1)*DIM 
        A2 = 1D0 - XA(IJ,1)*DIM 
      
      IF (XA(IJ,1).GT.0.001*DEL) THEN

      TEM = TEM + FF(K)*(XA(IJ,1) - XA(IJ-1,1))*DEL  
     >  + (FF(K+2) + FF(K+1)*DEL + FF(K)*DEL**2) 
     >  *(1D0/A2 - 1D0/A1) + (FF(K+1)*DEL  
     >  + FF(K)*2D0*DEL**2)*(LOG(A2) - LOG(A1)) 

        ELSEIF(XA(IJ,1).LT.0.001*DEL) THEN

      TEM = TEM + FF(K)*(XA(IJ,1) - XA(IJ-1,1))*DEL  
     >  + (FF(K+2) + FF(K+1)*DEL + FF(K)*DEL**2) 
     >  *((1D0-A2) + (1D0-A2)**2 - (1D0-A1) - (1D0-A1)**2) 
     >  - (FF(K+1)*DEL + FF(K)*2D0*DEL**2)*
     > ((1D0-A2) + (1D0-A2)**2/2D0 - (1D0-A1) - (1D0-A1)**2/2D0)

      ENDIF
           
	G(IJ) = TEM

	K = K - 3 
 
	ENDIF 
 
 171	CONTINUE 

	ENDIF 
 
C  
	IF (M .EQ. 14) THEN 
 
          K = 3*NX-5 
 
	G(1) = 0.0 

        G(ITC) = 0.0

        TEM1 = 0.0

	DO 571 IJ = 2, NX-1
 
	IF (XA(IJ,1) .GE. DEL) THEN  
 
	   GOTO 666 
 
	ELSE 
       
       TEM1 = TEM1 + FF(K)*(XA(IJ,3) - XA(IJ-1,3))/3D0  
     >  + FF(K+1)*(XA(IJ,2) - XA(IJ-1,2))/2D0  
     >  + FF(K+2)*(XA(IJ,1) - XA(IJ-1,1)) 
 
        G(IJ) = TEM1 

	K = K - 3 
 
	ENDIF 
 
 571	CONTINUE 

	ENDIF 
 
C 
	IF (M .EQ. 15) THEN

           G(1) = -F(1)*LOG(DEL)

	DO 203 I = 2, ITC-1

        AB2 = XA(I,1)*DIM 
	AB1 = 1D0 - AB2

           TEM1 = 0D0

           K = 3*NX-5 

C Integration in the BL region		

	    DO 113 JK = 2, I

	     AA1 = 1D0 - XA(JK-1,1)*DIM	 
       	     AA2 = 1D0 - XA(JK,1)*DIM 
	     BB1 = 1D0 - XA(JK-1,1)*XA(I,-1) 
	     BB2 = 1D0 - XA(JK,1)*XA(I,-1) 
             
C Endpoint singularity in integrand cancels, explicitly! 

	  IF (JK .EQ. I) THEN 

             IF (XA(JK,1).GT.0.01*DEL) THEN

	     TEM1 = TEM1 + FF(K)*DEL*AB1
     >   *(XA(JK,1) - XA(JK-1,1)) 
     >   + (FF(K+2) + FF(K+1)*DEL + FF(K)*DEL**2
     >   - F(I))*(LOG(AA2)-LOG(AA1))

             ELSE

        TEM1 = TEM1 + FF(K)*DEL*AB1
     >   *(XA(JK,1) - XA(JK-1,1)) 
     >   - (FF(K+2) + FF(K+1)*DEL + FF(K)*DEL**2
     >   - F(I))*((1D0-AA2) + (1D0-AA2)**2/2D0 - (1D0-AA1) 
     >   - (1D0-AA1)**2/2D0)  

        ENDIF

	ELSE 

C Value of the integration away from the endpoint singularity in the integrand.	 


           IF (XA(JK,1).GT.0.01*DEL) THEN

	     TEM1 = TEM1 + FF(K)*DEL*AB1*(XA(JK,1) - XA(JK-1,1))
     >   + (FF(K+1)*DEL + FF(K)*DEL**2 + FF(K+2) - F(I))
     >   *(LOG(AA2)-LOG(AA1)) + (F(I) - FF(K)*XA(I,2) 
     >   - FF(K+1)*XA(I,1) - FF(K+2))*(LOG(BB2) - LOG(BB1)) 

             ELSE

            TEM1 = TEM1 + FF(K)*DEL*AB1
     >   *(XA(JK,1) - XA(JK-1,1))
     >   - (FF(K+1)*DEL + FF(K)*DEL**2
     >   + FF(K+2) - F(I))*((1D0-AA2) + (1D0-AA2)**2/2D0 
     >   - (1D0-AA1) - (1D0-AA1)**2/2D0) + (F(I) 
     >   - FF(K)*XA(I,2) - FF(K+1)*XA(I,1) - FF(K+2)) 
     >   *(LOG(BB2) - LOG(BB1)) 

            ENDIF

	K = K - 3

	ENDIF

 113  CONTINUE

	G(I) = TEM1 + F(I)*LOG(AB2)       

 203    CONTINUE 

	ENDIF 

	IF (M .EQ. 16) THEN

           TEM = 0.0

	DO 670 I = NX-1, 1, -1
      
        IF (I.EQ.1) THEN

        TEM = TEM + 1D0/2D0*FF(K)*(XA(I+1,2) - XA(I,2))
     >  + FF(K+1)*(XA(I+1,1) - XA(I,1))
     >  + FF(K+2)*(XA(I+1,0))

        ELSE

         TEM = TEM + 1D0/2D0*FF(K)*(XA(I+1,2) - XA(I,2))
     >  + FF(K+1)*(XA(I+1,1) - XA(I,1))
     >  + FF(K+2)*(XA(I+1,0) - XA(I,0))

        ENDIF
 
        G(I) = TEM

	K = K + 3 

 670	CONTINUE

        ENDIF

C
        IF (M .EQ. 17) THEN

           TEM = 0.0
           G(1) = 0.0
           
	DO 680 I = NX-1, 1, -1 
      
        IF (I.EQ.1) THEN

        TEM = TEM + FF(K)*(XA(I+1,1) - XA(I,1))
     >  - FF(K+2)*(XA(I+1,-1))
     >  + FF(K+1)*(XA(I+1,0))

        ELSE

        TEM = TEM + FF(K)*(XA(I+1,1) - XA(I,1))
     >  - FF(K+2)*(XA(I+1,-1) - XA(I,-1))
     >  + FF(K+1)*(XA(I+1,0) - XA(I,0))

        ENDIF
 
        G(I) = TEM

	K = K + 3 

 680	CONTINUE

        ENDIF

C

        IF (M .EQ. 19) THEN

           TEM = 0.0

        DO 600 IP = 1,ITC

           AB1 = 1D0 - XA(IP,1)*DIM

	DO 610 I = NX-1, ITC, -1

           AK1 = XA(I,1)*DIM
           AK2 = XA(I+1,1)*DIM
           AA1 = AK1/AB1 - 1D0
           AA2 = AK2/AB1 - 1D0

           COMP = DEL*XA(I,-1)

           IF (COMP.GT.0.001.AND.COMP.LT.1D0) THEN
      
        TEM = TEM + (LOG(AA2) - LOG(AA1))*(FF(K)*DEL**2*AB1**2
     > + FF(K+1)*DEL*AB1 + FF(K+2)) - (XA(I+1,0) - XA(I,0))*
     > FF(K+2) + FF(K)*DEL*AB1*(XA(I+1,1) - XA(I,1))

        ELSEIF(COMP.EQ.1D0.AND.IP.GT.1) THEN

        TEM = TEM + (LOG(AK2-AB1) - LOG(1D0-AB1))*(FF(K)*DEL**2*AB1**2
     > + FF(K+1)*DEL*AB1 + FF(K+2)) - (XA(I+1,0) - XA(I,0))*
     > FF(K+2) + FF(K)*DEL*AB1*(XA(I+1,1) - XA(I,1))


        ELSEIF(COMP.EQ.1D0.AND.IP.EQ.1) THEN

        TEM = TEM + (LOG(AK2-AB1)+LOG(DEL))*(FF(K)*DEL**2*AB1**2
     > + FF(K+1)*DEL*AB1 + FF(K+2)) - (XA(I+1,0) - XA(I,0))*
     > FF(K+2) + FF(K)*DEL*AB1*(XA(I+1,1) - XA(I,1))

        ELSE

        TEM = TEM + (XA(I+1,0) - AB1/AK2 -(AB1/AK2)**2/2D0 
     > - XA(I,0) + AB1/AK1 + (AB1/AK1)**2/2D0)
     > *(FF(K)*DEL**2*AB1**2
     > + FF(K+1)*DEL*AB1 + FF(K+2)) - (XA(I+1,0) - XA(I,0))*
     > FF(K+2) + FF(K)*DEL*AB1*(XA(I+1,1) - XA(I,1))   

        ENDIF

	K = K + 3

 610	CONTINUE

        G(IP) = TEM

        TEM = 0
        K = 1
      
 600    CONTINUE

        ENDIF
C

        IF (M .EQ. 20) THEN

           TEM = 0.0

        DO 601 IP = 1,ITC

           AB1 = 1D0 - XA(IP,1)*DIM

	DO 611 I = NX-1, ITC, -1


           AK1 = XA(I,1)*DIM
           AK2 = XA(I+1,1)*DIM
           AA1 = AK1/AB1 - 1D0
           AA2 = AK2/AB1 - 1D0

           COMP = DEL*XA(I,-1)

           IF (COMP.GT.0.01.AND.COMP.LT.1D0) THEN
      
        TEM = TEM + (LOG(AA2) - LOG(AA1))*(FF(K)*DEL**2*AB1**2
     > + FF(K+1)*DEL*AB1 + FF(K+2)) + 1D0/2D0*FF(K)
     > *(XA(I+1,2)-XA(I,2)
     > + 2D0*DEL*AB1*(XA(I+1,1)-XA(I,1)))+FF(K+1)
     > *(XA(I+1,1) - XA(I,1))

        ELSEIF (COMP.EQ.1D0.AND.IP.GT.1) THEN

        TEM = TEM + (LOG(AK2-AB1)-LOG(1D0-AB1))*(FF(K)*DEL**2*AB1**2
     > + FF(K+1)*DEL*AB1 + FF(K+2)) + 1D0/2D0*FF(K)
     > *(XA(I+1,2)-XA(I,2)
     > + 2D0*DEL*AB1*(XA(I+1,1)-XA(I,1)))+FF(K+1)
     > *(XA(I+1,1) - XA(I,1))

        ELSEIF (COMP.EQ.1D0.AND.IP.EQ.1) THEN

        TEM = TEM + (LOG(AK2-1D0)+LOG(DEL))*(FF(K)*DEL**2
     > + FF(K+1)*DEL + FF(K+2)) + 1D0/2D0*FF(K)
     > *(XA(I+1,2)-XA(I,2)
     > + 2D0*DEL*(XA(I+1,1)-XA(I,1)))+FF(K+1)
     > *(XA(I+1,1) - XA(I,1))

        ELSE

       TEM = TEM + (XA(I+1,0) - AB1/AK2 -(AB1/AK2)**2/2D0 
     > - XA(I,0) + AB1/AK1 + (AB1/AK1)**2/2D0)*
     > (FF(K)*DEL**2*AB1**2 + FF(K+1)*DEL*AB1 + FF(K+2)) 
     > + 1D0/2D0*FF(K)*(XA(I+1,2)-XA(I,2)
     > + 2D0*DEL*AB1*(XA(I+1,1)-XA(I,1)))+FF(K+1)*(XA(I+1,1) 
     > - XA(I,1))   

       ENDIF

	K = K + 3

 611    CONTINUE
       
        G(IP) = TEM

        TEM = 0
        K = 1
      
 601    CONTINUE

        ENDIF

C
	IF (M .EQ. 21) THEN

           TEM = 0.0

	DO 770 I = NX-1, 1, -1
      
        TEM = TEM + 1D0/3D0*FF(K)*(XA(I+1,3) - XA(I,3))
     >  + FF(K+1)*(XA(I+1,2) - XA(I,2))/2D0
     >  + FF(K+2)*(XA(I+1,1) - XA(I,1))
 
        G(I) = TEM

	K = K + 3 

 770	CONTINUE

        ENDIF

C

        IF (M .EQ. 22) THEN

           TEM = 0.0

	DO 570 I = NX-1, 1, -1

           IF (XA(I,1).GT.0.1*DEL) THEN
      
        TEM = TEM + 1D0/3D0*FF(K)*(XA(I+1,3) - XA(I,3))
     >  + FF(K+1)*(XA(I+1,2) - XA(I,2))/2D0
     >  + FF(K+2)*(XA(I+1,1) - XA(I,1))

           ELSE

              A1 = XA(I,1)*DIM
              A2 = XA(I+1,1)*DIM

          TEM = TEM + 1D0/3D0*FF(K)*(A2**3 - A1**3)*DEL**3
     >  + FF(K+1)*(XA(I+1,2) - XA(I,2))/2D0
     >  + FF(K+2)*(XA(I+1,1) - XA(I,1))    
              
          ENDIF

        G(I) = TEM

	K = K + 3 

 570	CONTINUE

        ENDIF


 666    RETURN 


C              *********************************************         
       
      END         
    
      SUBROUTINE HINTEGN (I,NX,IV,D,FF,F,G)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      CHARACTER MSG*80 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (MXPN = MXF * 2 + 2) 
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN) 
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)           
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX) 
      COMMON / VARBAB / GB(NDG, NDH, MXX), H(NDH, MXX, M1:M2) 
C 
      DIMENSION   FF(3*MXX), G(NX), F(NX)

      DIM = 1./D

      IF (I.EQ.1) THEN
         
         TEM = 0.0

         K = 1

         DO 999 M = IV,1,-1

         DO 998 J = NX-1,M,-1

C     Cancellation of endpoint singularity!

         IF (J.EQ.M) THEN

          IF (J.EQ.1) THEN

          TEM = 1./2.*FF(K)*(XA(J+1,2)-XA(J,2))
     > + (XA(M,1)*FF(K) + FF(K+1))*(XA(J+1,1)-XA(J,1))
     > + F(M)*(XA(J+1,0)) + TEM

          ELSE

         TEM = 1./2.*FF(K)*(XA(J+1,2)-XA(J,2))
     > + (XA(M,1)*FF(K) + FF(K+1))*(XA(J+1,1)-XA(J,1))
     > + F(M)*(XA(J+1,0) - XA(J,0)) + TEM

         ENDIF

         K = K + 3

         ELSE

C     Normal integration

         AB2 = 1. - XA(M,1)*XA(J+1,-1)
         AB1 = 1. - XA(M,1)*XA(J,-1)

         TEM = TEM + 1./2.*FF(K)*(XA(J+1,2)-XA(J,2))
     > + (XA(M,1)*FF(K) + FF(K+1))*(XA(J+1,1)-XA(J,1))
     > + F(M)*(XA(J+1,0) - XA(J,0)) + (FF(K)*XA(M,2)
     > + FF(K+1)*XA(M,1) + FF(K+2) - F(M))*(XA(J+1,0) 
     > - XA(J,0) + LOG(AB2) - LOG(AB1))

         K = K + 3

         ENDIF

 998     CONTINUE

c         G(M) = TEM + F(M)*(LOG(1.-XA(M,1))-LOG(1.-D))
         G(M) = TEM + F(M)*(LOG(1.-XA(M,1)))

         TEM = 0.0

         K = 1

 999     CONTINUE

         ELSE
         
         G(NX) = 0.0
         TEM = 0.0
         K = 1 
        
         DO 997 MN = NX-1,IV+1,-1

         DO 996 J1 = NX-1,MN,-1

            IF (J1.EQ.MN) THEN

         TEM = 1./2.*FF(K)*(XA(J1+1,2)-XA(J1,2))
     > + (XA(MN,1)*FF(K) + FF(K+1))*(XA(J1+1,1)-XA(J1,1))
     > + F(MN)*(XA(J1+1,0) - XA(J1,0))+TEM

         K = K + 3

         ELSE

C     Normal integration

         AB2 = 1. - XA(MN,1)*XA(J1+1,-1)
         AB1 = 1. - XA(MN,1)*XA(J1,-1)

         TEM = TEM + 1./2.*FF(K)*(XA(J1+1,2)-XA(J1,2))
     > + (XA(MN,1)*FF(K) + FF(K+1))*(XA(J1+1,1)-XA(J1,1))
     > + F(MN)*(XA(J1+1,0) - XA(J1,0)) + (FF(K)*XA(MN,2)
     > + FF(K+1)*XA(MN,1) + FF(K+2) - F(MN))*(XA(J1+1,0) 
     > - XA(J1,0) + LOG(AB2) - LOG(AB1))

         K = K + 3

         ENDIF

 996     CONTINUE 

         G(MN) = TEM + F(MN)*LOG(1.-XA(MN,1))

         TEM = 0.0

         K = 1

 997     CONTINUE

      ENDIF

      RETURN
C
      END
C

      SUBROUTINE HINTEG (NX, F, H) 
C 
C       Computes the integral [yF(y)-xF(x)]/(y-x) * dy/y over the range [x, 1]; 
C       then add F(x) * Ln (1-x) to get Int 1/(1-x/y)(sub+)F(y)dy/y. 
C 
C       The input function F must be specified on an array of size NX over 
C       the range (0, 1] of the x variable.  In order to allow a possible 
c       singularity at x = 0, the first mesh-point is at x = 1/nx, not 0. 
C       The output function H is given on the same array as above. 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (MXPN = MXF * 2 + 2) 
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN) 
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2) 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
      COMMON / HINTEC / GH(NDG, MXX) 
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX) 
C 
      DIMENSION F(NX), H(NX), G(MXX) 
C 
      DZ = 1D0 / (NX-1) 
C                                       Loop to calculate the required integral 
C                                               |--|--|-----------|--|--| 
C                                         Iy:   I I+1 ...              nx 
C                                         Kz:   1  2   ...              np 
C 
C                                                ------------------------------ 
C                                                     Do each integral in turn. 
      DO 20 I = 1, NX-2 
C                                         Number of points in the I-th integral 
         NP = NX - I + 1 
C                                  Evaluate the first two bins of the integrals 
 
         TEM = GH(1,I)*F(I) + GH(2,I)*F(I+1) + GH(3,I)*F(I+2) 
 
C                                  Evaluate the integrand for the I-th integral 
         DO 30 KZ = 3, NP 
            IY = I + KZ - 1 
C                       DXDZ is the Jacobian due to the change of variable X->Z 
C 
            W = XA(I,1) / XA(IY,1) 
            G(KZ) = DXTZ(IY)*(F(IY)-W*F(I))/(1.-W) 
C 
   30    CONTINUE 
C 
         HTEM = SMPSNA (NP-2, DZ, G(3), ERR) 
         TEM1 = F(I) * ELY(I) 
         H(I) = TEM + HTEM + TEM1 

   20 CONTINUE 
C 
      H(NX-1) = F(NX) - F(NX-1) + F(NX-1) * (ELY(NX-1) - XA(NX-1,0)) 
      H(NX)   = 0 
C 
      RETURN 
C                        **************************** 
      END 

      Function xfrmzn(i,n,x)
C
C     function to calculate the x-grid in the ERBL region using 
C     a hyperbolic tangent thus ensuring a high sample rate around 
C     xmin and y=del! This routine also generates the first 10 points 
C     in the DGLAP region so as to ensure an accurate integration result
C     near y=del also in the DGLAP region. Arctanh(x) = 1/2*ln((1+x)/(1-x)) 
C     is defined through the logs!
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      PARAMETER (MXX = 1050) 
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX, DEL,IV 



      dim = 1./DEL

      if (i.eq.1) then

      MS = MOD(NX,4)

      IF (MS.EQ.0) THEN

      NC = NX/8

      ELSE

      NC = (NX/4+1)/2

      ENDIF
      
      IF (DEL.LT.0.1) THEN
         MS = MOD(NX,4)

      IF (MS.EQ.0) THEN

      NC = NX/8

      ELSE

      NC = (NX/4+1)/2

      ENDIF
      a = 10.
      s = 4.*((1.*(NX-1))/(1.*(NX-4)))
      ELSE
         MS = MOD(NX,2)

      IF (MS.EQ.0) THEN

      NC = NX/4

      ELSE

      NC = (NX/2+1)/2

      ENDIF
      a = 7.
      s = 2.*((1.*(NX-1))/(1.*(NX-2)))
      ENDIF

      xp = x*s

      a1 = a*(xp-1./2.)**2
      b = dexp(-a/4.)
      b1 = dexp(-a1)
    
      if (n.LE.NC) then

      tem2 = DEL/2.*(b1 + (4.*xp*(1-xp)-1.)*b)

      ELSE

      tem2 = - DEL/2.*(b1 + (4.*xp*(1-xp)-1.)*b) + DEL
 
      ENDIF

      if (del.lt.0.1.and.n.eq.NX/4) then
         tem2 = DEL
      elseif(del.ge.0.1.and.n.eq.NX/2) then
         tem2 = DEL
      endif

      xfrmzn = tem2

      else

          b = - LOG(dim)

      if (del.lt.0.1) then

      a = 1./4.*(NX-4)/(1.*(NX-1))
      a1 = LOG(a)

      else

      a = 1./2.*(NX-2)/(1.*(NX-1))   
      a1 = LOG(a)

      endif

      c = b/a1

      tem1 = x**c

      xfrmzn = tem1

      endif

      return
C
      entry dxdzn(i,n,x)

      dim = 1./DEL

      if (i.eq.1) then

      IF (DEL.LT.0.1) THEN
      MS = MOD(NX,4)

      IF (MS.EQ.0) THEN

      NC = NX/8

      ELSE

      NC = (NX/4+1)/2

      ENDIF
      a = 10.
      s = 4.*((1.*(NX-1))/(1.*(NX-4)))
      ELSE
      MS = MOD(NX,2)

      IF (MS.EQ.0) THEN

      NC = NX/4

      ELSE

      NC = (NX/2+1)/2

      ENDIF
      a = 7.
      s = 2.*((1.*(NX-1))/(1.*(NX-2)))
      ENDIF

      xp = x*s

      a1 = a*(xp-1./2.)**2
      b = dexp(-a/4.)
      b1 = dexp(-a1)

      if (n.LE.NC) then

      tem3 = DEL*(1.-2.*xp)*s*(a*b1 + 4.*b)/2.

      ELSE

      tem3 = - DEL*(1.-2.*xp)*s*(a*b1 + 4.*b)/2.

      ENDIF

      if (abs(tem3).lt.1E-10) tem3 = 0.0

      dxdzn = tem3

      else

      b = - LOG(dim)

      if (del.lt.0.1) then

      a = 1./4.*(NX-4)/(1.*(NX-1))
      a1 = LOG(a)

      else

      a = 1./2.*(NX-2)/(1.*(NX-1))   
      a1 = LOG(a)

      endif

      c = b/a1
      
      dxdzn = c*x**(c-1.)

      endif

      return     
c
      end

C

       SUBROUTINE PARPDF (IACT, NAME, VALUE, IRET) 
 

C ==================================================================== 
C      SUBROUTINE PARPDF (IACT, NAME, VALUE, IRET) 
C 
C               Actions: 0      type list of variables on unit VALUE. 
C                        1      set variable with name NAME to VALUE, if 
C                               it exists, Else set IRET to 0. 
C                        2      find value of variable.  If it does not exist, 
C                               set IRET to 0. 
C                        3      request values of all parameters from terminal. 
C                        4      type list of all values on unit VALUE. 
c 
C               IRET =   0      variable not found. 
C                        1      successful search. 
C                        2      variable found, but bad value. 
C                        3      bad value for IACT. 
C                        4      no variable search (i.e., IACT is 0, 3, or 4). 
C 
C               If necessary, VALUE is converted to integer by NINT 
C 
C              Use ILEVEL and START1 to start search for variable names close 
C              to previous name to ensure effiency when reading in many values. 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C 
      CHARACTER NAME*(*), Uname*10 
C 
      LOGICAL START1 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
C 
      DATA ILEVEL, LRET / 1, 1 / 
 
      JRET = IRET 
      CALL UPC (NAME, Ln, Uname) 
      IF (IACT .EQ. 0 .OR. IACT .EQ. 4) 
     >    IVALUE = NINT (VALUE) 
      START1 = (IACT .NE. 1) .AND. (IACT .NE. 2) 
      IF (START1)  ILEVEL = 1 
C 
      GOTO (1, 2), ILEVEL 
C 
    1 START1 = .TRUE. 
      ILEVEL = 0 
      CALL PARQCD (IACT, Uname(1:Ln), VALUE, JRET) 
              IF (JRET .EQ. 1)  GOTO 11 
              IF (JRET .EQ. 2)  GOTO 12 
              IF (JRET .EQ. 3)  GOTO 13 
              IF (JRET .GT. 4)  GOTO 15 
              ILEVEL =  ILEVEL + 1 
    2 CALL EVLPAR (IACT, Uname(1:Ln), VALUE, JRET) 
              IF (JRET .EQ. 1)  GOTO 11 
              IF (JRET .EQ. 2)  GOTO 12 
              IF (JRET .EQ. 3)  GOTO 13 
              IF (JRET .GT. 4)  GOTO 15 
              ILEVEL =  ILEVEL + 1 
C 
      IF (.NOT. START1) GOTO 1 
C 
      IF (JRET .EQ. 0)  GOTO 10 
C 
C                       Arrive here if IACT = 0, 3 and all is OK (IRET=4): 
    9 CONTINUE 
C                       WRITE (IVALUE, 100) 
      GOTO 14 
C                       Exits: 
   10 CONTINUE 
   11 CONTINUE 
   12 CONTINUE 
   13 CONTINUE 
   14 CONTINUE 
   15 CONTINUE 
C                                   LRET is used for debugging purpose only 
      IF (JRET .NE. 4) LRET = JRET 
      IF (LRET.EQ.0 .OR. LRET.EQ.2 .OR. LRET.EQ.3) THEN 
        PRINT *, 'Error in PARPDF: IRET, IACT, NAME, VALUE =', 
     >  LRET, IACT, NAME, VALUE 
        PAUSE 
      EndIf 
      IRET= JRET 
      RETURN 
C 
  100 FORMAT (/) 
C                        **************************** 
      END 

      SUBROUTINE EVLPAR (IACT, NAME, VALUE, IRET) 
C 
C               For STANDARD codes on Iact and Iret, see SUBROUTINE PARPDF 
C 
C               Additional options added for Version 7.2 and up of PDF package: 
C 
C          IACT =        5      find value of variable from the alternate set. 
C                        6      type list of all values from the alternate set. 
C 
c               These options must be called from EVLPAR directly, not from the 
C               front-end unit PARPDF because PARQCD cannot handle Iact > 4 
C 
C               NAME is assumed upper-case. 
C               If necessary, VALUE is converted to integer by NINT 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C 
      CHARACTER*(*) NAME 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
C 
      IRET = 1 
      IF     (IACT .EQ. 0) THEN 
              WRITE ( NINT (VALUE) , 101) 
  101         FORMAT (/ ' Initiation parameters:   Qini, Ipd0, Ihdn ' / 
     >                  ' Maximum Q, Order of Alpha:     Qmax, IKNL ' / 
     >                  ' X- mesh parameters   :   Xmin, Xcr,   Nx  ' / 
     >                  ' LnQ-mesh parameters  :         Nt,   Jt   ' / 
     >                  ' # of parton flavors  :         NfMx       ' /) 
              IRET = 4 
      ElseIF (IACT .EQ. 1) THEN 
              CALL EVLSET (NAME, VALUE, IRET) 
      ElseIF (IACT .EQ. 2) THEN 
              CALL EVLGET (NAME, VALUE, IRET) 
      ElseIF (IACT .EQ. 5) THEN 
              CALL EVLGT1 (NAME, VALUE, IRET) 
      ElseIF (IACT .EQ. 3) THEN 
              CALL EVLIN 
              IRET = 4 
      ElseIF (IACT .EQ. 4) THEN 
              CALL EVLOUT ( NINT (VALUE) ) 
              IRET = 4 
      ElseIF (IACT .EQ. 6) THEN 
              CALL EVLOT1 ( NINT (VALUE) ) 
              IRET = 4 
      Else 
              IRET = 3 
      EndIf 
C 
      RETURN 
C                        **************************** 
      END 
 
      SUBROUTINE EVLGET (NAME, VALUE, IRET) 
C 
C                                           Gets VALUE of variable named NAME. 
C 
C               IRET =   0      variable not found. 
C                        1      success. 
C 
C               NAME is assumed upper-case, and VALUE real. 
C               If necessary, VALUE is converted to integer by NINT 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      LOGICAL LSTX 
C 
      CHARACTER*(*) NAME 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (MXPN = MXF * 2 + 2) 
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN) 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
 
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX, DEL,IV 
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG 
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx 
C 
      IRET = 1 
C 
      IF     (NAME .EQ. 'QINI')  THEN 
          VALUE = QINI 
      ElseIF (NAME .EQ. 'IPD0')  THEN 
          VALUE = IPD0 
      ElseIF (NAME .EQ. 'IHDN') THEN 
          VALUE = IHDN 
      ElseIF (NAME .EQ. 'QMAX')  THEN 
          VALUE = QMAX 
      ElseIF (NAME .EQ. 'IKNL') THEN 
          VALUE = IKNL 
      ElseIF (NAME .EQ. 'XCR') THEN 
          VALUE = XCR 
      ElseIF (NAME .EQ. 'XMIN') THEN 
          VALUE = XMIN 
      ElseIF (NAME .EQ. 'DEL') THEN 
          VLAUE = DEL
      ElseIF (NAME .EQ. 'NX') THEN 
          VALUE = NX 
      ElseIF (NAME .EQ. 'NT') THEN 
          VALUE = NT 
      ElseIF (NAME .EQ. 'JT') THEN 
          VALUE = JT 
      ElseIF (NAME .EQ. 'NFMX') THEN 
          VALUE = NfMx 
      Else 
          IRET = 0 
      EndIf 
C 
      RETURN 
C                       ____________________________ 
C 
      ENTRY EVLOUT (NOUUT) 
C 
C                             Write current values of parameters to unit NOUUT 
C 
      WRITE (NOUUT, 131) QINI, IPD0, IHDN, QMAX, IKNL, 
     >      XMIN, XCR, NX, NT, JT, NfMx 
C 
  131 FORMAT ( / 
     >' Current parameters and values are: '// 
     >' Initiation parameters: Qini, Ipd0, Ihdn = ', F8.2, 2I8 // 
     >' Maximum Q, Order of Alpha:   Qmax, IKNL = ', 1PE10.2, I6// 
     >' X- mesh parameters   : Xmin, Xcr,  Nx   = ', 2(1PE10.2), I8// 
     >' LnQ-mesh parameters  : Nt,   Jt         = ', 2I8       // 
     >' # of parton flavors  : NfMx             = ',  I8       /) 
C 
      RETURN 
C                        **************************** 
      END 
 
      SUBROUTINE EVLGT1 (NAME, VALUE, IRET) 
C 
C                                       COPY OF EVLGET for the alternate set 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C 
      LOGICAL LSTX 
      CHARACTER*(*) NAME 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (MXPN = MXF * 2 + 2) 
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN) 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
 
      COMMON / X1RRAY / XCR, XMIN, XV(0:MXX), LSTX, NX 
      COMMON / Q1RAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG 
      COMMON / E1LPAR / AL, IKNL, IPD0, IHDN, NfMx 
C 
      IRET = 1 
C 
      IF     (NAME .EQ. 'QINI')  THEN 
          VALUE = QINI 
      ElseIF (NAME .EQ. 'IPD0')  THEN 
          VALUE = IPD0 
      ElseIF (NAME .EQ. 'IHDN') THEN 
          VALUE = IHDN 
      ElseIF (NAME .EQ. 'QMAX')  THEN 
          VALUE = QMAX 
      ElseIF (NAME .EQ. 'IKNL') THEN 
          VALUE = IKNL 
      ElseIF (NAME .EQ. 'XCR') THEN 
          VALUE = XCR 
      ElseIF (NAME .EQ. 'XMIN') THEN 
          VALUE = XMIN 
      ElseIF (NAME .EQ. 'DEL') THEN 
         VALUE = DEL
      ElseIF (NAME .EQ. 'NX') THEN 
          VALUE = NX 
      ElseIF (NAME .EQ. 'NT') THEN 
          VALUE = NT 
      ElseIF (NAME .EQ. 'JT') THEN 
          VALUE = JT 
      ElseIF (NAME .EQ. 'NFMX') THEN 
          VALUE = NfMx 
      Else 
          IRET = 0 
      EndIf 
C 
      RETURN 
C                       ____________________________ 
C 
      ENTRY EVLOT1 (NOUUT) 
C 
C                             Write current values of parameters to unit NOUUT 
C 
      WRITE (NOUUT, 131) QINI, IPD0, IHDN, QMAX, IKNL, 
     >      XMIN, XCR, NX, NT, JT, NfMx 
C 
  131 FORMAT ( / 
     >' Current parameters and values are: '// 
     >' Initiation parameters: Qini, Ipd0, Ihdn = ', F8.2, 2I8 // 
     >' Maximum Q, Kernel ID : Qmax, IKNL       = ', 1PE10.2, I6// 
     >' X- mesh parameters   : Xmin, Xcr,  Nx   = ', 2(1PE10.2), I8// 
     >' LnQ-mesh parameters  : Nt,   Jt         = ', 2I8       // 
     >' # of parton flavors  : NfMx             = ',  I8       /) 
C 
      RETURN 
C                        **************************** 
      END 
C 
 
      SUBROUTINE EVLIN 
C 
C                                       Solicits parameters in EVL calculations 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
 
      CALL EVLGET ('QINI', V1, IR1) 
      CALL EVLGET ('IPD0', V2, IR2) 
      CALL EVLGET ('IHDN', V3, IR3) 
 
    1 WRITE  (NOUT, 101) V1, V2, V3 
  101 FORMAT (/' Current values of parameters QINI, IPD0, IHDN: ', 
     >        1PE10.2, 2I6 / '$Type new values: ' ) 
      READ(NIN, *, ERR=302) V1, V2, V3 
C 
      CALL EVLSET ('QINI', V1, IRET1) 
      CALL EVLSET ('IPD0', V2, IRET2) 
      CALL EVLSET ('IHDN', V3, IRET3) 
      goto 301 
 302  write(nout,120) 
      goto 1 
C 
 301  IF ((IRET1.NE.1) .OR. (IRET2.NE.1) .OR. (IRET3.NE.3)) THEN 
   20   WRITE (NOUT, 120) 
        GOTO 1 
      EndIf 
 
      CALL EVLGET ('QMAX', V1, IR1) 
      CALL EVLGET ('NT',   V2, IR2) 
      CALL EVLGET ('JT',   V3, IR3) 
C 
    2 WRITE  (NOUT, 102)  V1, V2, V3 
  102 FORMAT (/' Current values of parameters QMAX, NT, JT: ', 
     >        1PE10.2, I6, I6 / '$Type new values: ' ) 
      READ(NIN, *, ERR=304) V1, V2, V3 
C 
      CALL EVLSET ('QMAX', V1, IRET1) 
      CALL EVLSET ('NT',   V2, IRET2) 
      CALL EVLSET ('JT',   V3, IRET3) 
      goto 303 
 304  write(nout,120) 
      goto 2 
C 
 303  IF ((IRET1.NE.1) .OR. (IRET2.NE.1) .OR. (IRET3.NE.1)) THEN 
   22   WRITE (NOUT, 120) 
        GOTO 2 
      EndIf 
C 
      CALL EVLGET ('XMIN', V1, IR1) 
      CALL EVLGET ('XCR',  V2, IR2) 
      CALL EVLGET ('NX',   V3, IR3) 
      CALL EVLGET ('NFMX',   V4, IR4) 
      CALL EVLGET ('IKNL', V5, IR5) 
 
    3 WRITE  (NOUT, 103)  V1, V2, V3, V4, V5 
  103 FORMAT(/' Current values of parameters XMIN, XCR, NX,NFMX,IKNL: ', 
     >  2(1PE12.3), 3I6 / '$Type new values: ' ) 
      READ(NIN, *, ERR=22) V1, V2, V3, V4, V5 
C 
      CALL EVLSET ('XMIN', V1, IRET1) 
      CALL EVLSET ('XCR',  V2, IRET2) 
      CALL EVLSET ('NX',   V3, IRET3) 
      CALL EVLSET ('NFMX',   V4, IRET4) 
      CALL EVLSET ('IKNL', V5, IRET5) 
C 
      IF ( (IRET1 .NE. 1) .OR. (IRET2 .NE. 1) .OR. (IRET3 .NE. 1) 
     > .OR. (IRET4 .NE. 1) .OR. (IRET5 .NE. 1) ) THEN 
   23   WRITE (NOUT, 120) 
        GOTO 3 
      EndIf 
C 
  120 FORMAT(' Bad values, Try again!' /) 
C 
      RETURN 
C                        **************************** 
      END 
 
      SUBROUTINE EVLSET (NAME, VALUE, IRET) 
C 
C                               Sets variable named NAME to VALUE. 
C 
C               IRET =   0      variable not found. 
C                        1      success. 
C                        2      variable found, but bad value. 
C 
C               NAME is assumed upper-case, and VALUE real. 
C               If necessary, VALUE is converted to integer by NINT 
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
      LOGICAL LSTX 
C 
      CHARACTER*(*) NAME 
C 
      PARAMETER (MXX = 1050, MXQ = 25, MXF = 6) 
      PARAMETER (MXPN = MXF * 2 + 2) 
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN) 
C 
      COMMON / IOUNIT / NIN, NOUT, NWRT 
 
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX, DEL,IV 
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG 
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx 
      
C 
      IRET = 1 
      IF     (NAME .EQ. 'QINI')  THEN 
          IF (VALUE .LE. 0) GOTO 12 
          QINI = VALUE 
      ElseIF (NAME .EQ. 'IPD0')  THEN 
          ITEM = NINT(VALUE) 
          ITM = MOD (ITEM, 20) 
          IF (ITM .LT. 0 .OR. ITM .GT. 9) GOTO 12 
          IPD0 = ITEM 
C 
      ElseIF (NAME .EQ. 'IHDN') THEN 
          ITEM = NINT(VALUE) 
          IF (ITEM .LT. -1 .OR. ITEM .GT. 5) GOTO 12 
          IHDN = ITEM 
      ElseIF (NAME .EQ. 'QMAX')  THEN 
          IF (VALUE .LE. QINI) GOTO 12 
          QMAX = VALUE 
      ElseIF (NAME .EQ. 'IKNL') THEN 
          ITMP = NINT(VALUE) 
          ITEM = ABS(ITMP) 
          IF (ITEM .NE. 1 .AND. ITEM .NE. 2) GOTO 12 
          IKNL = ITMP 
      ElseIF (NAME .EQ. 'XCR') THEN 
          IF (VALUE .LT. XMIN .OR. VALUE .GT. 10.) GOTO 12 
          XCR = VALUE 
          LSTX = .FALSE. 
      ElseIF (NAME .EQ. 'XMIN') THEN 
          IF (VALUE .LT. 1D-7 .OR. VALUE .GT. 1D0) GOTO 12 
          XMIN = VALUE 
          LSTX = .FALSE. 
      ElseIF (NAME .EQ. 'DEL') THEN 
          IF (VALUE .LT. 0D0 .OR. VALUE .GT. 1D0) GOTO 12 
          DEL = VALUE 
          LSTX = .FALSE. 
      ElseIF (NAME .EQ. 'NX') THEN 
          ITEM = NINT(VALUE) 
          IF (ITEM .LT. 10 .OR. ITEM .GT. MXX-1) GOTO 12 
          NX = ITEM 
          LSTX = .FALSE. 
      ElseIF (NAME .EQ. 'NT') THEN 
          ITEM = NINT(VALUE) 
          IF (ITEM .LT. 2 .OR. ITEM .GT. MXQ) GOTO 12 
          NT = ITEM 
      ElseIF (NAME .EQ. 'JT') THEN 
          ITEM = NINT(VALUE) 
          IF (ITEM .LT. 1 .OR. ITEM .GT. 5) GOTO 12 
          JT = ITEM 
      ElseIF (NAME .EQ. 'NFMX') THEN 
          ITEM = NINT(VALUE) 
          IF (ITEM .LT. 1 .OR. ITEM .GT. MXPN) GOTO 12 
          NfMx = ITEM 
      Else 
          IRET = 0 
      EndIf 
C 
      RETURN 
C                                                                  Error exit: 
   12 IRET = 2 
C 
      RETURN 
C                        **************************** 
      END  

C program qcdpac

      Function SetQCD ()


C   #  filelist 
C   g.f
C   alpi.f
C   alpior.f
C   alepi.f
C   anom.f
C   sud.f
C   evoluf.f
C   alfset.f
C   rtalf.f
C   setlam.f
C   cnvl1.f
C   zcnvlm.f
C   parqcd.f
C   qcdin.f
C   qcdout.f
C   qcdset.f
C   qcdget.f
C   namqcd.f
C   alambd.f
C   nfl.f
C   alamf.f
C   amhatf.f
C   amumin.f
C   amass.f
C   ch.f
C   nfltot.f
C   alphem.f
C   lamcwz.f
C   setl1.f
C   trnlam.f
C   zbrlam.f
C   alpqcd.f
C   qzbrnt.f
C   mtmul.f
C   qwarn.f

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      External DatQCD

      SetQCD = 0.

      Return
C                        ****************************
      END

      BLOCK DATA DATQCD
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /COMQCH/ VALQCH(9)
      COMMON /COMQMS/ VALQMS(9)
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      COMMON /COMALP/ ALPHA
      LOGICAL SET
C
      DATA AL, NF, NORDER, SET / .20, 5, 2, .FALSE. /
      DATA VALQCH/ 0.66666667, -0.33333333,
     >  -0.33333333, 0.66666667,
     >  -0.33333333, 0.66666667,
     >  3*0./
      DATA VALQMS/  2*0., .2, 1.6, 5., 180., 3*0./
      DATA ALPHA/  7.29927E-3 /
 
      END
C
C******************************
C
      FUNCTION G (AMU)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C               Returns effective coupling.
      PARAMETER (PI = 3.1415927)
      G=2.*PI*SQRT(ALPI(AMU))
      RETURN
      END
C
C******************
C
      FUNCTION ALPI (AMU)
C               Returns effective g**2/(4pi**2) = alpha/pi.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET
C                      Use the following as subroutine argument of type
C                      set by IMPLICIT statement:
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
 
      DATA IW1, IW2 / 2*0 /
C
      IF(.NOT.SET) CALL LAMCWZ
 
      NEFF = NFL(AMU)
      ALM  = ALAM(NEFF)
      ALPI = ALPQCD (NORDER, NEFF, AMU/ALM, IRT)
 
      IF (IRT .EQ. 1) THEN
         CALL QWARN (IW1, NWRT, 'AMU < ALAM in ALPI', 'MU', AMU,
     >              ALM, BIG, 1)
         WRITE (NWRT, '(A,I4,F15.3)') 'NEFF, LAMDA = ', NEFF, ALM
      ELSEIF (IRT .EQ. 2) THEN
         CALL QWARN (IW2, NWRT, 'ALPI > 1; Be aware!', 'ALPI', ALPI,
     >             D0, D1, 0)
         WRITE (NWRT, '(A,I4,2F15.3)') 'NF, LAM, MU= ', NEFF, ALM, AMU
      ENDIF
 
      RETURN
      END
C
C************
C
      FUNCTION ALPIOR (AMU, NL)
C               Returns effective g**2/(4pi**2) = alpha/pi.
C               Use formula with NL loops for beta function.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET
 
C                      Use the following as subroutine argument of type
C                      set by IMPLICIT statement:
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
 
      DATA IW1, IW2 / 2*0 /
 
      IF (.NOT.SET) CALL LAMCWZ
 
      NEFF = NFL(AMU)
      ALM  = ALAM(NEFF)
      ALPIOR = ALPQCD (NL, NEFF, AMU/ALM, IRT)
 
      IF (IRT .EQ. 1) THEN
         CALL QWARN (IW1, NWRT, 'AMU < ALAM in ALPIOR', 'MU', AMU,
     >              ALM, BIG, 1)
         WRITE (NWRT, '(A,I4,F15.3)') 'NEFF, LAMDA = ', NEFF, ALM
      ELSEIF (IRT .EQ. 2) THEN
         CALL QWARN (IW2,NWRT,'ALPIOR > 1; Be aware!','ALPIOR',ALPIOR,
     >             D0, D1, 0)
         WRITE (NWRT, '(A,I4,2F15.3)') 'NF, LAM, MU= ', NEFF, ALM, AMU
      ENDIF
      END
C
C*****************************************************************
C
      FUNCTION ALEPI (AMU, NEF)
C                   Returns ALPHA/PI using the Effective Lamda appropriate for
C                             NEF flavors without regard to the value of AMU.
C                   Appropriate for Renormalization Schemes with fixed NEF.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
 
      DATA IW1, IW2 / 2*0 /
C
      IF(.NOT.SET) CALL LAMCWZ
 
      ALM = ALAM(NEF)
      ALEPI = ALPQCD (NORDER, NEF, AMU/ALM, IRT)
 
      IF     (IRT .EQ. 1) THEN
         CALL QWARN (IW1, NWRT, 'AMU < ALAM in ALEPI', 'MU', AMU,
     >             ALM, BIG, 1)
         WRITE (NWRT, '(A,I4,F15.3)') 'NEFF, LAMDA = ', NEF, ALM
      ELSEIF (IRT .EQ. 2) THEN
         CALL QWARN (IW2, NWRT, 'ALPI > 1; Be aware!', 'ALEPI', ALEPI,
     >             D0, D1, 0)
         WRITE (NWRT, '(A,I4,2F15.3)') 'NF, LAM, MU= ', NEF, ALM, AMU
      ENDIF
 
      RETURN
      END
C
C                      **************************************
C
      FUNCTION ANOM(Q1,Q2,GARRAY,N)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Returns integral from Q1 to Q2 of
C                 (dmu/mu) * gam(g(mu))
C              where N terms in anomalous dim. gam(g) are used:
C                 gam(g) = sum(i=1 to N) of gamma(i)*((g**2/(4*pi**2))**i).
C                 GAMMA(i) = sum(j=1 to i) of gammnf(i,j)*(nlfl**(j-1))
C                   except GAMMA(1)=gammnf(1,1)+gammnf(1,2)*nlfl
C                 gammnf(i,j) = GARRAY(1+j+i*(i-1)/2)
C                   except gammnf(1,1) = GARRAY(1) & gammnf(1,2)= GARRAY(2)
C                       where nlfl is the number of "light" flavor.
       DIMENSION GAMMA(5),GARRAY(10)
       COMMON /QCDPAR/ AL, NF, NORDER, SET
       COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
       COMMON /IOUNIT/ NIN, NOUT, NWRT
       LOGICAL SET
C
      IF(.NOT.SET) CALL LAMCWZ
      ANOM=0.
      IF (N.LE.0.) RETURN
      IF ((Q1.LE.AMN).OR.(Q2.LE.AMN)) THEN
         WRITE (NOUT, *) 'Q1 OR/AND Q2 IS TOO SMALL IN ANOM'
         RETURN
      ENDIF
C
      NMMIN=NF-NFL(Q2)
      NMMAX=NF-NFL(Q1)
C
  90  T1= LOG(Q2)
      B1=FLOAT(33-2*(NF-NMMIN+1))/12.
      GAMMA(1)=GARRAY(1)+GARRAY(2)*(NF-NMMIN+1)
      DO 200 J=NMMIN+1,NMMAX+1
         T2=T1
         B1=B1+1./6.
         GAMMA(1)=GAMMA(1)-GARRAY(2)
         IF (J.EQ.(NMMAX+1)) THEN
            T1= LOG(Q1)
         ELSE
            T1= LOG(AMHAT(NF-J+1))
         ENDIF
         ANOM=ANOM+0.5*GAMMA(1)* LOG((T2- LOG(ALAM(NF+1-J)))/
     >        (T1- LOG(ALAM(NF+1-J))))/B1
 200  CONTINUE
      RETURN
      END
C
C********************
C
      FUNCTION SUD(Q1,Q2,GARRAY,N)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Returns integral from Q1 to Q2 of
C                 (dmu/mu) * gam(g(mu)) * ln(Q2/mu)
C              where N terms in anomalous dim. gam(g) are used:
C                 gam(g) = sum(i=1 to N) of GAMMA(i)*((g**2/(4*pi**2))**i).
C                 GAMMA(i) = sum(j=1 t0 i) of gammnf(i,j)*(nlfl**(j-1))
C                 gammnf(i,j) = GARRAY(j+i*(i-1)/2)
C                        where nlfl is the number of "light" flavors.
C
       DIMENSION GAMMA(5),GARRAY(10)
       COMMON / QCDPAR / AL, NF, NORDER, SET
       COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
       COMMON / IOUNIT / NIN, NOUT, NWRT
       LOGICAL SET
       GAMMNF(K,L)=GARRAY(L+K*(K-1)/2)
C
      IF(.NOT.SET) CALL LAMCWZ
      SUD=0.
      IF(N.LE.0.) RETURN
      IF((Q1.LE.AMN).OR.(Q2.LE.AMN)) THEN
         WRITE(NOUT, *) 'Q1 OR/AND Q2 IS TOO SMALL IN SUD'
         RETURN
      ENDIF
C
      NMMIN=NF-NFL(Q2)
      NMMAX=NF-NFL(Q1)
      T10=2.* LOG(Q2)
      B1=FLOAT(33-2*(NF-NMMIN+1))/12.
      B2B1SQ=FLOAT(153-19*(NF-NMMIN+1))/(24.*B1*B1)
      GAMMA(2)=GAMMNF(2,1)+GAMMNF(2,2)*(NF-NMMIN+1)
      DO 200 J=NMMIN+1,NMMAX+1
         T20=T10
         B2B1SQ=(B2B1SQ*B1*B1+19./24.)/(B1+1./6.)**2
         B1=B1+1./6.
         IF (J.EQ.(NMMAX+1)) THEN
            T10=2.* LOG(Q1)
         ELSE
            T10=2.* LOG(AMHAT(NF-J+1))
         ENDIF
         ALLAM=2.* LOG(ALAM(NF+1-J))
         TQ2=2.* LOG(Q2)-ALLAM
         T1=T10-ALLAM
         T2=T20-ALLAM
         ALNT2= LOG(T2)
         ALNT1= LOG(T1)
         SUD=SUD+GARRAY(1)*0.25/B1*(TQ2*(ALNT2-ALNT1
     1             +B2B1SQ*((ALNT2+1.)/T2-(ALNT1+1.)/T1))
     2             +T1-T2+B2B1SQ*(ALNT2**2-ALNT1**2)/2.)
C
         IF (N.GE.2) THEN
            GAMMA(2)=GAMMA(2)-GAMMNF(2,2)
            SUD=SUD+0.25*GAMMA(2)/(B1*B1)*(ALNT1-ALNT2
     >                                +TQ2*(1./T1-1./T2))
         ENDIF
200   CONTINUE
      RETURN
      END
C
C**********************
C
      SUBROUTINE EVOLUF(FQ2,FQ1,Q2,Q1,N)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                        ----------------------------
C?????????????????????????TO CHECK ????????????????????????????
C                        ----------------------------
C            Returns moments of parton distribution function at Q2 from Q1.
C            by integrating
C                D(F(Q,N,I))/D(LN(Q**2))=-GAMA(I,J)*F(Q,N,J)
C            where N is rank of the moments.
C            GAMA is aN NF X NF  matrix.
C            F(I=1) corresponds to gluon.
C            F(I=2 to NF+1) corresponds to NF quarks.
C
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      LOGICAL SET
      DIMENSION FQ1(11),FQ2(11),U(11,11),UINV(11,11),
     >          GAMAD(11,11),TEMP(11,11)
      DATA U,UINV,GAMAD/363*0./
C
      IF(.NOT.SET) CALL LAMCWZ
      IF((Q1.LE.AMN).OR.(Q2.LE.AMN)) GOTO 300
      IF(N.LT.2) GOTO 600
C
      NMMIN=NF-NFL(Q2)
      NMMAX=NF-NFL(Q1)
C
      T2=2.* LOG(Q1)
      B1=FLOAT(33-2*(NF-NMMAX-1))/4.
      SUM=0.
      DO 10 I=2,N
         SUM=SUM+1./FLOAT(I)
 10   CONTINUE
      GAGG=.75-9./FLOAT((N+1)*(N+2))-9./FLOAT(N*(N-1))
     >        +FLOAT(NF-NMMAX-1)/2.+SUM*9.
      GAFG=3./FLOAT(N)-6./FLOAT((N+1)*(N+2))
      GAGF=1./FLOAT(N+1)+2./FLOAT(N*(N-1))
      GAFF=1.-2./FLOAT(N*(N+1))+4.*SUM
      FQ2(1)=FQ1(1)
      DO 20 I=2,NF+1
         DO 30 J=3,NF+1
            UINV(J,I)=1./FLOAT(NF)
  30     CONTINUE
         FQ2(I)=FQ1(I)
         UINV(I,I)=UINV(I,I)-1.
         U(I,I)=-1.
         U(2,I)=1.
         U(I,1)=1.
         U(I,2)=1.
 20   CONTINUE
C
      DO 100 K=NMMAX+1,NMMIN+1,-1
         T1=T2
         B1=B1-.5
         GAGG=GAGG+.5
         IF(K.EQ.(NMMIN+1)) THEN
            T2=2.* LOG(Q2)
         ELSE
            T2=2.* LOG(AMHAT(NF-K+2))
         ENDIF
         A=SQRT((GAGG-GAFF)**2+4.*NF*GAGF*GAFG)
         U(1,1)=.5*(GAGG-GAFF+A)
         U(1,2)=.5*(GAGG-GAFF-A)
         DO 40 L=2,NF+1
            UINV(1,L)=-U(1,2)/A/FLOAT(NF)
            UINV(2,L)=U(1,1)/A/FLOAT(NF)
  40     CONTINUE
         UINV(1,1)=1./A
         UINV(1,2)=-1./A
         D=2* LOG(ALAM(NF+1-K))
         C=(T2-D)/(T1-D)
         DO 50 M=3,NF+1
            GAMAD(M,M)=C**(-GAFF/B1)
  50     CONTINUE
         GAMAD(1,1)=C**(-(U(1,1)+GAFF)/B1)
         GAMAD(2,2)=C**(-(U(1,2)+GAFF)/B1)
         CALL MTMUL(NF+1,NF+1,NF+1,U,GAMAD,TEMP)
         CALL MTMUL(NF+1,NF+1,NF+1,TEMP,UINV,TEMP)
         CALL MTMUL(NF+1,NF+1,1,TEMP,FQ2,FQ2)
 100  CONTINUE
      RETURN
C
300   WRITE (NOUT,400)
400   FORMAT('Q1 OR/AND Q2 IS TOO SMALL IN EVOLUF')
600   DO 500 I=1,11
500      FQ2(I)=0.
      RETURN
      END
C
C                             *************************
C
      SUBROUTINE ALFSET (QS, ALFS)
 
C                                 Given the value of Alpha(strong), ALFS, at
C                                 the scale QS, this routine determines the
C                                 effective # of flavors and Effective Lamda
C                                 at QS appropriate for the current value of
C                                 NORDER; and then call SET1 to setup the
C                                 whole package for subsequent use.
C                                 Calculates Alpha_s at loop order
C                                 specified by current value of NORDER.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      EXTERNAL RTALF
      COMMON / RTALFC / ALFST, JORD, NEFF
 
      DATA ALAM, BLAM, ERR / 0.01, 10.0, 0.02 /
 
      QST   = QS
      ALFST = ALFS
      CALL PARQCD (2, 'ORDR', ORDR, IR1)
      JORD  = ORDR
 
      NEFF = NFL(QS)
 
      EFLLN  = QZBRNT (RTALF, ALAM, BLAM, ERR, IR2)
      EFFLAM = QS / EXP (EFLLN)
 
      CALL SETL1 (NEFF, EFFLAM)
 
      END
 
C
C**************************************************************
C
      FUNCTION RTALF (EFLLN)
C  Auxiliary function for ALFSET, which solves equation RTALF=0.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (PI = 3.141592653589)
      COMMON / RTALFC / ALFST, JORD, NEFF
 
      EFMULM = EXP (EFLLN)
      TEM1 = PI / ALFST
      TEM2 = 1. / ALPQCD (JORD, NEFF, EFMULM, I)
 
      RTALF = TEM1 - TEM2
 
      END
C
C************************************************************
C
      SUBROUTINE SETLAM (NEF, WLAM, IRDR)
C     The values of LAMBDA=WLAM with NEF effective flavors is given. The
C     coupling is assumed to be given by the IRDR formula.
C     First lambda is converted to a value that gives approximately the
C     same value of alpha_s when the formula with the current value of
C     NORDER is used.  Then SETL1 is called to update the rest of the
C     internal arrays.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET
 
      IF ((NEF .LT. 0) .OR. (NEF .GT. NF)) THEN
         WRITE(NOUT,*)'NEF out of range in SETLAM, NEF, NF=', NEF,NF
         STOP
      ENDIF
C                         Adjust Lamda value if ORDER parameters don't match.
C                         NORDER is not changed.
      VLAM = WLAM
      IF (IRDR .NE. NORDER) CALL CNVL1 (IRDR, NORDER, NEF, VLAM)
      CALL SETL1 (NEF, VLAM)
      END
C
C************************************************************
C
      SUBROUTINE CNVL1 (IRDR, JRDR, NF, VLAM)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL ZCNVLM
C  Auxiliary routine for SETLAM.
C                  Given Lamda (NF) = VLAM at order IRDR, this subroutine
C                  finds the corresponding Lamda at order JRDR so that the
C                  resulting Alpha(Mu, NRDR) remains approximately the same
C                  within the range of Neff = NF.
 
      COMMON / LAMCNV / AMU, ULAM, NFL, IRD, JRD
      COMMON / IOUNIT / NIN, NOUT, NWRT
 
      DATA ALM, BLM, ERR, AMUMIN / 0.001, 2.0, 0.02, 1.5 /
 
      IRD = IRDR
      JRD = JRDR
      ULAM = VLAM
 
      CALL PARQCD(2, 'NFL', ANF, IRT)
      NTL = NFLTOT()
      IF (NF .GT. NTL) THEN
         WRITE (NOUT, *) ' NF .GT. NTOTAL in CNVL1; set NF = NTOTAL'
         WRITE (NOUT, *) ' NF, NTOTAL = ', NF, NTL
         NF = NTL
      ENDIF
C                                                First match at NFth threshold
      NFL = NF
      AMU = AMHATF(NF)
      AMU = MAX (AMU, AMUMIN)
      VLM1 = QZBRNT (ZCNVLM, ALM, BLM, ERR, IR1)
C                                          Match again at the next threshold
      IF (NF .LT. NTL) THEN
        AMU = AMHATF(NF+1)
        AMU = MAX (AMU, AMUMIN)
        VLM2 = QZBRNT(ZCNVLM, ALM, BLM, ERR, IR2)
      ELSE
        VLM2 = VLM1
      ENDIF
C                              Take the average and return new value of VLAM
      VLAM = (VLM1 + VLM2) / 2
 
      RETURN
C                        ****************************
      END
 
      FUNCTION ZCNVLM (VLAM)
C  Auxiliary function for CNVL1, which solves ZCNVLAM=0.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      COMMON / LAMCNV / AMU, ULAM, NFL, IRD, JRD
 
      ZCNVLM= ALPQCD (IRD,NFL,AMU/ULAM,I) - ALPQCD (JRD,NFL,AMU/VLAM,I)
 
      END
C
C**************************************************************
C
      SUBROUTINE PARQCD(IACT,NAME,VALUE,IRET)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C             Actions: 0     type list of variables on unit VALUE.
C                      1     set variable with name NAME to VALUE, if
C                             it exists, else set IRET to 0.
C                      2     find value of variable. If it does not exist,
C                             set IRET to 0.
C                      3     request values of all parameters from terminal.
C                      4     type list of all values on unit VALUE
C
C             IRET =   0     variable not found.
C                      1     successful search
C                      2     variable found, but bad value.
C                      3     bad value for IACT.
C                      4     no variable search (i.e., IACT is 0,3,or 4).
C
C             NAME is assumed upper-case.
C             if necessary, VALUE is converted to integer by NINT(VALUE)
C
      INTEGER IACT,IRET
      CHARACTER*(*) NAME
C
      IRET=1
      IF (IACT.EQ.0) THEN
         WRITE (NINT(VALUE), *)  'LAM(BDA), NFL, ORD(ER), Mi, ',
     >               '(i in 1 to 9), LAMi (i in 1 to NFL)'
         IRET=4
      ELSEIF (IACT.EQ.1) THEN
         CALL QCDSET (NAME,VALUE,IRET)
      ELSEIF (IACT.EQ.2) THEN
         CALL QCDGET (NAME,VALUE,IRET)
      ELSEIF (IACT.EQ.3) THEN
         CALL QCDIN
         IRET=4
      ELSEIF (IACT.EQ.4) THEN
         CALL QCDOUT(NINT(VALUE))
         IRET=4
      ELSE
         IRET=3
      ENDIF
 
      RETURN
      END
C
C*******************************************************************
C
      SUBROUTINE QCDIN
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                 Inputs QCD parameters.
C
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      DIMENSION VALMAS(9)
C              ==== MACHINE DEPENDENCE
C              ==== ASSUMES CHARACTER CODES FOR '0' TO '9' ARE CONSECUTIVE
C              NOTE THAT ICHAR('0') IS NON-STANDARD FORTRAN (ACCORDING TO
C              MICROSOFT).
      CHARACTER ONECH*(1)
      ONECH = '0'
      IASC0 = ICHAR(ONECH)
C
      CALL QCDGET ('LAM',ALAM,IRET1)
      CALL QCDGET ('NFL',ANF,IRET2)
      CALL QCDGET ('ORDER',ORDER,IRET3)
      NF = NINT(ANF)
      NORDER = NINT(ORDER)
 1    WRITE (NOUT, *) 'LambdaMSBAR, # Flavors, loop order ?'
      READ (NIN,*, IOSTAT = IRET) ALAM, NF, NORDER
      ORDER = NORDER
      ANF = NF
      IF (IRET .LT. 0) GOTO 22
      IF (IRET .EQ. 0) THEN
         CALL QCDSET ('LAM',ALAM,IRET1)
         CALL QCDSET ('NFL',ANF,IRET2)
         CALL QCDSET ('ORDER',ORDER,IRET3)
         ENDIF
      IF ((IRET.NE.0) .OR. (IRET1.NE.1) .OR. (IRET2.NE.1)
     >     .OR. (IRET3.NE.1)) THEN
         WRITE (NOUT, *) 'Bad value(s), try again.'
         GOTO 1
         ENDIF
      DO 20, I = 1, NF
         CALL QCDGET('M'//CHAR(I+IASC0),VALMAS(I),IRET1)
 10      WRITE (NOUT, '(1X,A,I2,A)') 'Mass of Quark', I, '?'
         READ (NIN,*, IOSTAT=IRET) VALMAS(I)
         IF (IRET .LT. 0) GOTO 22
         IF (IRET .EQ. 0)
     >      CALL QCDSET('M'//CHAR(I+IASC0),VALMAS(I),IRET1)
         IF ((IRET .NE. 0) .OR. (IRET1 .NE. 1)) THEN
            WRITE (NOUT, *) 'Bad value, try again.'
            GOTO 10
            ENDIF
 20      CONTINUE
      RETURN
C
 22   WRITE (NOUT, *) 'END OF FILE ON INPUT'
      WRITE (NOUT, *)
      RETURN
      END
C
C *****
C
      SUBROUTINE QCDOUT(NOUT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Prints out the values of parameters to unit NOUT.
C
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON /COMQMS/ VALQMS(9)
      LOGICAL SET
C
      IF (.NOT. SET) CALL LAMCWZ
      WRITE (NOUT,110) AL, NF, NORDER
 110   FORMAT(
     1 ' Lambda (MSBAR) =',G13.5,', NFL (total # of Flavors) =',I3,
     2 ', Order (loops) =', I2)
      WRITE (NOUT,120) (I,VALQMS(I),I=1,NF)
 120   FORMAT (3(' M', I1, '=', G13.5, :, ','))
      IF (NHQ .GT. 0)
     1   WRITE (NOUT,130) (I, ALAMF(I), I = NF-NHQ, NF)
 130   FORMAT (: ' ! Effective lambda given number of light quarks:'/
     >    (2(' ! ', I1, ' quarks => lambda = ', G13.5 : '; ')) )
      RETURN
      END
C
C***************************************************************
C
      SUBROUTINE QCDSET (NAME,VALUE,IRET)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  Assign the variable whose name is specified by NAME the value VALUE
C             IRET=  0      variable not found.
C                    1      success.
C                    2      variable found, but bad value.
C             NAME is assumed upper-case and VALUE is real.
C             If necessary, VALUE is converted to integer by NINT(VALUE).
C
      CHARACTER*(*) NAME
      COMMON / COMQMS / VALMAS(9)
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET
      PARAMETER (PI=3.1415927, EULER=0.57721566)
C
      IVALUE = NINT(VALUE)
      ICODE  = NAMQCD(NAME)
      IF (ICODE .EQ. 0) THEN
         IRET=0
      ELSE
         IRET = 1
         SET = .FALSE.
         IF (ICODE .EQ. 1) THEN
            IF (VALUE.LE.0) GOTO 12
            AL=VALUE
         ELSEIF (ICODE .EQ. 2) THEN
            IF ( (IVALUE .LT. 0) .OR. (IVALUE .GT. 9)) GOTO 12
            NF = IVALUE
         ELSEIF ((ICODE .GE. 3) .AND. (ICODE .LE. 11))  THEN
            IF (VALUE .LT. 0) GOTO 12
            VALMAS(ICODE - 2) = VALUE
         ELSEIF ((ICODE .GE. 13) .AND. (ICODE .LE. 13+NF))  THEN
            IF (VALUE .LE. 0) GOTO 12
            CALL SETL1 (ICODE-13, VALUE)
         ELSEIF (ICODE .EQ. 24)  THEN
            IF ((IVALUE .LT. 1) .OR. (IVALUE .GT. 2)) GOTO 12
            NORDER = IVALUE
         ENDIF
         IF (.NOT. SET) CALL LAMCWZ
      ENDIF
      RETURN
C
C              Illegal value
 12   IRET=2
      RETURN
      END
C
C************************************************
C
      SUBROUTINE QCDGET(NAME,VALUE,IRET)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  Sets VALUE to the value of variable named NAME.
C            IRET=  0          variable not found.
C                   1          success.
C
C            NAME is assumed to be an upper-case character variable.
C
      CHARACTER*(*) NAME
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      COMMON / COMQMS / VALQMS(9)
      LOGICAL SET
      PARAMETER (PI=3.1415927, EULER=0.57721566)
C
      ICODE = NAMQCD(NAME)
      IRET = 1
      IF (ICODE .EQ. 1) THEN
         VALUE = AL
      ELSEIF (ICODE .EQ. 2) THEN
         VALUE = NF
      ELSEIF ((ICODE .GE. 3) .AND. (ICODE .LE. 12))  THEN
         VALUE = VALQMS(ICODE - 2)
      ELSEIF ((ICODE .GE. 13) .AND. (ICODE .LE. 13+NF))  THEN
         VALUE = ALAM(ICODE - 13)
      ELSEIF (ICODE .EQ. 24) THEN
         VALUE = NORDER
      ELSE
         IRET=0
      ENDIF
      END
C
C *****
C
      FUNCTION NAMQCD(NNAME)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                   Find integer code corresponding to NAME. If no
C                   match,then return zero.
C
      CHARACTER NNAME*(*), NAME*8
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      LOGICAL SET
C              ==== MACHINE DEPENDENCE
C              ==== ASSUMES CHARACTER CODES FOR '0' TO '9' ARE CONSECUTIVE
      CHARACTER ONECH*(1)
      ONECH = '0'
      IASC0 = ICHAR(ONECH)
C
C                Use local variable to avoid problems with short
C                passed argument:
      NAME = NNAME
      NAMQCD=0
      IF ( (NAME .EQ. 'ALAM') .OR. (NAME .EQ. 'LAMB') .OR.
     1        (NAME .EQ. 'LAM') .OR. (NAME .EQ. 'LAMBDA') )
     2             NAMQCD=1
      IF ( (NAME .EQ. 'NFL') .OR. (NAME(1:3) .EQ. '#FL') .OR.
     1        (NAME .EQ. '# FL') )
     2             NAMQCD=2
      DO 10 I=1, 9
         IF (NAME .EQ. 'M'//CHAR(I+IASC0))
     1             NAMQCD=I+2
10       CONTINUE
      DO 20 I= 0, NF
         IF (NAME .EQ. 'LAM'//CHAR(I+IASC0))
     1             NAMQCD=I+13
20       CONTINUE
      IF (NAME(:3).EQ.'ORD' .OR. NAME(:3).EQ.'NRD') NAMQCD = 24
      RETURN
      END
C
C
C***************************************************************
C
      FUNCTION ALAMBD(AMU)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      ALAMBD = ALAMF(NFL(AMU))
      RETURN
      END
C
C***************************************************************
C
      FUNCTION NFL(AMU)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C           NFL returns the number of 'light' flavors.
C
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      LOGICAL SET
C
      IF (.NOT. SET) CALL LAMCWZ
      NFL = NF - NHQ
      IF ((NFL .EQ. NF) .OR. (AMU .LE. AMN)) GOTO 20
      DO 10 I = NF - NHQ + 1, NF
         IF (AMU .GE. AMHAT(I)) THEN
            NFL = I
         ELSE
            GOTO 20
         ENDIF
10       CONTINUE
20    RETURN
      END
C
C***************************************************************
C
      FUNCTION ALAMF(N)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                 Returns the value of LambdaCWZ in the energy range
C                 with N "light" quarks.
 
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      COMMON / IOUNIT / NIN, NOUT, NWRT
      LOGICAL SET
C
      IF (.NOT.SET) CALL LAMCWZ
      IF ((N.LT.0) .OR. (N.GT.9)) THEN
         WRITE (NOUT, *) ' N IS OUT OF RANGE IN ALAMF'
         ALAMF=0.
      ELSE
         ALAMF = ALAM(MAX(N, NF-NHQ))
      ENDIF
      RETURN
      END
C
C***************************************************************
C
      FUNCTION AMHATF(I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                Returns the boundary in mass scale between the regions
C                with I & (I-1) effective "light" quarks.
C
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      COMMON / IOUNIT / NIN, NOUT, NWRT
      LOGICAL SET
C
      IF (.NOT.SET) CALL LAMCWZ
      IF ((I.LE.0).OR.(I.GT.9)) THEN
         WRITE (NOUT,*) 'I IS OUT OF RANGE IN AMHATF'
         AMHATF = 0
      ELSE
         AMHATF = AMHAT(I)
      ENDIF
      RETURN
      END
C
C********************************************************************
C
      FUNCTION AMUMIN()
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Returns the minimum mu allowed by perturbative theory.
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NFL, NORDER, SET
      LOGICAL SET
      IF (.NOT.SET) CALL LAMCWZ
      AMUMIN = AMN
      RETURN
      END
C
C************************************************************************
C
      FUNCTION AMASS(I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                     Returns mass of parton I.
C
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      COMMON /COMQMS/ VALMAS(9)
      LOGICAL SET
      AMASS = 0.
      IF (IABS(I) .GT. NF) THEN
          WRITE (NOUT, 20)
  20      FORMAT(' I IS OUT OF RANGE IN FUNCTION AMASS')
      ELSEIF (I .NE. 0)  THEN
          AMASS = VALMAS(IABS(I))
      ENDIF
      RETURN
      END
C
C***********************************************************
C
      FUNCTION CH(I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                     Returns charge of parton I.
C                        I=0 is gluon.
C                        I>0 is quark.
C                        I<0 is antiquark.
C                     See  BLOCK DATA for code.
C
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      COMMON /COMQCH/ VALCH(9)
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      LOGICAL SET
      CH=0.
      IF (IABS(I).GT.NF) THEN
          WRITE (NOUT, *) 'I IS OUT OF RANGE IN FUNCTION CH'
      ELSEIF (I.NE.0)  THEN
          CH = VALCH(IABS(I))
          IF (I.LT.0)  CH = -CH
      ENDIF
      RETURN
      END
C
C**********************************************************************
C
      FUNCTION NFLTOT()
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C          Returns the total number of flavors.
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      LOGICAL SET
      NFLTOT=NF
      RETURN
      END
C
C***********************************************************
C
      FUNCTION ALPHEM()
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Returns the value of electromagnetic interaction coupling.
      COMMON /COMALP/A
      ALPHEM=A
      RETURN
      END
C
C===================================================================
C===================================================================
C
      SUBROUTINE LAMCWZ
C                       Set /CWZPRM/ from /QCDPAR/
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET
      CALL SETL1 (NF, AL)
      END
C
C***********************
C
      SUBROUTINE SETL1  (NEF, VLAM)
C     Given LAMDA = VLAM for NEF flavors:
C                    (i) fills the array  ALAM (0:NF) with effective LAMDA's;
C                    (ii) fills the array AMHAT (NF) with threshold masses;
C                    (iii) count the # of "heavy quarks" (QMS > EFFLAM);
C                    (iv) fix the parameter AMN defined as MAX (ALAM),
C                         times safety factor;
C                    (v) set AL in / QCDPAR / equal to ALAM (NF);
C                    (vi) let SET = .TRUE.
C       Uses formula with NORDER (1 or 2) -- see /QCDPAR/
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      LOGICAL SET
 
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      COMMON / COMQMS / QMS(9)
      COMMON / IOUNIT / NIN, NOUT, NWRT
 
      IF (NEF .LT. 0 .OR. NEF .GT. NF) THEN
        WRITE(NOUT,*)'NEF out of range in SETL1, NEF, NF =',NEF,NF
        STOP
      ENDIF
C             Mass Thresholds are given by the Quark masses in the CWZ scheme
      AMHAT(0) = 0.
      DO 5 N = 1, NF
         AMHAT(N) = QMS(N)
    5    CONTINUE
      ALAM(NEF) = VLAM
      DO 10 N = NEF, 1, -1
         CALL TRNLAM(NORDER, N, -1, IR1)
   10    CONTINUE
      DO 20 N = NEF, NF-1
         CALL TRNLAM(NORDER, N, 1, IR1)
   20    CONTINUE
C=========================                Find first light quark:
      DO 30, N = NF, 1, -1
         IF ((ALAM(N) .GE. 0.7 * AMHAT(N))
     >       .OR. (ALAM(N-1) .GE. 0.7 * AMHAT(N)))THEN
            NHQ = NF - N
            GOTO 40
            ENDIF
   30    CONTINUE
      NHQ = NF
   40 CONTINUE
      DO 50, N = NF-NHQ, 1, -1
         AMHAT(N) = 0
         ALAM(N-1) = ALAM(N)
   50    CONTINUE
C========================               Find minimum mu
      AMN = ALAM(NF)
      DO 60, N = 0, NF-1
         IF (ALAM(N) .GT. AMN)  AMN = ALAM(N)
   60    CONTINUE
      AMN = AMN * 1.0001
      AL = ALAM(NF)
      SET = .TRUE.
      RETURN
      END
C
C**************************************************************
C
      SUBROUTINE TRNLAM (IRDR, NF, IACT, IRT)
 
C     This routine transforms LAMDA (N) to LAMDA (N+IACT) where IACT = 1/-1
C     The transformation is obtained by requiring the coupling constant to
C                be continuous at the scale Mu = Mass of the (N+1)th quark.
 
C                                         IRT is an return code.
C                                            (0 for OK)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / TRNCOM / VMULM, JRDR, N, N1
 
      EXTERNAL ZBRLAM
 
      DATA ALM0, BLM0, RERR / 0.01, 10.0, 0.0001 /
      DATA IW1, IW2, IW3, IW4, IR1, SML / 4*0, 0, 1.E-5 /
 
      IRT = 0
 
      N = NF
      JRDR = IRDR
      JACT = IACT
      VLAM = ALAM(N)
 
      IF (JACT .GT. 0) THEN
         N1 = N + 1
         THMS = AMHAT(N1)
         ALM = LOG (THMS/VLAM)
         BLM = BLM0
      ELSE
         N1 = N -1
         THMS = AMHAT(N)
         ALM = ALM0
         THMS = MAX (THMS, SML)
         BLM = LOG (THMS/VLAM)
      ENDIF
C                          Fix up for light quark:
      IF (VLAM .GE. 0.7 * THMS) THEN
         IF (JACT . EQ. 1) THEN
            AMHAT(N1) = 0
         ELSE
            AMHAT(N) = 0
         ENDIF
         IRT = 4
         ALAM(N1) = VLAM
         RETURN
      ENDIF
 
C             QZBRNT is the root-finding function to solve ALPHA(N) = ALPHA(N1)
C             Since 1/Alpha is roughly linear in Log(Mu/Lamda), we use the
C             former in ZBRLAM and the latter as the function variable.
      IF (ALM .GE. BLM) THEN
         WRITE (NOUT, *) 'TRNLAM has ALM >= BLM: ', ALM, BLM
         WRITE (NOUT, *) 'I do not know how to continue'
         STOP
         ENDIF
      VMULM = THMS/VLAM
      ERR = RERR * LOG (VMULM)
      WLLN = QZBRNT (ZBRLAM, ALM, BLM, ERR, IR1)
      ALAM(N1) = THMS / EXP (WLLN)
 
      IF (IR1 .NE. 0) THEN
         WRITE (NOUT, *) 'QZBRNT failed to find VLAM in TRNLAM; ',
     >        'NF, VLAM =', NF, VLAM
         WRITE (NOUT, *) 'I do not know how to continue'
        STOP
      ENDIF
      RETURN
      END
C                             *************************
 
      FUNCTION ZBRLAM (WLLN)
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / TRNCOM / VMULM, JRDR, N, N1
 
      WMULM = EXP (WLLN)
      TEM1 = 1./ ALPQCD(JRDR, N1, WMULM, I)
      TEM2 = 1./ ALPQCD(JRDR, N,  VMULM, I)
 
      ZBRLAM = TEM1 - TEM2
 
      END
 
 
C************************************************************
 
      FUNCTION ALPQCD (IRDR, NF, RML, IRT)
 
C                                 Returns the QCD alpha/pi for RML = MU / LAMDA
C                                 using the standard perturbative formula for
C                                 NF flavors and to IRDR th order in 1/LOG(RML)
 
C                                 Return Code:  IRT
C                                                0:   O.K.
C                                                1:   Mu < Lamda; returns 99.
C                                                2:   Alpha > 1 ; be careful!
C                                                3:   IRDR out of range
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
      PARAMETER (CG = 3.0, TR = 0.5, CF = 4.0/3.0)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT
 
      DATA IW1 / 0/
 
      IRT = 0
 
      IF (IRDR .LT. 1 .OR. IRDR .GT. 2) THEN
        WRITE(NOUT, *) 'Order parameter out of range in ALPQCD;
     >  IRDR = ', IRDR
        IRT = 3
        STOP
      ENDIF
 
      B0 = (11.* CG  - 2.* NF) / 3.

      B1 = (34.* CG**2 - 10.* CG * NF - 6.* CF * NF) / 3.

      RM2 = RML ** 2
 
      IF (RM2 .LE. 1.) THEN
         IRT = 1
C         CALL QWARN(IW1, NWRT,
C     >      'RM2 (=MU/LAMDA) < 1. not allowed in ALPQCD',
C     >      'RM2', RM2, D1, BIG, 1)
         ALPQCD = 99
         RETURN
      ENDIF
 
      ALN = LOG (RM2)
      AL = 4./ B0 / ALN
 
      IF (IRDR .GE. 2) AL = AL * (1.- B1 * LOG(ALN) / ALN / B0**2)
 
      IF (AL .GE. 1.) THEN
         IRT = 2
C     >   CALL QWARN(IW2, NWRT, 'ALPQCD > 1 in ALPQCD', 'ALPQCD', AL,
C     >              D0, D1, 1)
      ENDIF
 
      ALPQCD = AL
 
      RETURN
      END
C
C
C===================================================================
C===================================================================
C
C            Utilities for warning messages and interpolation here
C            to make self-contained package.
C
      FUNCTION QZBRNT(FUNC, X1, X2, TOLIN, IRT)
 
C                          Return code  IRT = 1 : limits do not bracket a root;
C                                             2 : function call exceeds maximum
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      PARAMETER (ITMAX = 1000, EPS = 3.E-12)

      external func
 
      TOL = ABS(TOLIN)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.)  THEN
        WRITE (NOUT, *) 'Root must be bracketed for QZBRNT.'
        IRT = 1
      ENDIF
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          QZBRNT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
      WRITE (NOUT, *) 'QZBRNT exceeding maximum iterations.'
      IRT = 2
      QZBRNT=B
      RETURN
      END
 
C
C**************************************************
C
      SUBROUTINE MTMUL(L,M,N,A,B,C)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  TO DO MATRIX MULTIPLICATION
      DIMENSION A(11,11),B(11,11),C(11,11)
      DO 20 I=1,L
      DO 20 K=1,N
      C(I,K)=0.
      DO 20 J=1,M
20       C(I,K)=C(I,K)+A(I,J)*B(J,K)
      RETURN
      END
C
C***********************************************************
C
      SUBROUTINE QWARN (IWRN, NWRT1, MSG, NMVAR, VARIAB,
     >                  VMIN, VMAX, IACT)
 
C     Subroutine to handle warning messages.  Writes the (warning) message
C     and prints out the name and value of an offending variable to SYS$OUT
C     the first time, and to output file unit # NWRT1 in subsequent times.
C
C     The switch IACT decides whether the limits (VMIN, VMAX) are active or
C     not.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
 
      CHARACTER*(*) MSG, NMVAR
 
      IW = IWRN
      VR = VARIAB
 
      WRITE (NWRT1,'(I5, 3X,A/ 1X,A,'' = '',1PD16.7)') IW, MSG,
     >                  NMVAR, VR
 
      IF  (IW .EQ. 0) THEN
         WRITE (NOUT, '(1X, A/1X, A,'' = '', 1PD16.7/A,I4)')
     >      MSG, NMVAR, VR,
     >      ' Complete set of warning messages on file unit #', NWRT1
         IF (IACT .EQ. 1) THEN
         WRITE (NOUT,'(1X,A/2(1PD15.3))')'The limits are: ', VMIN,VMAX
         WRITE (NWRT1,'(1X,A/2(1PD15.3))')'The limits are: ', VMIN,VMAX
         ENDIF
      ENDIF
 
      IWRN = IW + 1
 
      RETURN
C                         *************************
      END

C Subject: utlpac

      Function SetUTL ()


      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      External DatUTL

      SetUTL = 0.

      Return
C                        ****************************
      END

      BLOCK DATA DATUTL

      COMMON / IOUNIT / NIN, NOUT, NWRT

      DATA NIN, NOUT, NWRT / 5, 6, 56 /

C                         *************************
      END
C                         ==========================
C
      SUBROUTINE ADZ2AL (F,I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E15)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.

      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)

        CALL SGL2NT (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL SGL2NT (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
      ELSE
C                                                                   Closed endS
        FW(I) = F(W)
        TEM = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
C                                       Preliminary error Simpson - trapezoidal:
        ER  = DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.
      ENDIF

      RESULT(I) = TEM
      ERR   (I) = ABS (ER)

      RETURN
C                        ****************************
      END

      FUNCTION ADZ2NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
      DATA SMLL / 1E-20 /

      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in ADZ2NT!'

      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ2NT call',
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ2NT call',
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB

      NUMINT = 3
      DX = (B-A)/ NUMINT
      DO 10  I = 1, NUMINT
          IF (I .EQ. 1)  THEN
             U(1) = A
             IF (IACTA .EQ. 0) THEN
               FU(1) = F(U(1))
             ELSE
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2.)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I))
             ELSE
               IB = I
               FB = F(B-DX/2.)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I))
          ENDIF
          CALL ADZ2AL(F,I)
   10     CONTINUE
       CALL TOT2LZ
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
          IF (ERS .GT. TARGET)  THEN
              NUMOLD = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*NUMOLD .GT. TARGET) CALL ADZ2PL(F,I,IER)
   40             CONTINUE
              IF (IER.EQ.0 .AND. NUMINT.NE.NUMOLD)  GOTO 30
              ENDIF
      ADZ2NT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      SUBROUTINE ADZ2PL (F, I, IER)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
      DATA TINY / 1.D-20 /

      IF (NUMINT .GE. MAXINT)  THEN
          IER = 1
          RETURN
          ENDIF
      NUMINT = NUMINT + 1
C                                                         New interval NUMINT
      IF (I .EQ. IB) IB = NUMINT
      U(NUMINT) = (U(I) + V(I)) / 2.
      V(NUMINT) = V(I)

      FU(NUMINT) = FW(I)
      FV(NUMINT) = FV(I)
C                                                             New interval I
       V(I) =  U(NUMINT)
      FV(I) = FU(NUMINT)
C                                                    Save old Result and Error
      OLDRES = RESULT(I)
      OLDERR = ERR(I)

      CALL ADZ2AL (F, I)
      CALL ADZ2AL (F, NUMINT)
C                                                               Update result
      DELRES = RESULT(I) + RESULT(NUMINT) - OLDRES
      RES = RES + DELRES
C                                  Good error estimate based on Simpson formula
      GODERR = ABS(DELRES)
C                                                             Update new global
      ERS = ERS + GODERR - OLDERR
C                                  Improve local error estimates proportionally
      SUMERR = ERR(I) + ERR(NUMINT)
      IF (SUMERR .GT. TINY) THEN
         FAC = GODERR / SUMERR
      ELSE
         FAC = 1.
      ENDIF

      ERR(I)      = ERR(I) * FAC
      ERR(NUMINT) = ERR(NUMINT) * FAC

      RETURN
C                        ****************************
      END

      SUBROUTINE ADZCAL (F,I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E15)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.

      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)

        CALL SGLINT (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL SGLINT (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
      ELSE
C                                                                   Closed endS
        FW(I) = F(W)
        TEM = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
C                                       Preliminary error Simpson - trapezoidal:
        ER  = DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.
      ENDIF

      RESULT(I) = TEM
      ERR   (I) = ABS (ER)

      RETURN
C                        ****************************
      END

      FUNCTION ADZINT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)

C     Adaptive integration routine which allows the integrand to be
C     indeterminant at the lower and/or the upper ends of integration.

C     Can self-adjust to any integrable singularity at the ends and compute
C     the closest approximant, hence achieve the required accuracy efficiently
C     (provided the switch(s) IACTA (IACTB) are set to 2).

C     Input switches for end-treatment:
C        IACTA = 0 :   Use closed lower-end algorithm
C                1 :   Open lower-end -- use open quadratic approximant
C                2 :   Open lower-end -- use adaptive singular approximant

C        IACTB = 0, 1, 2   (same as above, for the upper end)

C                Integral of F(X) from A to B, with error
C                less than ABS(AERR) + ABS(RERR*INTEGRAL)
C                Best estimate of error returned in ERREST.
C                Error code is IER:               0 :  o.k.
C                1 :  maximum calls to function reached before the
C                     error criteria are met;
C                2 :  IACTA out of range, set to 1;
C                3 :  IACTB out of range, set to 1.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      DATA SMLL / 1E-20 /

      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in ADZINT!'

      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call',
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call',
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB

      NUMINT = 3
      DX = (B-A)/ NUMINT
      DO 10  I = 1, NUMINT
          IF (I .EQ. 1)  THEN
             U(1) = A
             IF (IACTA .EQ. 0) THEN
               FU(1) = F(U(1))
             ELSE
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2.)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I))
             ELSE
               IB = I
               FB = F(B-DX/2.)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I))
          ENDIF
          CALL ADZCAL(F,I)
   10     CONTINUE
       CALL TOTALZ
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
          IF (ERS .GT. TARGET)  THEN
              NUMOLD = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*NUMOLD .GT. TARGET) CALL ADZSPL(F,I,IER)
   40         CONTINUE
              IF (IER.EQ.0 .AND. NUMINT.NE.NUMOLD)  GOTO 30
              ENDIF
      ADZINT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      SUBROUTINE ADZSPL (F, I, IER)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      DATA TINY / 1.D-20 /

      IF (NUMINT .GE. MAXINT)  THEN
          IER = 1
          RETURN
          ENDIF
      NUMINT = NUMINT + 1
C                                                         New interval NUMINT
      IF (I .EQ. IB) IB = NUMINT
      U(NUMINT) = (U(I) + V(I)) / 2.
      V(NUMINT) = V(I)

      FU(NUMINT) = FW(I)
      FV(NUMINT) = FV(I)
C                                                             New interval I
       V(I) =  U(NUMINT)
      FV(I) = FU(NUMINT)
C                                                    Save old Result and Error
      OLDRES = RESULT(I)
      OLDERR = ERR(I)

      CALL ADZCAL (F, I)
      CALL ADZCAL (F, NUMINT)
C                                                               Update result
      DELRES = RESULT(I) + RESULT(NUMINT) - OLDRES
      RES = RES + DELRES
C                                  Good error estimate based on Simpson formula
      GODERR = ABS(DELRES)
C                                                             Update new global
      ERS = ERS + GODERR - OLDERR
C                                  Improve local error estimates proportionally
      SUMERR = ERR(I) + ERR(NUMINT)
      IF (SUMERR .GT. TINY) THEN
         FAC = GODERR / SUMERR
      ELSE
         FAC = 1.
      ENDIF

      ERR(I)      = ERR(I) * FAC
      ERR(NUMINT) = ERR(NUMINT) * FAC

      RETURN
C                        ****************************
      END

      FUNCTION ASK(QUERY)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      LOGICAL ASK
      CHARACTER QUEND*11,  QUERY*(*), CH*1
C
      PARAMETER (QUEND= ' (Y OR N)? ' )
C
      CALL RTB (QUERY, LEN)
1     WRITE(NOUT, 90) QUERY(1:LEN), QUEND
      READ(NIN, 91) CH
      CALL UPCASE (CH)
      ASK = .FALSE.
      IF ( (CH.EQ.'Y') .OR. (CH.EQ.' ') ) THEN
         ASK = .TRUE.
         RETURN
      ELSE IF (CH.EQ.'N') THEN
         RETURN
      ELSE
         WRITE(NOUT, 92)
      ENDIF
      GOTO 1
 90     FORMAT ('$', 3A)
 91     FORMAT (A1)
 92     FORMAT (' BAD ANSWER--TRY AGAIN ')
C               *************************
      END
      FUNCTION BTA (X, Y)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      BTA = GAMA (X) * GAMA (Y) / GAMA (X+Y)
C
      RETURN
      END
C                       ****************************
C
      SUBROUTINE CHEPT (I)
C
C                               THIS ROUTINE PROVIDES CHECK-POINTS IN FORTRAN
C                               PROGRAMS FOR DEBUGGING PURPOSES
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      WRITE (NOUT, 900) I
  900 FORMAT ( 1X, 'CHECK-POINT  ', I2, '  FOUND HERE', /)
C
      RETURN
C               *************************
      END
C
        CHARACTER*(*) FUNCTION UCASE(A)
C               Converts A to all upper case.
        CHARACTER*(*) A
        UCASE = A
        CALL UC(UCASE)
        RETURN
        END
C                       ----------------------------
C
      SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION COVAR(NCVM,NCVM),LISTA(MFIT)
      DO 12 J=1,MA-1
        DO 11 I=J+1,MA
          COVAR(I,J)=0.
11      CONTINUE
12    CONTINUE
      DO 14 I=1,MFIT-1
        DO 13 J=I+1,MFIT
          IF(LISTA(J).GT.LISTA(I)) THEN
            COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
          ELSE
            COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
          ENDIF
13      CONTINUE
14    CONTINUE
      SWAP=COVAR(1,1)
      DO 15 J=1,MA
        COVAR(1,J)=COVAR(J,J)
        COVAR(J,J)=0.
15    CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
        COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16    CONTINUE
      DO 18 J=2,MA
        DO 17 I=1,J-1
          COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18    CONTINUE
      RETURN
C                        ****************************
      END

      FUNCTION EULER()
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (E=0.57721566)
      ENTRY EULERF()
      EULER = E
      RETURN
      END
C                       ****************************
C
      FUNCTION FINTR2 (FF,  X0, DX, NX, Y0, DY, NY,  XV, YV,  ERR, IR)

C     Two variable interpolation  --  Double Precision Version
C
c     Given array FF defined on evenly spaced lattice (0:Nx; 0:Ny), this
c     function routine calculates an interpolated value for the function at an
C     interrior point XV, YV.
c     The lowest value for the variable x is x0, the mesh-size is dx and the
c     array-size is 0:NX;  similarly for y:
C
C            TX         0                Tx              Nx
C            IZ         0  1  2   I-1  I         Nx-2    Nx
C                       |--|--| ... |--|--|--| ... |--|--|
C             X        X0                X               XM
C            XX                        0  1  2

C     It uses (MX-1)th ((MY-1)th) order polynomial fits to MX (MY) neighoring
C     points in the X (Y) direction.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MX=3, MY=3, MDX=MX/2, MDY=MY/2, AMDX=MDX, AMDY=MDY)
      DIMENSION FF (0:NX, 0:NY), XX(MX), YY(MY), ZZ(MX,MY)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DATA SMLL, HUGE / 1.D-6, 1.D30 /
      DATA XX, YY / 0.0, 1.0, 2.0, 0.0, 1.0, 2.0 /
      DATA  IW1, IW2, IW3, IW4, IW5, IW6, IW7 / 7 * 0 /
C                                                                     Initiation
      IR = 0
      X = XV
      Y = YV
      ERR = 0.
      ANX = NX
      ANY = NY
      FINTR2 = 0.
C                                                                          Tests
      IF (NX .LT. 1) THEN
         CALL WARNI(IW1, NWRT, 'Nx < 1, error in FINTR2.',
     >              'NX', NX, 1, 256, 1)
         IR = 1
         RETURN
      ELSE
         MNX = MIN (NX+1, MX)
      ENDIF
C
      IF (DX .LE. 0) THEN
         CALL WARNR(IW3, NWRT, 'DX < 0, error in FINTR2.',
     >              'DX', DX, D0, D1, 1)
         IR = 2
         RETURN
      ENDIF
C
      IF (NY .LT. 1) THEN
         CALL WARNI(IW2, NWRT, 'NY < 1, error in FINTR2.',
     >              'NY', NY, 1, 256, 1)
         IR = 3
         RETURN
      ELSE
         MNY = MIN(NY+1, MY)
      ENDIF
C
      IF (DY .LE. 0) THEN
         CALL WARNR(IW4, NWRT, 'DY < 0, error in FINTR2.',
     >              'DY', DY, D0, D1, 1)
         IR = 4
         RETURN
      ENDIF
C                                        Set up interpolation point and 3-array
      XM = X0 + DX * NX
      IF (X .LT. X0 .OR. X .GT. XM) THEN
        CALL WARNR(IW5,NWRT,'X out of range in FINTR2,
     >     Extrapolation used.','X',X,X0,XM,1)
      ENDIF
C
      TX = (X - X0) / DX
      IF (TX .LE. AMDX) THEN
        IX = 0
      ELSEIF (TX .GE. ANX-AMDX) THEN
        IX = NX - MX + 1
      ELSE
        IX = TX - MDX + 1
      ENDIF
      DDX = TX - IX
C                                  Set up interpolation point and 3-array for Y
      YM = Y0 + DY * NY
      ytiny = 1.e-3
c ytiny prevents complaints whose source is in small round off errors
c and not anything real
      IF (Y .LT. (Y0-ytiny) .OR. Y .GT. (YM+ytiny)) THEN
        CALL WARNR(IW6,NWRT,'Y out of range in FINTR2,
     >     Extrapolation used.','Y',Y,Y0,YM,1)
      ENDIF
C                                        Set up interpolation point and 3-array
      TY = (Y - Y0) / DY
      IF (TY .LE. AMDY) THEN
        IY = 0
      ELSEIF (TY .GE. ANY-AMDY) THEN
        IY = NY - MY + 1
      ELSE
        IY = TY - MDY + 1
      ENDIF
      DDY = TY - IY
C                                        POLIN2 is taken from "Numerical Recipe"
      DO 15 JX = 1, MNX
      DO 15 JY = 1, MNY
      ZZ (JX, JY) = FF (IX+JX-1, IY+JY-1)
   15 CONTINUE

      CALL POLIN2 (XX, YY, ZZ, MNX, MNY, DDX, DDY, TEM, ERR)

      FINTR2 = TEM
C
      RETURN
C
      END
C                       ****************************

	FUNCTION FINTRPO (FF,  X0, DX, NX,  XV,  ERR, IR)

C     Single variable interpolation  --  Double Precision Version
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MX = 3)
C
c     Given array FF defined on evenly spaced lattice (0:Nx), this function
c     routine calculates an interpolated value for the function at an
C     interrior point XV.
c     The lowest value for the variable x is x0, the mesh-size is dx and the
c     array-size is 0:NX;
C
C            TX         0                Tx              Nx
C            IZ         0  1  2   I-1  I         Nx-2    Nx
C                       |--|--| ... |--|--|--| ... |--|--|
C             X        X0                X               XM
C            XX                        0  1  2

C     It uses 2nd order polynomial fit from three neighoring points.

      DIMENSION FF (0:NX), XX(MX)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DATA SML, HUGE, XX / 1.D-7, 1.D30, 0., 1.0, 2.0 /
      DATA  IW1, IW2, IW3, IW4, IW5, IW6, IW7 / 7 * 0 /
C                                                                     Initiation
      IR = 0
      X = XV
      ERR = 0.
      ANX = NX
      FINTRPO = 0.
C                                                                          Tests
      IF (NX .LT. 1) THEN
         CALL WARNI(IW1, NWRT, 'Nx < 1, error in FINTRP.',
     >              'NX', NX, 1, 256, 1)
         IR = 1
         RETURN
      ELSE
         MNX = MIN(NX+1, MX)
      ENDIF
C
      IF (DX .LE. 0) THEN
         CALL WARNR(IW3, NWRT, 'DX < 0, error in FINTRP.',
     >              'DX', DX, D0, D1, 1)
         IR = 2
         RETURN
      ENDIF
C
      XM = X0 + DX * NX
      IF (X .LT. X0-SML .OR. X .GT. XM+SML) THEN
        CALL WARNR(IW5,NWRT,'X out of range in FINTRP,
     >      Extrapolation used.','X',X,X0,XM,1)
      IR = 3
      ENDIF
C                                        Set up interpolation point and 3-array
      TX = (X - X0) / DX
      IF (TX .LE. 1.) THEN
        IX = 0
      ELSEIF (TX .GE. ANX-1.) THEN
        IX = NX - 2
      ELSE
        IX = TX
      ENDIF
      DDX = TX - IX
C                                        POLINT is taken from "Numerical Recipe"
C      CALL POLINT (XX, FF(IX), MNX, DDX, TEM, ERR)
      CALL RATINTO (XX, FF(IX), MNX, DDX, TEM, ERR)

      FINTRPO = TEM
C
      RETURN
C
      END
C                       ****************************


      FUNCTION FINTRP (D,IV,IT,IK,FF,XI,X0,NDX,NX,XV,ERR,IR)

C     Single variable interpolation  --  Double Precision Version
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MX = 7, MXX = 1050)
C
c     Given array FF defined on evenly spaced lattice (0:Nx), this function
c     routine calculates an interpolated value for the function at an
C     interrior point XV.
c     The lowest value for the variable x is x0, the mesh-size is dx and the
c     array-size is 0:NX;
C
C            TX         0                Tx              Nx
C            IZ         0  1  2   I-1  I         Nx-2    Nx
C                       |--|--| ... |--|--|--| ... |--|--|
C             X        X0                X               XM
C            XX                        0  1  2

C     It uses 2nd order polynomial fit from three neighoring points.

      DIMENSION FF(0:MXX,0:MXX), XX(10), F1(10), XI(0:MXX),
     > XX1(10), F11(10)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DATA SML, HUGE/ 1.D-7, 1.D30/
      DATA  IW1, IW2, IW3, IW4, IW5, IW6, IW7 / 7 * 0 /
C                                                                     Initiation
      IR = 0
      X = XV
      ERR = 0.
      ANX = NX
      FINTRP = 0.
      MNX = MX

C                                        Set up interpolation point and 3-array

	IF (NDX.EQ.1 .AND. XI(IT).LE.D.AND.IT.LT.NX-7) THEN

	XX(1) = D
	XX(2) = XI(IT+1)
	XX(3) = XI(IT+2)
	XX(4) = XI(IT+3)
	XX(5) = XI(IT+4)
	XX(6) = XI(IT+5)
	XX(7) = XI(IT+6)
	F1(1) = FF(IK,IT)
	F1(2) = FF(IK,IT+1)
	F1(3) = FF(IK,IT+2)
	F1(4) = FF(IK,IT+3)
	F1(5) = FF(IK,IT+4)
	F1(6) = FF(IK,IT+5)
	F1(7) = FF(IK,IT+6)

	ELSEIF(NDX.EQ.1.AND.XI(IT).GT.D.AND.IT.LT.NX-7) THEN
      
	XX(1) = XI(IT)
	XX(2) = XI(IT+1)
	XX(3) = XI(IT+2)
	XX(4) = XI(IT+3)
	XX(5) = XI(IT+4)
	XX(6) = XI(IT+5)
	XX(7) = XI(IT+6)
	F1(1) = FF(IK,IT)
	F1(2) = FF(IK,IT+1)
	F1(3) = FF(IK,IT+2)
	F1(4) = FF(IK,IT+3)
	F1(5) = FF(IK,IT+4)
	F1(6) = FF(IK,IT+5)
	F1(7) = FF(IK,IT+6)

        ELSEIF(NDX.EQ.1 .AND. IT.GE.NX-7) THEN
      
	XX(1) = XI(NX-7)
	XX(2) = XI(NX-6)
	XX(3) = XI(NX-5)
	XX(4) = XI(NX-4)
	XX(5) = XI(NX-3)
	XX(6) = XI(NX-2)
	XX(7) = XI(NX-1)
	F1(1) = FF(IK,NX-7)
	F1(2) = FF(IK,NX-6)
	F1(3) = FF(IK,NX-5)
	F1(4) = FF(IK,NX-4)
	F1(5) = FF(IK,NX-3)
	F1(6) = FF(IK,NX-2)
	F1(7) = FF(IK,NX-1)


	ELSEIF (NDX.GT.1 .AND. IT.GE.1 .AND. IT.LT.NX-IV-7) THEN

	XX(1) = XI(IT+IV)
	XX(2) = XI(IT+IV+1)
	XX(3) = XI(IT+IV+2)
	XX(4) = XI(IT+IV+3)
	XX(5) = XI(IT+IV+4)
	XX(6) = XI(IT+IV+5)
	XX(7) = XI(IT+IV+6)
	F1(1) = FF(IT,IK)
	F1(2) = FF(IT+1,IK)
	F1(3) = FF(IT+2,IK)
	F1(4) = FF(IT+3,IK)
	F1(5) = FF(IT+4,IK)
	F1(6) = FF(IT+5,IK)
	F1(7) = FF(IT+6,IK)

	ELSEIF (NDX.GT.1 .AND. IT.GE.NX-IV-7) THEN
	
	XX(1) = XI(NX-7)
	XX(2) = XI(NX-6)
	XX(3) = XI(NX-5)
	XX(4) = XI(NX-4)
	XX(5) = XI(NX-3)
	XX(6) = XI(NX-2)
	XX(7) = XI(NX-1)
	F1(1) = FF(NX-IV-7,IK)
	F1(2) = FF(NX-IV-6,IK)
	F1(3) = FF(NX-IV-5,IK)
	F1(4) = FF(NX-IV-4,IK)
	F1(5) = FF(NX-IV-3,IK)
	F1(6) = FF(NX-IV-2,IK)
	F1(7) = FF(NX-IV-1,IK)

	ENDIF


C                                        POLINT is taken from "Numerical Recipe"
C      CALL POLINT (XX, F1, MNX, X, TEM, ERR)

	
      CALL RATINT (XX, F1, MNX, X, TEM, ERR)
	
         IF (X.LT.XX(1) .OR. X.GT.XX(MNX)) THEN

         GOTO 8

        ENDIF

        KO = 1
        IC = 1

	DO 21 KI = 1,MNX-1
	IF (X.GT.XX(KI).AND.X.LT.XX(KI+1)) THEN
	KO = KI
	ENDIF
 21	CONTINUE

        XMP = F1(KO)
        XMP1 = F1(KO+1)

 7      IF (XMP.LT.XMP1) THEN

	IF (TEM.LT.XMP.OR.TEM.GT.XMP1) THEN
c           print *,1,X,XX(KO),XX(KO+1)
           GOTO 6

        ELSE

           GOTO 8

	ENDIF

	ELSEIF (XMP.GT.XMP1) THEN

        IF (TEM.LT.XMP1.OR.TEM.GT.XMP) THEN

           GOTO 6

           ELSE

              GOTO 8

	ENDIF

        ENDIF

 6    XMI = (XX(KO) + XX(KO+1))/2.

      TEM1 = (F1(KO)+F1(KO+1))/2.

      IC = IC + 1
      K = 1

      DO 100 I = 1,MNX

         IF (I.EQ.KO) THEN

         XX1(K) = XX(I)
         F11(K) = F1(I)
         XX1(K+1) = XMI
         F11(K+1) = TEM1
         K = K + 2

         ELSE

         XX1(K) = XX(I)
         F11(K) = F1(I)

         K = K + 1

         ENDIF

 100     CONTINUE

         DO 101 P = 1,MNX

         XX(P) = XX1(P)
         F1(P) = F11(P)

 101     CONTINUE


         CALL RATINT (XX, F1, MNX, X, TEM, ERR)


	DO 22 KI = 1,MNX-1
	IF (X.GT.XX(KI).AND.X.LT.XX(KI+1)) THEN
	KO = KI
	ENDIF
 22     CONTINUE


        IF (IC.GT.5) THEN

           TEM = (XMP + XMP1)/2.

           GOTO 8

        ENDIF

         GOTO 7
C
 8      FINTRP = TEM

      RETURN
C
      END

      FUNCTION FINTRP01 (D,IT,IK,FF,XI,X0,NX,XV,ERR,IR)

C     Single variable interpolation  --  Double Precision Version
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MX = 7, MXX = 1050)
C
c     Given array FF defined on evenly spaced lattice (0:Nx), this function
c     routine calculates an interpolated value for the function at an
C     interrior point XV.
c     The lowest value for the variable x is x0, the mesh-size is dx and the
c     array-size is 0:NX;
C
C            TX         0                Tx              Nx
C            IZ         0  1  2   I-1  I         Nx-2    Nx
C                       |--|--| ... |--|--|--| ... |--|--|
C             X        X0                X               XM
C            XX                        0  1  2

C     It uses 2nd order polynomial fit from three neighoring points.

      DIMENSION FF(0:MXX,0:MXX), XX(10), F1(10), XI(0:MXX),
     > XX1(10), F11(10)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DATA SML, HUGE/ 1.D-7, 1.D30/
      DATA  IW1, IW2, IW3, IW4, IW5, IW6, IW7 / 7 * 0 /
C                                                                     Initiation
      IR = 0
      X = XV
      ERR = 0.
      ANX = NX
      FINTRP = 0.
      MNX = MX

C                                        Set up interpolation point and 3-array
      
        IF(XI(IT+6).LE.D) THEN

	XX(1) = XI(IT)
	XX(2) = XI(IT+1)
	XX(3) = XI(IT+2)
	XX(4) = XI(IT+3)
	XX(5) = XI(IT+4)
	XX(6) = XI(IT+5)
	XX(7) = XI(IT+6)
	F1(1) = FF(IK,IT)
	F1(2) = FF(IK,IT+1)
	F1(3) = FF(IK,IT+2)
	F1(4) = FF(IK,IT+3)
	F1(5) = FF(IK,IT+4)
	F1(6) = FF(IK,IT+5)
	F1(7) = FF(IK,IT+6)

        ELSE
      
	XX(1) = XI(IT-5)
	XX(2) = XI(IT-4)
	XX(3) = XI(IT-3)
	XX(4) = XI(IT-2)
	XX(5) = XI(IT-1)
	XX(6) = XI(IT)
	XX(7) = D
	F1(1) = FF(IK,IT-5)
	F1(2) = FF(IK,IT-4)
	F1(3) = FF(IK,IT-3)
	F1(4) = FF(IK,IT-2)
	F1(5) = FF(IK,IT-1)
	F1(6) = FF(IK,IT)
	F1(7) = FF(IK,IT+1)	

	ENDIF


C                                        POLINT is taken from "Numerical Recipe"
C      CALL POLINT (XX, F1, MNX, X, TEM, ERR)

	
      CALL RATINT (XX, F1, MNX, X, TEM, ERR)
	
         IF (X.LT.XX(1) .OR. X.GT.XX(MNX)) THEN

         GOTO 8

        ENDIF

        KO = 1
        IC = 1

	DO 21 KI = 1,MNX-1
	IF (X.GT.XX(KI).AND.X.LT.XX(KI+1)) THEN
	KO = KI
	ENDIF
 21	CONTINUE

        XMP = F1(KO)
        XMP1 = F1(KO+1)


 7      IF (XMP.LT.XMP1) THEN

	IF (TEM.LT.XMP.OR.TEM.GT.XMP1) THEN

           GOTO 6

        ELSE

           GOTO 8

	ENDIF

	ELSEIF (XMP.GT.XMP1) THEN

        IF (TEM.LT.XMP1.OR.TEM.GT.XMP) THEN

           GOTO 6

           ELSE

              GOTO 8

	ENDIF

        ENDIF

 6    XMI = (XX(KO) + XX(KO+1))/2.

      TEM1 = (F1(KO)+F1(KO+1))/2.

      IC = IC + 1
      K = 1

      DO 100 I = 1,MNX

         IF (I.EQ.KO) THEN

         XX1(K) = XX(I)
         F11(K) = F1(I)
         XX1(K+1) = XMI
         F11(K+1) = TEM1
         K = K + 2

         ELSE

         XX1(K) = XX(I)
         F11(K) = F1(I)

         K = K + 1

         ENDIF

 100     CONTINUE

         DO 101 P = 1,MNX

         XX(P) = XX1(P)
         F1(P) = F11(P)

 101     CONTINUE


         CALL RATINT (XX, F1, MNX, X, TEM, ERR)


	DO 22 KI = 1,MNX-1
	IF (X.GT.XX(KI).AND.X.LT.XX(KI+1)) THEN
	KO = KI
	ENDIF
 22     CONTINUE


        IF (IC.GT.5) THEN

           TEM = (XMP + XMP1)/2.

           GOTO 8

        ENDIF

         GOTO 7
C
 8      FINTRP01 = TEM


      RETURN
C
      END

C                       ****************************

	FUNCTION FINTRP1(D,IV,IT,IK,FFX,XI,X0,DX,NX,XV,ERR,IR)

C     Single variable interpolation  --  Double Precision Version
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MX = 7, MXX = 1050)
C
c     Given array FF defined on evenly spaced lattice (0:Nx), this function
c     routine calculates an interpolated value for the function at an
C     interrior point XV.
c     The lowest value for the variable x is x0, the mesh-size is dx and the
c     array-size is 0:NX;
C
C            TX         0                Tx              Nx
C            IZ         0  1  2   I-1  I         Nx-2    Nx
C                       |--|--| ... |--|--|--| ... |--|--|
C             X        X0                X               XM
C            XX                        0  1  2

C     It uses 2nd order polynomial fit from three neighoring points.

      DIMENSION FFX(MXX), XX(MX), F2(MX), XI(0:MXX),
     > XX1(10), F11(10)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DATA SML, HUGE/ 1.D-7, 1.D30/
      DATA  IW1, IW2, IW3, IW4, IW5, IW6, IW7 / 7 * 0 /
C                                                                     Initiation
      IR = 0
      X = XV
      ERR = 0.
      ANX = NX
      FINTRP1 = 0.
C                                                                          Tests
      IF (NX .LT. 1) THEN
         CALL WARNI(IW1, NWRT, 'Nx < 1, error in FINTRP.',
     >              'NX', NX, 1, 256, 1)
         IR = 1
         RETURN
      ELSE
         MNX = MIN(NX+1, MX)
      ENDIF
C
      IF (DX .LE. 0) THEN
         CALL WARNR(IW3, NWRT, 'DX < 0, error in FINTRP.',
     >              'DX', DX, D0, D1, 1)
         IR = 2
         RETURN
      ENDIF
C
      XM = -X0 + DX * NX


      IF (X .LT. X0-SML .OR. X .GT. XM+SML) THEN
        CALL WARNR(IW5,NWRT,'X out of range in FINTRP,
     >      Extrapolation used.','X',X,X0,XM,1)
      IR = 3
      ENDIF
C                                        Set up interpolation point and 3-array

      IF (IT.EQ.IV .AND. IK.EQ.4) THEN

	IT1 = IT - IV
	XX(1) = D
        XX(2) = XI(IV+1)
	XX(3) = XI(IV+2)
	XX(4) = XI(IV+3)
	XX(5) = XI(IV+4)
        XX(6) = XI(IV+5)
        XX(7) = XI(IV+6)
	F2(1) = FFX(IT1+1)
	F2(2) = FFX(IT1+2)
	F2(3) = FFX(IT1+3)
	F2(4) = FFX(IT1+4)
	F2(5) = FFX(IT1+5)
        F2(6) = FFX(IT1+6)
        F2(7) = FFX(IT1+7)

	ELSEIF (IK.EQ.0 .AND. IT.GE.IV+1 .AND. IT.LT.NX-5) THEN

	IT1 = IT - IV
	XX(1) = XI(IT)
	XX(2) = XI(IT+1)
	XX(3) = XI(IT+2)
	XX(4) = XI(IT+3)
	XX(5) = XI(IT+4)
        XX(6) = XI(IT+5)
        XX(7) = XI(IT+6)
	F2(1) = FFX(IT1)
	F2(2) = FFX(IT1+1)
	F2(3) = FFX(IT1+2)
	F2(4) = FFX(IT1+3)
	F2(5) = FFX(IT1+4)
        F2(6) = FFX(IT1+5)
        F2(7) = FFX(IT1+6)


	ELSEIF (IK.EQ.0 .AND. IT.GE.NX-5) THEN
	
	IT1 = IT - IV

        XX(1) = XI(IT-5)
        XX(2) = XI(IT-4)
	XX(3) = XI(IT-3)
	XX(4) = XI(IT-2)
	XX(5) = XI(IT-1)
	XX(6) = XI(IT)
	XX(7) = XI(IT+1)
        F2(1) = FFX(IT1-5)
        F2(2) = FFX(IT1-4)
	F2(3) = FFX(IT1-3)
	F2(4) = FFX(IT1-2)
	F2(5) = FFX(IT1-1)
	F2(6) = FFX(IT1)
	F2(7) = FFX(IT1+1)


	ELSEIF (IK.EQ.3 .AND. IT.LT.NX-5) THEN

	XX(1) = XI(IT)
	XX(2) = XI(IT+1)
	XX(3) = XI(IT+2)
	XX(4) = XI(IT+3)
	XX(5) = XI(IT+4)
        XX(6) = XI(IT+5)
        XX(7) = XI(IT+6)
	F2(1) = FFX(IT)
	F2(2) = FFX(IT+1)
	F2(3) = FFX(IT+2)
	F2(4) = FFX(IT+3)
	F2(5) = FFX(IT+4)
        F2(6) = FFX(IT+5)
        F2(7) = FFX(IT+6)

	ELSEIF (IK.EQ.3 .AND. IT.GE.NX-5) THEN
	
        XX(1) = XI(IT-5)
        XX(2) = XI(IT-4)
	XX(3) = XI(IT-3)
	XX(4) = XI(IT-2)
	XX(5) = XI(IT-1)
	XX(6) = XI(IT)
	XX(7) = XI(IT+1)
        F2(1) = FFX(IT-5)
        F2(2) = FFX(IT-4)
	F2(3) = FFX(IT-3)
	F2(4) = FFX(IT-2)
	F2(5) = FFX(IT-1)
	F2(6) = FFX(IT)
	F2(7) = FFX(IT+1)

	ELSEIF (IK.EQ.4 .AND. IT.GE.IV+1 .AND. IT.LT.NX-5) THEN

	IT1 = IT - IV
	XX(1) = XI(IT)
	XX(2) = XI(IT+1)
	XX(3) = XI(IT+2)
	XX(4) = XI(IT+3)
	XX(5) = XI(IT+4)
        XX(6) = XI(IT+5)
        XX(7) = XI(IT+6)
	F2(1) = FFX(IT1+1)
	F2(2) = FFX(IT1+2)
	F2(3) = FFX(IT1+3)
	F2(4) = FFX(IT1+4)
	F2(5) = FFX(IT1+5)
        F2(6) = FFX(IT1+6)
        F2(7) = FFX(IT1+7)

	ELSEIF (IK.EQ.4 .AND. IT.GE.NX-5) THEN
	
	IT1 = IT - IV
        XX(1) = XI(IT-5)
        XX(2) = XI(IT-4)
	XX(3) = XI(IT-3)
	XX(4) = XI(IT-2)
	XX(5) = XI(IT-1)
	XX(6) = XI(IT)
	XX(7) = XI(IT+1)
        F2(1) = FFX(IT1-4)
        F2(2) = FFX(IT1-3)
	F2(3) = FFX(IT1-2)
	F2(4) = FFX(IT1-1)
	F2(5) = FFX(IT1)
	F2(6) = FFX(IT1+1)
	F2(7) = FFX(IT1+2)

	ENDIF
	
      CALL RATINT (XX, F2, MNX, X, TEM, ERR)


       IF (X.LT.XX(1) .OR. X.GT.XX(MNX)) THEN

         GOTO 8

       ENDIF

      KO = 1
      IC = 1
	DO 21 KI = 1,MNX-1
	IF (X.GT.XX(KI).AND.X.LT.XX(KI+1)) THEN
	KO = KI
	ENDIF
 21	CONTINUE

        XMP = F2(KO)
        XMP1 = F2(KO+1)


 7      IF (XMP.LT.XMP1) THEN

	IF (TEM.LT.XMP.OR.TEM.GT.XMP1) THEN

           GOTO 6

        ELSE

           GOTO 8

	ENDIF

	ELSEIF (XMP.GT.XMP1) THEN

        IF (TEM.LT.XMP1.OR.TEM.GT.XMP) THEN

           GOTO 6

           ELSE

              GOTO 8

	ENDIF

        ENDIF

 6    XMI = (XX(KO) + XX(KO+1))/2.

      TEM1 = (F2(KO)+F2(KO+1))/2.
      IC = IC + 1
      K = 1

      DO 100 I = 1,MNX

         IF (I.EQ.KO) THEN

         XX1(K) = XX(I)
         F11(K) = F2(I)
         XX1(K+1) = XMI
         F11(K+1) = TEM1
         K = K + 2

         ELSE

         XX1(K) = XX(I)
         F11(K) = F2(I)

         K = K + 1

         ENDIF

 100     CONTINUE

         DO 101 P = 1,MNX

         XX(P) = XX1(P)
         F2(P) = F11(P)

 101     CONTINUE


         CALL RATINT (XX, F2, MNX, X, TEM, ERR)


	DO 22 KI = 1,MNX-1
	IF (X.GT.XX(KI).AND.X.LT.XX(KI+1)) THEN
	KO = KI
	ENDIF
 22     CONTINUE

         IF (IC.GT.5) THEN

           TEM = (XMP + XMP1)/2.

           GOTO 8

        ENDIF

         GOTO 7
C
 8      FINTRP1 = TEM

C
      RETURN
C
      END
C                       ****************************

	FUNCTION FINTRP2 (ID,D,IV,IT,IK,FF,XI,X0,DX,NX,XV,ERR,IR)

C     Single variable interpolation  --  Double Precision Version
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MX = 7, MXX = 1050)
C
c     Given array FF defined on evenly spaced lattice (0:Nx), this function
c     routine calculates an interpolated value for the function at an
C     interrior point XV.
c     The lowest value for the variable x is x0, the mesh-size is dx and the
c     array-size is 0:NX;
C
C            TX         0                Tx              Nx
C            IZ         0  1  2   I-1  I         Nx-2    Nx
C                       |--|--| ... |--|--|--| ... |--|--|
C             X        X0                X               XM
C            XX                        0  1  2

C     It uses 2nd order polynomial fit from three neighoring points.

      DIMENSION FF(0:MXX,0:MXX), XX(10), F1(10), XI(0:MXX),
     > XX1(10),F11(10)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DATA SML, HUGE/ 1.D-7, 1.D30/
      DATA  IW1, IW2, IW3, IW4, IW5, IW6, IW7 / 7 * 0 /
C                                                                     Initiation
      IR = 0
      X = XV
      ERR = 0.
      ANX = NX
      FINTRP2 = 0.0
C                                                                          Tests
      IF (NX .LT. 1) THEN
         CALL WARNI(IW1, NWRT, 'Nx < 1, error in FINTRP.',
     >              'NX', NX, 1, 256, 1)
         IR = 1
         RETURN
      ELSE
         MNX = MIN(NX+1, MX)
      ENDIF
C
      IF (DX .LE. 0) THEN
         CALL WARNR(IW3, NWRT, 'DX < 0, error in FINTRP.',
     >              'DX', DX, D0, D1, 1)
         IR = 2
         RETURN
      ENDIF
C
      XM = D


      IF (X .LT. X0-SML .OR. X .GT. XM+SML) THEN
        CALL WARNR(IW5,NWRT,'X out of range in FINTRP,
     >      Extrapolation used.','X',X,X0,XM,1)
      IR = 3
      ENDIF
C                                        Set up interpolation point and 3-array

	M1 = MNX

	IU = IT - IK

        IF (ID.EQ.1 .AND. IT.LT.8) THEN

        XX(1) = XI(1)
	XX(2) = XI(2)
	XX(3) = XI(3)
	XX(4) = XI(4)
	XX(5) = XI(5)
	XX(6) = XI(6)
	XX(7) = XI(7)
	F1(1) = FF(IK,1)
	F1(2) = FF(IK,2)
	F1(3) = FF(IK,3)
	F1(4) = FF(IK,4)
	F1(5) = FF(IK,5)
	F1(6) = FF(IK,6)
	F1(7) = FF(IK,7)
           
        ELSEIF(ID.EQ.1.AND.IT.GE.8.AND.IU.LT.-6) THEN

	XX(1) = XI(IT-3)
	XX(2) = XI(IT-2)
	XX(3) = XI(IT-1)
	XX(4) = XI(IT)
	XX(5) = XI(IT+1)
	XX(6) = XI(IT+2)
	XX(7) = XI(IT+3)
	XX(8) = XI(IT+4)
	F1(1) = FF(IK,IT-3)
	F1(2) = FF(IK,IT-2)
	F1(3) = FF(IK,IT-1)
	F1(4) = FF(IK,IT)
	F1(5) = FF(IK,IT+1)
	F1(6) = FF(IK,IT+2)
	F1(7) = FF(IK,IT+3)
	F1(8) = FF(IK,IT+4)
	M1 = 8

	ELSEIF (ID.EQ.1.AND.IT.GE.8.AND.IU.GE.-6) THEN

	XX(1) = XI(IK-6)
	XX(2) = XI(IK-5)
	XX(3) = XI(IK-4)
	XX(4) = XI(IK-3)
	XX(5) = XI(IK-2)
	XX(6) = XI(IK-1)
	XX(7) = XI(IK)
	F1(1) = FF(IK,IK-6)
	F1(2) = FF(IK,IK-5)
	F1(3) = FF(IK,IK-4)
	F1(4) = FF(IK,IK-3)
	F1(5) = FF(IK,IK-2)
	F1(6) = FF(IK,IK-1)
	F1(7) = FF(IK,IK)


	ELSEIF(ID.EQ.2 .AND. IT.LT.IV-5) THEN
      
	XX(1) = XI(IT)
	XX(2) = XI(IT+1)
	XX(3) = XI(IT+2)
	XX(4) = XI(IT+3)
	XX(5) = XI(IT+4)
	XX(6) = XI(IT+5)
	XX(7) = XI(IT+6)
	F1(1) = FF(IK,IT)
	F1(2) = FF(IK,IT+1)
	F1(3) = FF(IK,IT+2)
	F1(4) = FF(IK,IT+3)
	F1(5) = FF(IK,IT+4)
	F1(6) = FF(IK,IT+5)
	F1(7) = FF(IK,IT+6)

	ELSEIF (ID.EQ.2 .AND. IT.GE.IV-5) THEN

	XX(1) = XI(IV-6)
	XX(2) = XI(IV-5)
	XX(3) = XI(IV-4)
	XX(4) = XI(IV-3)
	XX(5) = XI(IV-2)
	XX(6) = XI(IV-1)
	XX(7) = XI(IV)
	F1(1) = FF(IK,IV-6)
	F1(2) = FF(IK,IV-5)
	F1(3) = FF(IK,IV-4)
	F1(4) = FF(IK,IV-3)
	F1(5) = FF(IK,IV-2)
	F1(6) = FF(IK,IV-1)
	F1(7) = FF(IK,IV)

	ENDIF

	
      CALL RATINT (XX, F1, M1, X, TEM, ERR)
	      

       IF (X.LT.XX(1) .OR. X.GT.XX(M1)) THEN

         GOTO 8

       ENDIF
      

       KO = 1
       IC = 1
	DO 21 KI = 1,M1-1
	IF (X.GT.XX(KI).AND.X.LT.XX(KI+1)) THEN
	KO = KI
	ENDIF
 21	CONTINUE

        XMP = F1(KO)
        XMP1 = F1(KO+1)

 7      IF (XMP.LT.XMP1) THEN

	IF (TEM.LT.XMP.OR.TEM.GT.XMP1) THEN

           GOTO 6

        ELSE

           GOTO 8

	ENDIF

	ELSEIF (XMP.GT.XMP1) THEN

        IF (TEM.LT.XMP1.OR.TEM.GT.XMP) THEN

           GOTO 6

           ELSE

              GOTO 8

	ENDIF

        ENDIF

 6    XMI = (XX(KO) + XX(KO+1))/2.

      TEM1 = (F1(KO)+F1(KO+1))/2.
      IC = IC + 1
      K = 1

      DO 100 I = 1,M1

         IF (I.EQ.KO) THEN

         XX1(K) = XX(I)
         F11(K) = F1(I)
         XX1(K+1) = XMI
         F11(K+1) = TEM1
         K = K + 2

         ELSE

         XX1(K) = XX(I)
         F11(K) = F1(I)

         K = K + 1

         ENDIF

 100     CONTINUE

         DO 101 P = 1,M1

         XX(P) = XX1(P)
         F1(P) = F11(P)

 101     CONTINUE


         CALL RATINT (XX, F1, M1, X, TEM, ERR)


         IF (IC.GT.5) THEN

           TEM = (XMP + XMP1)/2.

           GOTO 8

        ENDIF

        DO 22 KI = 1,M1-1
	IF (X.GT.XX(KI).AND.X.LT.XX(KI+1)) THEN
	KO = KI
	ENDIF
 22     CONTINUE

         GOTO 7
C
 8      FINTRP2 = TEM

C
      RETURN
C
      END
C                       ****************************

	FUNCTION FINTRP3 (ID,D,IV,IT,IK,FF,XI,X0,DX,NX,XV,ERR,IR)

C     Single variable interpolation  --  Double Precision Version
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MX = 7, MXX = 1050)
C
c     Given array FF defined on evenly spaced lattice (0:Nx), this function
c     routine calculates an interpolated value for the function at an
C     interrior point XV.
c     The lowest value for the variable x is x0, the mesh-size is dx and the
c     array-size is 0:NX;
C
C            TX         0                Tx              Nx
C            IZ         0  1  2   I-1  I         Nx-2    Nx
C                       |--|--| ... |--|--|--| ... |--|--|
C             X        X0                X               XM
C            XX                        0  1  2

C     It uses 2nd order polynomial fit from three neighoring points.

      DIMENSION FF(MXX), XX(10), F1(10), XI(0:MXX),
     > XX1(10), F11(10)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DATA SML, HUGE/ 1.D-7, 1.D30/
      DATA  IW1, IW2, IW3, IW4, IW5, IW6, IW7 / 7 * 0 /
C                                                                     Initiation
      IR = 0
      X = XV
      ERR = 0.
      ANX = NX
      FINTRP3 = 0.0
C                                                                          Tests
      IF (NX .LT. 1) THEN
         CALL WARNI(IW1, NWRT, 'Nx < 1, error in FINTRP.',
     >              'NX', NX, 1, 256, 1)
         IR = 1
         RETURN
      ELSE
         MNX = MIN(NX+1, MX)
      ENDIF
C
      IF (DX .LE. 0) THEN
         CALL WARNR(IW3, NWRT, 'DX < 0, error in FINTRP.',
     >              'DX', DX, D0, D1, 1)
         IR = 2
         RETURN
      ENDIF
C
      XM = D


      IF (X .LT. X0-SML .OR. X .GT. XM+SML) THEN
        CALL WARNR(IW5,NWRT,'X out of range in FINTRP,
     >      Extrapolation used.','X',X,X0,XM,1)
      IR = 3
      ENDIF
C                                        Set up interpolation point and 3-array
	
        M1 = MNX

	IU = IT - IK

	IF (ID.EQ.1 .AND. IT.LT.6) THEN

        XX(1) = XI(0)
	XX(2) = XI(1)
	XX(3) = XI(2)
	XX(4) = XI(3)
	XX(5) = XI(4)
	XX(6) = XI(5)
	XX(7) = XI(6)
	F1(1) = FF(1)
	F1(2) = FF(2)
	F1(3) = FF(3)
	F1(4) = FF(4)
	F1(5) = FF(5)
	F1(6) = FF(6)
	F1(7) = FF(7)
           
        ELSEIF(ID.EQ.1.AND.IT.GE.6.AND.IU.LT.-6) THEN

	XX(1) = XI(IT-3)
	XX(2) = XI(IT-2)
	XX(3) = XI(IT-1)
	XX(4) = XI(IT)
	XX(5) = XI(IT+1)
	XX(6) = XI(IT+2)
	XX(7) = XI(IT+3)
	XX(8) = XI(IT+4)
	F1(1) = FF(IT-2)
	F1(2) = FF(IT-1)
	F1(3) = FF(IT)
	F1(4) = FF(IT+1)
	F1(5) = FF(IT+2)
	F1(6) = FF(IT+3)
	F1(7) = FF(IT+4)
	F1(8) = FF(IT+5)
	M1 = 8

	ELSEIF (ID.EQ.1.AND.IT.GE.6.AND.IU.GE.-6) THEN

	XX(1) = XI(IK-6)
	XX(2) = XI(IK-5)
	XX(3) = XI(IK-4)
	XX(4) = XI(IK-3)
	XX(5) = XI(IK-2)
	XX(6) = XI(IK-1)
	XX(7) = XI(IK)
	F1(1) = FF(IK-5)
	F1(2) = FF(IK-4)
	F1(3) = FF(IK-3)
	F1(4) = FF(IK-2)
	F1(5) = FF(IK-1)
	F1(6) = FF(IK)
	F1(7) = FF(IK+1)

        ELSEIF(ID.EQ.2 .AND. IU.LT.6) THEN
      
	XX(1) = XI(IK)
	XX(2) = XI(IK+1)
	XX(3) = XI(IK+2)
	XX(4) = XI(IK+3)
	XX(5) = XI(IK+4)
	XX(6) = XI(IK+5)
	XX(7) = XI(IK+6)
	F1(1) = FF(IK)
	F1(2) = FF(IK+1)
	F1(3) = FF(IK+2)
	F1(4) = FF(IK+3)
	F1(5) = FF(IK+4)
	F1(6) = FF(IK+5)
	F1(7) = FF(IK+6)

	ELSEIF(ID.EQ.2.AND.IU.GE.6.AND.IT.LT.IV-5) THEN
      
	XX(1) = XI(IT-3)
	XX(2) = XI(IT-2)
	XX(3) = XI(IT-1)
	XX(4) = XI(IT)
	XX(5) = XI(IT+1)
	XX(6) = XI(IT+2)
	XX(7) = XI(IT+3)
	F1(1) = FF(IT-3)
	F1(2) = FF(IT-2)
	F1(3) = FF(IT-1)
	F1(4) = FF(IT)
	F1(5) = FF(IT+1)
	F1(6) = FF(IT+2)
	F1(7) = FF(IT+3)

	ELSEIF (ID.EQ.2 .AND. IT.GE.IV-5) THEN

	XX(1) = XI(IV-5)
	XX(2) = XI(IV-4)
	XX(3) = XI(IV-3)
	XX(4) = XI(IV-2)
	XX(5) = XI(IV-1)
	XX(6) = XI(IV)
	XX(7) = D
	F1(1) = FF(IV-5)
	F1(2) = FF(IV-4)
	F1(3) = FF(IV-3)
	F1(4) = FF(IV-2)
	F1(5) = FF(IV-1)
	F1(6) = FF(IV)
	F1(7) = FF(IV+1)

	ENDIF


C                                        POLINT is taken from "Numerical Recipe"
C      CALL POLINT (XX, F1, MNX, X, TEM, ERR)

	
      CALL RATINT (XX, F1, M1, X, TEM, ERR)

      IF (X.LT.XX(1) .OR. X.GT.XX(M1)) THEN

         GOTO 8

      ENDIF

      KO = 1
      IC = 1
	DO 21 KI = 1,M1-1
	IF (X.GT.XX(KI).AND.X.LT.XX(KI+1)) THEN
	KO = KI
	ENDIF
 21	CONTINUE

        XMP = F1(KO)
        XMP1 = F1(KO+1)

 7      IF (XMP.LT.XMP1) THEN

	IF (TEM.LT.XMP.OR.TEM.GT.XMP1) THEN

           GOTO 6

        ELSE

           GOTO 8

	ENDIF

	ELSEIF (XMP.GT.XMP1) THEN

        IF (TEM.LT.XMP1.OR.TEM.GT.XMP) THEN

           GOTO 6

           ELSE

              GOTO 8

	ENDIF

        ENDIF

 6    XMI = (XX(KO) + XX(KO+1))/2.

      TEM1 = (F1(KO)+F1(KO+1))/2.
      IC = IC + 1 
      K = 1

      DO 100 I = 1,M1

         IF (I.EQ.KO) THEN

         XX1(K) = XX(I)
         F11(K) = F1(I)
         XX1(K+1) = XMI
         F11(K+1) = TEM1
         K = K + 2

         ELSE

         XX1(K) = XX(I)
         F11(K) = F1(I)

         K = K + 1

         ENDIF

 100     CONTINUE

         DO 101 P = 1,M1

         XX(P) = XX1(P)
         F1(P) = F11(P)

 101     CONTINUE


         CALL RATINT (XX, F1, M1, X, TEM, ERR)


	DO 22 KI = 1,M1-1
	IF (X.GT.XX(KI).AND.X.LT.XX(KI+1)) THEN
	KO = KI
	ENDIF
 22     CONTINUE

         IF (IC.GT.5) THEN

           TEM = (XMP + XMP1)/2.

           GOTO 8

        ENDIF

         GOTO 7
C
 8      FINTRP3 = TEM

C
      RETURN
C
      END
C                       ****************************

      FUNCTION FNTR2P (FF, X0,DX,NX, P0,DP,NP, X, P)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
c     Given function ff defined on the two-dimensional array, this function
c     routine calculates an interpolated value for the function at any
c     given point.
c     The lowest value for the variable x is x0, the mesh-size is dx and the
c     array-size is 0:NX;  the minimum value for the variable p is p0, the
c     increment it dp, and the array-size is 0:NP.
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DIMENSION FF(0:NX, 0:NP)
C
      IF (X.EQ.X0+DX*NX.AND.P.EQ.P0+DP*NP) THEN
                        FNTR2P = FF(NX,NP)
                        RETURN
      END IF
C
      TX = (X - X0) / DX
      IX = TX
      DDX = TX - IX
      IF (IX .LT. 0 .OR. IX .GT. NX) GOTO 90
C
      TP = (P - P0) / DP
      IP = TP
      DDP = TP - IP
      IF (IP .LT. 0 .OR. IP .GT. NP) GOTO 90
C
      FNTR2P = (1.-DDX) * (1.-DDP) * FF(IX  , IP  )
     >         +    DDX  * (1.-DDP) * FF(IX+1, IP  )
     >         +(1.-DDX) *     DDP  * FF(IX  , IP+1)
     >         +    DDX  *     DDP  * FF(IX+1, IP+1)
C
      RETURN
C
   90 WRITE (NOUT, 990) X, P
  990 FORMAT (' value(s) of x, and/or p out of range in interp'/
     > / '  x =', F10.3, ' p =', F10.3 / )
C
      STOP
C                       *******************************
      END

      FUNCTION GAMA (X)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      Z = X
      CALL LOGGAM (Z, U)
      GAMA = EXP (U)
C
      RETURN
      END
C                       ****************************
C
      SUBROUTINE GAMMA(Z,G)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      X = Z
      CALL LOGGAM(X,U)
      G = EXP(U)
C
      RETURN
      END
C                       ****************************
C
      FUNCTION GAMMLN(XX)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
C                        ****************************
      END

      FUNCTION GAMMQ(A,X)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IF(X.LT.0..OR.A.LE.0.)PAUSE
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMQ,A,X,GLN)
      ENDIF
      RETURN
C                        ****************************
      END

      FUNCTION GAUSSO(F,XL,XR,AERR,RERR,ERR,IRT)
C                                                       FUNCTION XINT(XL,F,XR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
        DIMENSION XLIMS(100), R(93), W(93)
        INTEGER PTR(4),NORD(4),NIN,NOUT,NWRT

        external f

        COMMON/IOUNIT/NIN,NOUT,NWRT
        DATA PTR,NORD/4,10,22,46,  6,12,24,48/
        DATA R/.2386191860,.6612093865,.9324695142,
     1 .1252334085,.3678314990,.5873179543,.7699026742,.9041172563,
     1 .9815606342,.0640568929,.1911188675,.3150426797,.4337935076,
     1 .5454214714,.6480936519,.7401241916,.8200019860,.8864155270,
     1 .9382745520,.9747285560,.9951872200,.0323801710,.0970046992,
     1 .1612223561,.2247637903,.2873624873,.3487558863,.4086864820,
     1 .4669029048,.5231609747,.5772247261,.6288673968,.6778723796,
     1 .7240341309,.7671590325,.8070662040,.8435882616,.8765720203,
     1 .9058791367,.9313866907,.9529877032,.9705915925,.9841245837,
     1 .9935301723,.9987710073,.0162767488,.0488129851,.0812974955,
     1 .1136958501,.1459737146,.1780968824,.2100313105,.2417431561,
     1 .2731988126,.3043649444,.3352085229,.3656968614,.3957976498,
     1 .4254789884,.4547094222,.4834579739,.5116941772,.5393881083,
     1 .5665104186,.5930323648,.6189258401,.6441634037,.6687183100,
     1 .6925645366,.7156768123,.7380306437,.7596023411,.7803690438,
     1 .8003087441,.8194003107,.8376235112,.8549590334,.8713885059,
     1 .8868945174,.9014606353,.9150714231,.9277124567,.9393703398,
     1 .9500327178,.9596882914,.9683268285,.9759391746,.9825172636,
     1 .9880541263,.9925439003,.9959818430,.9983643759,.9996895039/
        DATA W/.4679139346,.3607615730,.1713244924,
     1 .2491470458,.2334925365,.2031674267,.1600783285,.1069393260,
     1 .0471753364,.1279381953,.1258374563,.1216704729,.1155056681,
     1 .1074442701,.0976186521,.0861901615,.0733464814,.0592985849,
     1 .0442774388,.0285313886,.0123412298,.0647376968,.0644661644,
     1 .0639242386,.0631141923,.0620394232,.0607044392,.0591148397,
     1 .0572772921,.0551995037,.0528901894,.0503590356,.0476166585,
     1 .0446745609,.0415450829,.0382413511,.0347772226,.0311672278,
     1 .0274265097,.0235707608,.0196161605,.0155793157,.0114772346,
     1 .0073275539,.0031533461,.0325506145,.0325161187,.0324471637,
     1 .0323438226,.0322062048,.0320344562,.0318287589,.0315893308,
     1 .0313164256,.0310103326,.0306713761,.0302999154,.0298963441,
     1 .0294610900,.0289946142,.0284974111,.0279700076,.0274129627,
     1 .0268268667,.0262123407,.0255700360,.0249006332,.0242048418,
     1 .0234833991,.0227370697,.0219666444,.0211729399,.0203567972,
     1 .0195190811,.0186606796,.0177825023,.0168854799,.0159705629,
     1 .0150387210,.0140909418,.0131282296,.0121516047,.0111621020,
     1 .0101607705,.0091486712,.0081268769,.0070964708,.0060585455,
     1 .0050142027,.0039645543,.0029107318,.0018539608,.0007967921/
        DATA TOLABS,TOLREL,NMAX/1.E-35,5.E-4,100/
C
C
        TOLABS=AERR
        TOLREL=RERR

        GAUSSO=0.
        NLIMS=2
        XLIMS(1)=XL
        XLIMS(2)=XR
C
10      AA=(XLIMS(NLIMS)-XLIMS(NLIMS-1))/2D0
        BB=(XLIMS(NLIMS)+XLIMS(NLIMS-1))/2D0
        TVAL=0.
        DO 15 I=1,3
15      TVAL=TVAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
        TVAL=TVAL*AA
        DO 25 J=1,4
        VAL=0.
        DO 20 I=PTR(J),PTR(J)-1+NORD(J)
20      VAL=VAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
        VAL=VAL*AA
        TOL=MAX(TOLABS,TOLREL*ABS(VAL))
        IF (ABS(TVAL-VAL).LT.TOL) THEN
                GAUSSO=GAUSSO+VAL
                NLIMS=NLIMS-2
                IF (NLIMS.NE.0) GO TO 10
                RETURN
                END IF
25      TVAL=VAL
        IF (NMAX.EQ.2) THEN
                GAUSSO=VAL
                RETURN
                END IF
        IF (NLIMS.GT.(NMAX-2)) THEN
                WRITE(NOUT,50) GAUSSO,NMAX,BB-AA,BB+AA
                RETURN
                END IF
        XLIMS(NLIMS+1)=BB
        XLIMS(NLIMS+2)=BB+AA
        XLIMS(NLIMS)=BB
        NLIMS=NLIMS+2
        GO TO 10
C
50      FORMAT (' GAUSSO FAILS, GAUSS,NMAX,XL,XR=',G15.7,I5,2G15.7)
C                        ****************************
        END

      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
C                        ****************************
      END

      FUNCTION GAUSS(IT,IT1,F,XL,XR,AERR,RERR,ERR,IRT)
C                                                       FUNCTION XINT(XL,F,XR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
        DIMENSION XLIMS(100), R(93), W(93)
        INTEGER PTR(4),NORD(4),NIN,NOUT,NWRT

        external f

        COMMON/IOUNIT/NIN,NOUT,NWRT
        DATA PTR,NORD/4,10,22,46,  6,12,24,48/
        DATA R/.2386191860,.6612093865,.9324695142,
     1 .1252334085,.3678314990,.5873179543,.7699026742,.9041172563,
     1 .9815606342,.0640568929,.1911188675,.3150426797,.4337935076,
     1 .5454214714,.6480936519,.7401241916,.8200019860,.8864155270,
     1 .9382745520,.9747285560,.9951872200,.0323801710,.0970046992,
     1 .1612223561,.2247637903,.2873624873,.3487558863,.4086864820,
     1 .4669029048,.5231609747,.5772247261,.6288673968,.6778723796,
     1 .7240341309,.7671590325,.8070662040,.8435882616,.8765720203,
     1 .9058791367,.9313866907,.9529877032,.9705915925,.9841245837,
     1 .9935301723,.9987710073,.0162767488,.0488129851,.0812974955,
     1 .1136958501,.1459737146,.1780968824,.2100313105,.2417431561,
     1 .2731988126,.3043649444,.3352085229,.3656968614,.3957976498,
     1 .4254789884,.4547094222,.4834579739,.5116941772,.5393881083,
     1 .5665104186,.5930323648,.6189258401,.6441634037,.6687183100,
     1 .6925645366,.7156768123,.7380306437,.7596023411,.7803690438,
     1 .8003087441,.8194003107,.8376235112,.8549590334,.8713885059,
     1 .8868945174,.9014606353,.9150714231,.9277124567,.9393703398,
     1 .9500327178,.9596882914,.9683268285,.9759391746,.9825172636,
     1 .9880541263,.9925439003,.9959818430,.9983643759,.9996895039/
        DATA W/.4679139346,.3607615730,.1713244924,
     1 .2491470458,.2334925365,.2031674267,.1600783285,.1069393260,
     1 .0471753364,.1279381953,.1258374563,.1216704729,.1155056681,
     1 .1074442701,.0976186521,.0861901615,.0733464814,.0592985849,
     1 .0442774388,.0285313886,.0123412298,.0647376968,.0644661644,
     1 .0639242386,.0631141923,.0620394232,.0607044392,.0591148397,
     1 .0572772921,.0551995037,.0528901894,.0503590356,.0476166585,
     1 .0446745609,.0415450829,.0382413511,.0347772226,.0311672278,
     1 .0274265097,.0235707608,.0196161605,.0155793157,.0114772346,
     1 .0073275539,.0031533461,.0325506145,.0325161187,.0324471637,
     1 .0323438226,.0322062048,.0320344562,.0318287589,.0315893308,
     1 .0313164256,.0310103326,.0306713761,.0302999154,.0298963441,
     1 .0294610900,.0289946142,.0284974111,.0279700076,.0274129627,
     1 .0268268667,.0262123407,.0255700360,.0249006332,.0242048418,
     1 .0234833991,.0227370697,.0219666444,.0211729399,.0203567972,
     1 .0195190811,.0186606796,.0177825023,.0168854799,.0159705629,
     1 .0150387210,.0140909418,.0131282296,.0121516047,.0111621020,
     1 .0101607705,.0091486712,.0081268769,.0070964708,.0060585455,
     1 .0050142027,.0039645543,.0029107318,.0018539608,.0007967921/
        DATA TOLABS,TOLREL,NMAX/1.E-35,5.E-4,100/
C
C
        TOLABS=AERR
        TOLREL=RERR

        GAUSS=0.
        NLIMS=2
        XLIMS(1)=XL
        XLIMS(2)=XR
C
10      AA=(XLIMS(NLIMS)-XLIMS(NLIMS-1))/2D0
        BB=(XLIMS(NLIMS)+XLIMS(NLIMS-1))/2D0
        TVAL=0.
        DO 15 I=1,3
15      TVAL=TVAL+W(I)*(F(IT,IT1,BB+AA*R(I))+F(IT,IT1,BB-AA*R(I)))
        TVAL=TVAL*AA
        DO 25 J=1,4
        VAL=0.
        DO 20 I=PTR(J),PTR(J)-1+NORD(J)
20      VAL=VAL+W(I)*(F(IT,IT1,BB+AA*R(I))+F(IT,IT1,BB-AA*R(I)))
        VAL=VAL*AA
        TOL=MAX(TOLABS,TOLREL*ABS(VAL))
        IF (ABS(TVAL-VAL).LT.TOL) THEN
                GAUSS=GAUSS+VAL
                NLIMS=NLIMS-2
                IF (NLIMS.NE.0) GO TO 10
                RETURN
                END IF
25      TVAL=VAL
        IF (NMAX.EQ.2) THEN
                GAUSS=VAL
                RETURN
                END IF
        IF (NLIMS.GT.(NMAX-2)) THEN
                WRITE(NOUT,50) GAUSS,NMAX,BB-AA,BB+AA
                RETURN
                END IF
        XLIMS(NLIMS+1)=BB
        XLIMS(NLIMS+2)=BB+AA
        XLIMS(NLIMS)=BB
        NLIMS=NLIMS+2
        GO TO 10
C
50      FORMAT (' GAUSS FAILS, GAUSS,NMAX,XL,XR=',G15.7,I5,2G15.7)
C                        ****************************
        END

 
C                      ================================
C                              LIBMSC.FOR
C                       ----------------------------

      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*LOG(X)-GLN)*G
      RETURN
C                        ****************************
      END

      SUBROUTINE GLFIT
     > (X,Y,SIG,NDATA,A,MA,LISTA,MFIT,COVAR,NCVM,CHISQ,FUNCS)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                             Fit data to linear combination of given functions
C                                            From LFIT of Numerical Recipes
      PARAMETER (MMAX=50)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),LISTA(MFIT),
     *    COVAR(NCVM,NCVM),BETA(MMAX),AFUNC(MMAX)

      external funcs

      KK=MFIT+1
      DO 12 J=1,MA
        IHIT=0
        DO 11 K=1,MFIT
          IF (LISTA(K).EQ.J) IHIT=IHIT+1
11      CONTINUE
        IF (IHIT.EQ.0) THEN
          LISTA(KK)=J
          KK=KK+1
        ELSE IF (IHIT.GT.1) THEN
          PAUSE 'Improper set in LISTA in GLFIT'
        ENDIF
12    CONTINUE
      IF (KK.NE.(MA+1)) PAUSE 'Improper set in LISTA in GLFIT'
      DO 14 J=1,MFIT
        DO 13 K=1,MFIT
          COVAR(J,K)=0.
13      CONTINUE
        BETA(J)=0.
14    CONTINUE
      DO 18 I=1,NDATA
        CALL FUNCS(X(I),AFUNC,MA)
        YM=Y(I)
        IF(MFIT.LT.MA) THEN
          DO 15 J=MFIT+1,MA
            YM=YM-A(LISTA(J))*AFUNC(LISTA(J))
15        CONTINUE
        ENDIF
        SIG2I=1./SIG(I)**2
        DO 17 J=1,MFIT
          WT=AFUNC(LISTA(J))*SIG2I
          DO 16 K=1,J
            COVAR(J,K)=COVAR(J,K)+WT*AFUNC(LISTA(K))
16        CONTINUE
          BETA(J)=BETA(J)+YM*WT
17      CONTINUE
18    CONTINUE
      IF (MFIT.GT.1) THEN
        DO 21 J=2,MFIT
          DO 19 K=1,J-1
            COVAR(K,J)=COVAR(J,K)
19        CONTINUE
21      CONTINUE
      ENDIF
      CALL GAUSSJ(COVAR,MFIT,NCVM,BETA,1,1)
      DO 22 J=1,MFIT
        A(LISTA(J))=BETA(J)
22    CONTINUE
      CHISQ=0.
      DO 24 I=1,NDATA
        CALL FUNCS(X(I),AFUNC,MA)
        SUM=0.
        DO 23 J=1,MA
          SUM=SUM+A(J)*AFUNC(J)
23      CONTINUE
        CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
24    CONTINUE
      CALL COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      RETURN
C                        ****************************
      END

      SUBROUTINE GRDATD (NCUR, NX, NY, NUMPT, XPT, YPT)
C                                       Loads Array with data points
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*1 A, SYMBOL(10)
      COMMON / GRAPH / A(130,130)
      COMMON / RANGED / XMIN, XMAX, DX, YMIN, YMAX, DY
      DIMENSION XPT(NUMPT), YPT(NUMPT, NCUR)
      DATA SYMBOL /'1','2','3','4','5','6','7','8','9','0'/
      DO 20 I = 1, NUMPT
         IX =  (XPT(I) - XMIN) / DX + 0.5
         IF (IX.LT.1) IX = 1
         IF (IX.GT.NX) IX = NX
         DO 10 J = 1, NCUR
            YY = YPT(I, J)
            IY = (YMAX - YY) / DY + 0.5
            IF (IY.LT.1) IY = 1
            IF (IY.GT.NY) IY = NY

            A(IX,IY) = SYMBOL(J)
10       CONTINUE
20    CONTINUE
      RETURN
C                       ****************************
      END
C
      SUBROUTINE GRNULL(NX,NY)
C                                               LOADS ARRAY WITH BLANK SPACES
      CHARACTER*1  A, B
      COMMON / GRAPH / A(130,130)
      DATA B/' '/
      DO 10 I = 1,NX
         DO 10 J = 1,NY
            A (I, J) = B
10    CONTINUE
      RETURN
C                       ****************************
      END
C
      SUBROUTINE GRPRNT(NX,NY)

      CHARACTER*1 A, AA, BB(0:9), CC
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / GRAPH / A(130,130)
      DATA AA,CC /'-','|'/
      DATA BB /'0','1','2','3','4','5','6','7','8','9' /
C
      WRITE(NOUT,99) AA,(AA, I=1,NX),AA
      DO 30 I = 1,NY
         IPRIME = MOD(NY+1 - I,10)
         WRITE(NOUT,99) BB(IPRIME), (A(J,I), J = 1,NX), CC
30    CONTINUE
      WRITE(NOUT,99) AA,(BB(MOD(J,10)), J = 1,NX), AA
99    FORMAT(' ', 130A1)
      RETURN
C                       ----------------------------
C
      ENTRY GRFILE (NNX, NNY, NUNIT)
C
C      OPEN (NUNIT, STATUS='OLD')
      WRITE (NUNIT, 9) AA,(AA, I=1,NNX),AA
C
      DO 31 I = 1, NNY
         IPRIME = MOD( NNY+1 - I,10)
         WRITE (NUNIT, 9) BB(IPRIME), (A(J,I), J = 1,NNX), CC
   31 CONTINUE
C
      WRITE (NUNIT, 9) AA,(BB(MOD(J,10)), J = 1,NNX), AA
   9  FORMAT(' ', 130A1)
C
      RETURN
C                       ****************************
      END
C                 ===========================================
C                                 LIBDINT.FOR
C     DEC 02 89
C                        ----------------------------
C                        ****************************
C           MODIFIED VERSION USING NEW ADAPTIVE END-POINT PROCEDURE

C                        ----------------------------
C List of Subprograms:

C     FUNCTION   ADZINT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZSPL (F, I, IER)
C     SUBROUTINE ADZCAL (F,I)
C     SUBROUTINE SGLINT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOTALZ
C     FUNCTION   INTUSZ ()
C
C     COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
C    > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), NUMINT,
C    > ICTA, ICTB, FA, FB, IB
C                        ****************************

      SUBROUTINE GSER(GAMSER,A,X,GLN)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)PAUSE
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
C                        ****************************
      END

      SUBROUTINE INPT2I (V1, V2, PROMPT, VN1, VX1, VN2, VX2, IR)
C
C                       Given 'Prompt', reads 2 Integer input numbers from
C                       terminal as the values for the REAL dummy variables
C                       V1 - V2.   The values are checked against
C                       the allowed ranges (VnI, VxI) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
      IMPLICIT INTEGER (T, V)
      CHARACTER*(*) PROMPT, REPLY*4
      LOGICAL TEST
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DIMENSION T(2), TN(2), TX(2)
C
      TEST = .TRUE.
      IF (IR .EQ. 99) TEST = .FALSE.
      IR = 0

      T(1) = V1
      T(2) = V2
      TN(1) = VN1
      TN(2) = VN2
      TX(1) = VX1
      TX(2) = VX2
      CALL RTB (PROMPT, LEN)
C
    1 KPRT = 0
      PRINT '(/1X, A/1X, A, 2I10)', PROMPT (1:LEN),
     >' or type / to keep the defaults :   ', (T(I), I=1,2)
C
      READ (NIN, *, ERR=99) (T(I), I=1,2)
C
      IF (.NOT.TEST) GOTO 8
      DO 10 I = 1, 2
      IF (T(I) .LT. TN(I) .OR. T(I) .GT. TX(I)) THEN
         PRINT '(1X, A, I2, A, I10, A, I10, a)',
     >  'Input #', I, ' outside the range [', TN(I), ', ', TX(I), ']'
         KPRT = 1
      ENDIF
   10 CONTINUE
C
      IF (KPRT .EQ. 1) THEN
    2    PRINT *,
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 99) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 V1 = T(1)
      V2 = T(2)
C
      RETURN
C
   99 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C
C               *************************
      END
C
      SUBROUTINE INPT2R (V1, V2, PROMPT, VN1, VX1, VN2, VX2, IR)
C
C                       Given 'Prompt', reads 2 REAL*8 input numbers from
C                       terminal as the values for the REAL dummy variables
C                       V1 - V2.   The values are checked against
C                       the allowed ranges (VnI, VxI) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*(*) PROMPT, REPLY*4
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DIMENSION T(2), TN(2), TX(2)
C
      ITST = 1
      IF (IR .EQ. 99) ITST = 0
      IR = 0

      T(1) = V1
      T(2) = V2
      TN(1) = VN1
      TN(2) = VN2
      TX(1) = VX1
      TX(2) = VX2
      CALL RTB (PROMPT, LEN)
C
    1 KPRT = 0
      PRINT '(/1X, A/A, 2(1pE15.3))', PROMPT (1:LEN),
     >      ' or type / to keep the defaults :   ', (T(I), I=1,2)
C
      READ (NIN, *, ERR=99) (T(I), I=1,2)
C
      IF (ITST .EQ. 0) GOTO 8
      DO 10 I = 1, 2
      IF (T(I) .LT. TN(I) .OR. T(I) .GT. TX(I)) THEN
         PRINT '(1X, A, I2, A, 1pE10.3, A, E10.3, a)',
     >  'Input #', I, ' outside the range [', TN(I), ', ', TX(I), ']'
         KPRT = 1
      ENDIF
   10 CONTINUE
C
      IF (KPRT .EQ. 1) THEN
    2    PRINT *,
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 99) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 V1 = T(1)
      V2 = T(2)
C
      RETURN
C
   99 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C
C                       ****************************
      END
C
      SUBROUTINE INPT3I (V1, V2, V3, PROMPT,
     >                   VN1, VX1, VN2, VX2, VN3, VX3, IR)
C
C                       Given 'Prompt', reads 3 Integer input numbers from
C                       terminal as the values for the REAL dummy variables
C                       V1 - V3.   The values are checked against
C                       the allowed ranges (VnI, VxI) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
      IMPLICIT INTEGER (T, V)
      CHARACTER*(*) PROMPT, REPLY*4
      LOGICAL TEST
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DIMENSION T(3), TN(3), TX(3)
C
      TEST = .TRUE.
      IF (IR .EQ. 99) TEST = .FALSE.
      IR = 0

      IR = 0
      T(1) = V1
      T(2) = V2
      T(3) = V3
      TN(1) = VN1
      TN(2) = VN2
      TN(3) = VN3
      TX(1) = VX1
      TX(2) = VX2
      TX(3) = VX3
      CALL RTB (PROMPT, LEN)
C
    1 KPRT = 0
      PRINT '(/1X, A/1X, A, 3I10)', PROMPT (1:LEN),
     >' or type / to keep the defaults :   ', (T(I), I=1,3)
C
      READ (NIN, *, ERR=99) (T(I), I=1,3)
C
      IF (.NOT.TEST) GOTO 8
      DO 10 I = 1, 3
      IF (T(I) .LT. TN(I) .OR. T(I) .GT. TX(I)) THEN
         PRINT '(1X, A, I2, A, I10, A, I10, a)',
     >  'Input #', I, ' outside the range [', TN(I), ', ', TX(I), ']'
         KPRT = 1
      ENDIF
   10 CONTINUE
C
      IF (KPRT .EQ. 1) THEN
    2    PRINT *,
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 99) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 V1 = T(1)
      V2 = T(2)
      V3 = T(3)
C
      RETURN
C
   99 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C
C               *************************
      END
C
      SUBROUTINE INPT3R (V1, V2, V3, PROMPT,
     >VN1, VX1, VN2, VX2, VN3, VX3, IR)
C
C                       Given 'Prompt', reads 3 REAL*8 input numbers from
C                       terminal as the values for the REAL dummy variables
C                       V1 - V3.   The values are checked against
C                       the allowed ranges (VnI, VxI) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*(*) PROMPT, REPLY*4
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DIMENSION T(3), TN(3), TX(3)
C
      ITST = 1
      IF (IR .EQ. 99) ITST = 0
      IR = 0

      T(1) = V1
      T(2) = V2
      T(3) = V3
      TN(1) = VN1
      TN(2) = VN2
      TN(3) = VN3
      TX(1) = VX1
      TX(2) = VX2
      TX(3) = VX3
      CALL RTB (PROMPT, LEN)
C
    1 KPRT = 0
      PRINT '(/1X, A/ A, 3(1pE15.3))', PROMPT (1:LEN),
     >      ' or type / to keep the defaults : ', (T(I), I=1,3)
C
      READ (NIN, *, ERR=99) (T(I), I=1,3)
C
      DO 10 I = 1, 3
      IF (T(I) .LT. TN(I) .OR. T(I) .GT. TX(I)) THEN
         PRINT '(1X, A, I2, A, 1pE10.3, A, E10.3, a)',
     >  'Input #', I, ' outside the range [', TN(I), ', ', TX(I), ']'
         KPRT = 1
      ENDIF
   10 CONTINUE
C
      IF (KPRT .EQ. 1) THEN
    2    PRINT *,
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 99) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 V1 = T(1)
      V2 = T(2)
      V3 = T(3)
C
      RETURN
C
   99 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C
C                       ****************************
      END
C
      SUBROUTINE INPT5R (V1, V2, V3, V4, V5, PROMPT,
     >VN1, VX1, VN2, VX2, VN3, VX3, VN4, VX4, VN5, VX5, IR)
C
C                       Given 'Prompt', reads 5 REAL*8 input numbers from
C                       terminal as the values for the REAL dummy variables
C                       V1 - V5.   The values are checked against
C                       the allowed ranges (VnI, VxI) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*(*) PROMPT, REPLY*4
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DIMENSION T(5), TN(5), TX(5)
C
      ITST = 1
      IF (IR .EQ. 99) ITST = 0
      IR = 0

      T(1) = V1
      T(2) = V2
      T(3) = V3
      T(4) = V4
      T(5) = V5
      TN(1) = VN1
      TN(2) = VN2
      TN(3) = VN3
      TN(4) = VN4
      TN(5) = VN5
      TX(1) = VX1
      TX(2) = VX2
      TX(3) = VX3
      TX(4) = VX4
      TX(5) = VX5
      CALL RTB (PROMPT, LEN)
C
    1 KPRT = 0
      PRINT '(/1X, A/ A / 5(1pE15.3))', PROMPT (1:LEN),
     >      ' or type / to keep the defaults : ', (T(I), I=1,5)
C
      READ (NIN, *, ERR=99) (T(I), I=1,5)
C
      IF (ITST .EQ. 0) GOTO 8
      DO 10 I = 1, 5
      IF (T(I) .LT. TN(I) .OR. T(I) .GT. TX(I)) THEN
         PRINT '(1X, A, I2, A, 1pE10.3, A, E10.3, a)',
     >  'Input #', I, ' outside the range [', TN(I), ', ', TX(I), ']'
         KPRT = 1
      ENDIF
   10 CONTINUE
C
      IF (KPRT .EQ. 1) THEN
    2    PRINT *,
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 99) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 V1 = T(1)
      V2 = T(2)
      V3 = T(3)
      V4 = T(4)
      V5 = T(5)
C
      RETURN
C
   99 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C                       ****************************
      END

      SUBROUTINE INPUTC (NAME, PROMPT, IR)
C
C                       Given 'Prompt', reads character input from
C                       terminal as the value for the char variable
C                       'Name'.  A return of '/' leave 'Name' unchanged
C
      CHARACTER*(*) NAME, PROMPT, TEM*80
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      IR = 0

      CALL RTB (PROMPT, LEN)
    1 PRINT '(/ 1X, A/ 2A)', PROMPT (1:LEN),
     >      ' or type / to keep the default :   ', NAME
C
      READ (NIN, '(A)', ERR=10) TEM
C
      IF (TEM(1:1) .NE. '/') NAME = TEM
C
      RETURN
C
   10 WRITE (NOUT, *) 'Error in INPUTC, Try again!'
      GOTO 1
C
C               *************************
      END

      SUBROUTINE INPUTI (IVALUE, PROMPT, IMIN, IMAX, IR)
C
C                       Given 'Prompt', reads INTEGER value input
C                       from terminal as the value for the dummy
C                       variable IVALUE.   The answer is checked against
C                       the allowed range (Imin, Imax) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
C                       If IR = 99 upon input, the limits are ignored.

      CHARACTER*(*) PROMPT, REPLY*1
      LOGICAL TEST
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      TEST = .TRUE.
      IF (IR .EQ. 99) TEST = .FALSE.
      IR = 0

      CALL RTB (PROMPT, LEN)

      ITEM = IVALUE
    1 PRINT '(/1X, A/A, I8)', PROMPT (1:LEN),
     >      ' or type / to keep the default :   ', ITEM
C
      READ (NIN, *, ERR=10) ITEM
C
      IF (.NOT.TEST) GOTO 8
      IF (ITEM .LT. IMIN .OR. ITEM .GT. IMAX) THEN
         PRINT *, 'Input outside the range [', IMIN, ', ', IMAX, ']'
    2    PRINT *,
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 10) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 IVALUE = ITEM
C
      RETURN
C
   10 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C
C               *************************
      END
C
      SUBROUTINE INPUTR (VALUE, PROMPT, VMIN, VMAX, IR)
C
C                       Given 'Prompt', reads DOUBLE PRECISION input
C                       from terminal as the value for the dummy
C                       variable VALUE.   The answer is checked against
C                       the allowed range (Vmin, Vmax) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
C          An input value of IR = 99 will cause the limit test to be by-passed.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*(*) PROMPT, REPLY*4
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      ITST = 1
      IF (IR .EQ. 99) ITST = 0
      IR = 0

      CALL RTB (PROMPT, LEN)
      TEM = VALUE
    1 PRINT '(/1X, A/A, 1PE15.4)', PROMPT (1:LEN),
     >      ' or type / to keep the default :   ', TEM
C
      READ (NIN, *, ERR=10) TEM
C
      IF (ITST .EQ. 0) GOTO 8
      IF (TEM .LT. VMIN .OR. TEM .GT. VMAX) THEN
         PRINT *, 'Input outside the range [', VMIN, ', ', VMAX, ']'
    2    PRINT *,
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 10) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 VALUE = TEM
C
      RETURN
C
   10 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C
C                       ****************************
      END
C
      SUBROUTINE INQUIR (PROMPT, REPLY, LEN)
      CHARACTER*(*) PROMPT, REPLY
C
      COMMON /IOUNIT/NIN, NOUT, NWRT
C
      WRITE (NOUT, 100) PROMPT
      READ (NIN, 110) REPLY
      CALL RTB (REPLY, LEN)
      RETURN
C
 100  FORMAT ('$', A, ': ')
 110  FORMAT (A)
C               *************************
      END
C
      SUBROUTINE INQURE (NAME, PROMPT, IR)
C
C                       Given 'Prompt', reads character input from
C                       terminal as the value for the char variable
C                       'Name'.  A return of '/' leave 'Name' unchanged
C
      CHARACTER*(*) NAME, PROMPT, TEM*80
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
    1 WRITE (NOUT, '(/ ''$'', 2A)') PROMPT, ' or type / to keep the
     > default :   '
C
      READ (NIN, '(A)', ERR=10) TEM
C
      IF (TEM(1:1) .NE. '/') NAME = TEM
C
      RETURN
C
   10 WRITE (NOUT, *) 'Error in INQURE, Try again!'
      GOTO 1
C
C               *************************
      END
C
      FUNCTION INT2SZ ()
C                    Return number of intervals used in last call to ADZ2NT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
      INT2SZ = NUMINT
      RETURN
C                        ****************************
      END
C
      FUNCTION INTUSZ ()
C                    Return number of intervals used in last call to ADZINT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      INTUSZ = NUMINT
      RETURN
C                        ****************************
      END
C
      SUBROUTINE LINFIT (X,Y,NDATA,SIG,MWT, A,B,SIGA,SIGB,CHI2,Q)

C                                           Routine FIT from Numerical Recipes
C     Linear fit to Y(i) = A * X(i) + B

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION X(NDATA), Y(NDATA), SIG(NDATA)

      SX=0.
      SY=0.
      ST2=0.
      B=0.
      IF(MWT.NE.0) THEN
        SS=0.
        DO 11 I=1,NDATA
          WT=1./(SIG(I)**2)
          SS=SS+WT
          SX=SX+X(I)*WT
          SY=SY+Y(I)*WT
11      CONTINUE
      ELSE
        DO 12 I=1,NDATA
          SX=SX+X(I)
          SY=SY+Y(I)
12      CONTINUE
        SS=FLOAT(NDATA)
      ENDIF
      SXOSS=SX/SS
      IF(MWT.NE.0) THEN
        DO 13 I=1,NDATA
          T=(X(I)-SXOSS)/SIG(I)
          ST2=ST2+T*T
          B=B+T*Y(I)/SIG(I)
13      CONTINUE
      ELSE
        DO 14 I=1,NDATA
          T=X(I)-SXOSS
          ST2=ST2+T*T
          B=B+T*Y(I)
14      CONTINUE
      ENDIF
      B=B/ST2
      A=(SY-SX*B)/SS
      SIGA=SQRT((1.+SX*SX/(SS*ST2))/SS)
      SIGB=SQRT(1./ST2)
      CHI2=0.
      IF(MWT.EQ.0) THEN
        DO 15 I=1,NDATA
          CHI2=CHI2+(Y(I)-A-B*X(I))**2
15      CONTINUE
        Q=1.
        SIGDAT=SQRT(CHI2/(NDATA-2))
        SIGA=SIGA*SIGDAT
        SIGB=SIGB*SIGDAT
      ELSE
        DO 16 I=1,NDATA
          CHI2=CHI2+((Y(I)-A-B*X(I))/SIG(I))**2
16      CONTINUE
        AAA = 0.5*(NDATA-2)
        BBB = 0.5*CHI2
        Q=GAMMQ(AAA, BBB)
      ENDIF
      RETURN
      END
      SUBROUTINE LOGGAM(X, U)
C
C       SUBROUTINE LOGGAM--TRANSCRIBED FROM NYU FAP ROUTINE OF MAX
C       GOLDSTEIN WRITTEN FOR FORTRAN IV
C       CALL IS CALL LOGGAM(X, U) WHERE X IS REAL
C        AND U IS THE FUNCTION VALUE
C        MUTILATION TO REAL VALUES COURTESY OF PORTER JOHNSON JUNE 1982
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION H(7)
C
      DATA H/ 2.69488974, 1.517473649, 1.011523068, 0.525606469,
     1          0.2523809524, 0.033333333, 0.0833333333 /
C
      DATA E2/ 1.5709632679 /, E8/3.14159265359/
C
      B2 = 0.D0
      J = 2
      X2 = X
3797    IF(X) 2794, 9999, 100
100   T = X**2
5793    B7 =  T
C
C       REAL PART OF LOG
C
      T1 = 5D-1 * LOG(B7)
      IF(X .GE. 2.D0) GO TO 3793
      B2 = B2 + T1
      X  = X  + 1.D0
      J  = 1
      GOTO 3797
C
3793    T3 = T1* (X - 5D-1) -X + 0.9189385332
      T4 = X
      T1 = B7
      DO 200 I = 1,7
      T = H(I)/T1
      T4 = T*T4 + X
200     T1 = T4**2
C
      T3 = T4 -X + T3
      GO TO (8795,4794) , J
C
8795    T3 = T3 - B2
C
4794    IF(X2) 4796, 4795, 4795
4795    U = T3
      X = X2
      RETURN
C
4796    U = T3 -E4
      X = X2
      RETURN
C
C               X IS NEGATIVE
C
2794    E4 = 0.D0
      IE6 = 0
5797    E4 = E4 + LOG(ABS(X))
      IE6 = IE6 + 1
      X = X + 1.D0
      IF( X .LT. 0.D0) GO TO 5797
      GO TO 3797
C
9999    WRITE(1,9990) X2
9990    FORMAT(' ATTEMPTED TO TAKE LOGGAMOF X = ',F12.5)
C
      STOP
      END
C                       *********************************
C
      FUNCTION NEXTUN()
C                                    Returns an unallocated FORTRAN i/o unit.
      LOGICAL EX
C
      DO 1  N = 10, 300
          INQUIRE (UNIT=N, OPENED=EX)
          IF (.NOT. EX) then
             nextun = n
             RETURN
          end if
   1  CONTINUE
      RETURN
C               *************************
      END
C
      FUNCTION PI()
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (P=3.1415927)
      ENTRY PIFN()
      PI = P
      RETURN
      END
C                       ****************************
C
      SUBROUTINE PLT4D (NCUR, NUMPT, XPT, YPT, NAMEX, NAMEY, NUU)
C
C       03 10 87

C                       Main Switching routine for plotting
C                       Determines ranges of variables
C                       Sorts data points in descending order in y
C                       Performs re-scaling if  desired
C                       Switches to selected plotting device
C                       Returns to calling program
C
C                       Mxxpt is the maximum limit of Ncur * Numpt

C       Input parameters:

C       Ncur    :       Number of curves = Number of columns of Y-array

C       Numpt   :       Number of x-points

C       Xpt     :       Array of x - values   [ Xpt (1 : Numpt) ]

C       Ypt     :       2-dim array of y - values    [Numpt x Ncur]

C       NameX   :       character const. = Name of x - variable

C       NameY   :       character array  = Names of y - variables

C       Nuu     :       Unit number to which records of plots and
C                       tables are to be written (on demand).
C
C                       Nuu = 0  means that no output file exists
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MXCUR = 20, MXXPT = 2500)

      CHARACTER*130 COMENT, HEADNG
      CHARACTER*2 PLOTER, AL*4, FLNM*15
      CHARACTER*(*) NAMEX, NAMEY(NCUR), NMY*20
      CHARACTER*11 NAMX, NAMY(MXCUR), NAM
      LOGICAL ASK, LDPD, ZX, ZY

      DIMENSION SCALE(MXCUR)
      DIMENSION XPT(NUMPT), YPT(NUMPT, NCUR)
      DIMENSION SX (MXXPT), SY (MXXPT)

      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / RANGED / XMIN, XMAX, DX, YMIN, YMAX, DY

      DATA FLNM, BIGNUM, XSML, YSML / 'P00', 1E30, 1E-5, 1E-5 /
      DATA PLOTER, SML / 'CH', 1E-5 /
C                                Initiations
      NU = NUU
      AL = 'None'
      IACT = 0

      DO 19 I = 1, MXCUR
19       SCALE(I)= 1.
C                                Process the labels of fit the heading format
      DO 50 J = 0, NCUR
         IF (J .NE. 0) THEN
            NMY = NAMEY(J)
         ELSE
            NMY = NAMEX
         ENDIF

         CALL TRMSTR (NMY, KNM)
         KPD = (11-KNM)/2
         IF (KPD .LE. 0) THEN
            NAM = NMY
         ELSE
            DO 55 I = 1, KPD+1
               NAM(I:I) = ' '
   55       CONTINUE
            NAM(KPD+1:) = NMY
         ENDIF

         IF (J .EQ. 0) THEN
            NAMX = NAM
         ELSE
            NAMY(J) = NAM
         ENDIF
   50 CONTINUE
C                    Make sure that there is enough room for data points

      IF ( (NCUR * NUMPT)  .GT. MXXPT) THEN
         PRINT *,
     >   'Input data points exceed maximum array size in PLOTT routine'
         STOP
      ENDIF
C                       SORTD does two main tasks:
C                       It computes the limits of x- and y- coordinates
C                       and gives the results in  RANGED Common.
C                       It sorts the Y- values into descending order,
C                       and outputs the results as ONE-DIM arrays.

    1 CALL SORTD (NCUR, NUMPT, XPT, YPT, SX, SY, IACT)
      PRINT *, 'The limits for the X,Y-axes and the curves are:'

      WRITE (NOUT, 931) XMAX, XMIN, YMAX, YMIN
  931 FORMAT (/' Xmax = ', T10, 1PE10.2, 10X, 'Xmin = ', T40, E10.2
     >        /' Ymax = ', T10, 1PE10.2, 10X, 'Ymin = ', T40, E10.2 )
      WRITE (NOUT, *) '------------------------------------------------'
C
      WRITE (NOUT, 933)
     >    ( J, SY(1 + (J-1)*NUMPT), J, SY(NUMPT*J), J = 1, NCUR )
  933 FORMAT(' Y',I2,'mx = ' , T10, 1PE10.2, 10X,
     >          'Y',I2,'mn = ' , T40, 1PE10.2 )
      WRITE (NOUT, *)

      IF (ASK('Do you wish to re-scale any of the curves')) THEN
         PRINT *, 'Enter the scaling factors, in order: ',
     >                   'Terminate with / <cr>'
         READ *, SCALE
         DO 71 JCUR = 1, NCUR
            DO 81 IX = 1, NUMPT
                YPT(IX, JCUR) = YPT(IX, JCUR) * SCALE(JCUR)
                SY (IX + (JCUR-1)*NUMPT) =
     >                  SY (IX + (JCUR-1)*NUMPT) * SCALE(JCUR)
81           CONTINUE
71       CONTINUE
         IACT = 1
         GOTO 1
      ENDIF

      YX = YMAX
      YN = YMIN
701   WRITE (NOUT, *) 'Choose the Y-range for the graph; defaults ',
     >  'are:  Ymax, Ymin  ='
      WRITE (NOUT, *) YMAX, YMIN

      READ (NIN, *, ERR=701) YX, YN
      IF (YN .GE. YX) THEN
         WRITE (NOUT, *) 'Ymin > Ymax not allowed, try again!'
         GOTO 701
      ENDIF
C      If (Yn .Lt. Ymin-Abs(Ymin*Sml) .or. Yx .Gt. Ymax+Abs(Ymax*Sml))
C     >   Then
C         Write (Nout, *) 'New range cannot exceed original range;'
C         Write (Nout, *) 'Defaults are:', Ymax, Ymin, '   Try again!'
C         Goto 701
C      Endif

      DO 72 JC = 1, NCUR
         DO 82 IX = 1, NUMPT
            YPT (IX, JC) = MAX (YPT(IX, JC), YN)
            YPT (IX, JC) = MIN (YPT(IX, JC), YX)
            ST = SY (IX + (JC - 1) * NUMPT)
            ST = MAX (ST, YN)
            ST = MIN (ST, YX)
            SY (IX + (JC - 1) * NUMPT) = ST
  82     CONTINUE
  72  CONTINUE

      YMIN = YN
      YMAX = YX

   7  CALL PLTCHD (NCUR, NUMPT, XPT, YPT, NAMEX, NAMEY, SCALE, AL, NU)

C                               Restore the scale factors and Log flag
      DO 29 I = 1, MXCUR
   29    SCALE(I)= 1.

      IF (ASK('Do you wish to print out the results in tabular form'))
     >   THEN
         NCR = MIN (NCUR, 6)
         WRITE (NOUT, 1003) NAMX, (NAMY(I), I = 1, NCR)
 1003    FORMAT(/ 3X, 6(1X, A11, 1X) /)

         DO 543 IPT = 1, NUMPT
            WRITE (NOUT, 1004) XPT(IPT), (YPT(IPT, ICUR), ICUR=1, NCR)
  543    CONTINUE
 1004    FORMAT ( 1PE13.4, 6E13.4 )

         IF (NCUR .GT. 6) THEN
            WRITE (NOUT, 1005) NAMX, (NAMY(I), I = 6, NCUR)
 1005       FORMAT(// 2X, 6(1X, A11, 1X) /)

            DO 544 IPT = 1, NUMPT
               WRITE (NOUT, 1006) XPT(IPT), (YPT(IPT,ICUR),ICUR=6,NCUR)
 1006          FORMAT( 1PE13.4, 6E13.4 )
  544       CONTINUE
         ENDIF
         IF (NU .NE. 0) THEN
            LDPD=ASK('Do you wish to dump this table to the extl file')
            IF (LDPD) THEN
  401          WRITE (NOUT, *) 'Input Heading for the table: '
               WRITE (NOUT, *)
               READ (NIN, '(A)', IOSTAT=IOS, ERR=501) HEADNG
  501          IF (IOS .GT. 0) THEN
                  PRINT *, 'Read Error, IOSTAT = ', IOS, ' Try again!'
                  GOTO 401
               ENDIF
               WRITE (NU, 1013) HEADNG, NAMX, (NAMY(I), I = 1, NCR)
 1013          FORMAT(/ 1X, A / 2X, 6(1X, A11, 1X) /)

               DO 556 IPT = 1, NUMPT
                  WRITE (NU, 1014) XPT(IPT), (YPT(IPT,ICUR), ICUR=1,NCR)
 1014             FORMAT( 1PE13.4, 6E13.4 )
  556          CONTINUE

               IF (NCUR .GT. 6) THEN
                  WRITE (NU, 1015) NAMX, (NAMY(I), I = 6, NCUR)
 1015             FORMAT(/ 2X, 6(1X, A11, 1X) /)
C
                  DO 557 IPT = 1, NUMPT
                     WRITE (NU,1016)
     >                  XPT(IPT),(YPT(IPT,ICUR),ICUR=6,NCUR)
 1016                FORMAT( 1PE13.4, 6E13.4 )
  557             CONTINUE
               ENDIF
            ENDIF
         ENDIF
      ENDIF
C                                       Initiation for the next section
      AL = 'None'
      ZX = .TRUE.
      ZY = .TRUE.

      IF ( ASK('Do you wish to convert either axis to Log-scale') )
     >          THEN
   46    WRITE (NOUT, 941)
941      FORMAT(' Type   LY   for log-scale y-axis,' /
     >                 '   "    LX   for log-scale x-axis,' /
     >                 '   "    LL   for log-log plot.    ')
         READ  (NIN, '(A)')  AL
         CALL UPCASE (AL)

C           XSML = MAX(XMIN, XSML)
C           YSML = MAX(YSML, YMIN)
         WRITE  (NOUT, 1101) XMIN, YMIN
 1101    FORMAT (
     >     ' If x or y becomes negative, taking log is illegal;' /
     >     ' Alternatively, if Xmin or Ymin is very close to zero, ',
     >         'the logarithm'/
     >    ' may extend too far into the negative range than desirable.'
     >    / ' Current values of XMIN and YMIN are: ',2(1PE15.2))

         CALL INPT2R (XSML, YSML, 'Specify cutoffs for X & Y',
     >          SML*SML, XMAX, 1.D-20, YMAX, J9)

C     >   /' Specify desired cutoff for  X  and  Y ;    Defaults are: '/
C     >   2(1PE15.2) )
C         READ (NIN, *) XSML, YSML

         IF (AL .EQ. 'LX' .OR. AL .EQ. 'LL') THEN
            IF (XMIN .LT. 0) PRINT *, 'Negative argument in Log X.'
            DO 61 IX = 1, NUMPT
               XX = XPT(IX)
               IF (XX .LT. XSML) THEN
                  IF (ZX) THEN
                     PRINT *, 'X < XSMALL; X set = ', XSML
                     ZX = .FALSE.
                  ENDIF
                  XX = XSML
               ENDIF
               XPT (IX) = LOG10 (XX)
61          CONTINUE
            IACT = 1
         ENDIF
         IF (AL .EQ. 'LY' .OR. AL .EQ. 'LL') THEN
94          IF (YMIN .LE. 0)  THEN
               PRINT *,
     >          'Negative argument in Log Y, converting to Log AbsY'
            ENDIF
            DO 602 JC = 1, NCUR
               DO 62 IX = 1, NUMPT
                  YY = ABS( YPT(IX,JC) )
                  YS = ABS( SY(IX+(JC-1)*NUMPT) )
                  IF (YY .LT. YSML) THEN
                     IF (ZY) THEN
                        PRINT *, 'Y < YSMALL; set Y = ', YSML
                        ZY = .FALSE.
                     ENDIF
                     YY = YSML
                  ENDIF
                  IF (YS .LT. YSML) YS = YSML
                  YPT(IX,JC) = LOG10 (YY)
                  SY(IX+(JC-1)*NUMPT) = LOG10 (YS)
62             CONTINUE
602         CONTINUE
            IACT = 1
         ENDIF

         IF (IACT .EQ. 0) THEN
            PRINT *, 'Illegal Choice of LX, LY and LL, Try again!'
            GOTO 46
         ELSE
            GOTO 1
         ENDIF
      ENDIF

  705 IF (ASK('Add a comment to the above results for the record'))THEN
         WRITE (NOUT, *)
         READ (NIN, '(A)', ERR=707) COMENT
         WRITE (NU, *) COMENT
         WRITE (NOUT, *)
      ENDIF

      RETURN

707   WRITE (NOUT, *) 'Error reading input; Try again!'
      GOTO 705
C                       ****************************
      END

      SUBROUTINE PLTCHD
     >(NCUR, NUMPT, XPT, YPT, NAMEX, NAMEY, SCALE, AL, NUU)
C
C                               Produces Plots of Ncur curves in ASC charaters
C
C               Input parameters:
C
C               Ncur, Numpt, Xpt, Ypt, NameX, NameY, Nuu:    Same as in PLOT2
C
C               Scale   :       Number array of Scale factor for the cureves
C
C               Al      :       Flag for Logrithmic X- and/or Y-axis
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*(*) AL
      CHARACTER*(*) NAMEX, NAMEY(NCUR)
      CHARACTER HEADNG*130, FLNM*15
      LOGICAL LDP, ASK
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / RANGED / XMIN, XMAX, DX, YMIN, YMAX, DY
C
      DIMENSION SCALE(10)
      DIMENSION XPT(NUMPT), YPT(NUMPT, NCUR)
C
      DATA NUNIT, FLNM, NX, NY / 80, 'P00', 76, 20 /
C
      NU = NUU
C
      WRITE (NOUT, 99)
99    FORMAT(/ ' The following functions are being plotted' / )
      DO 30 I = 1, NCUR
         IF (I .EQ. 10) THEN
            J = 0
         ELSE
            J = I
         ENDIF
         WRITE (NOUT, 98)  NAMEY(I), J
98       FORMAT(' The function ', A,' --- symbol is ', I1)
30    CONTINUE
C
40    WRITE(NOUT,999) NUMPT
999   FORMAT( 38X, ' Number of data points is = ', I4, //
     >  ' Put in Nx (x-mesh), Ny (y-mesh), ',
     >  'or  Type  /  to use current default values.' / )
      PRINT *, 'Default:  Nx =', NX,  '         Ny =', NY
      READ (NIN,*) NX, NY
      IF ( (NX.LT.10) .OR. (NX.GT.130) ) GOTO 101
C
      DX = (XMAX - XMIN) / FLOAT(NX-1)
      DY = (YMAX - YMIN) / FLOAT(NY-1)
      IF (DY.LE.0) GOTO 50
C
      CALL GRNULL (NX, NY)
      CALL GRDATD (NCUR, NX, NY, NUMPT, XPT, YPT)
      CALL GRPRNT (NX, NY)
C
      WRITE(NOUT,997) XMIN, NAMEX, XMAX,   YMIN, YMAX
997     FORMAT( 1PE12.3, '<= X (', A, ') <= ', 1PE12.3,
     >          T50, 1PE12.3, '<= Y <= ', 1PE12.3 )
C
      IF (NU .NE. 0) THEN
C
         LDP = ASK('Do you wish to dump this graph to the extl file')
         IF (LDP) THEN
  401       WRITE (NOUT, *) 'Input Heading for the Plot:'
            READ  (NIN, '(A)', ERR=401) HEADNG
C
            WRITE (NU, *) HEADNG
            CALL GRFILE (NX, NY, NU)
            WRITE (NU, 977) XMIN, NAMEX, XMAX,   YMIN, YMAX
977         FORMAT( 1PE12.3, '<= X (', A, ') <= ', 1PE12.3,
     >          T49, 1PE12.3, '<= Y <= ', 1PE12.3 )
            IF (AL .EQ. 'LX' .OR. AL .EQ. 'LL')
     >          WRITE (NU, *) 'X-axis in Ln scale'
            IF (AL .EQ. 'LY' .OR. AL .EQ. 'LL') THEN
               DO 97 I = 1, NCUR
                  NAMEY(I) = 'Ln '//NAMEY(I)
97             CONTINUE
            ENDIF
            WRITE (NU, '(/)')
            WRITE (NU, 198) (I, NAMEY(I), SCALE(I), I = 1, NCUR)
198         FORMAT ( 1X, I2, ' :  ', A, '*', F8.1 )
         ENDIF
      ENDIF
C
50    CONTINUE
C
      IF (ASK (
     >    'Do you wish to change the meshsize for the same plot'))
     >  GOTO 40
C
      RETURN
C
101   WRITE (NOUT, 901)
901   FORMAT ( ' Input mesh-size outside the acceptable range!'//
     >         ' 10 < Nx < 130 .    Try Again!'/)
      GOTO 40
C                       ****************************
      END
C
      SUBROUTINE PLTIND (MXCUR, MXPT, NCUR, NUMPT, YP, YPT)
C               Converts the two-dim array Yp (Mxpt, Mxcur) with actual non-
C               zero dimension (Numpt, Ncur) into an one-dim array Ypt with
C               no unwanted zeros in the middle.   Helps avoid confusion and
C               mis-handling of data in the plotting routines which follow.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION YP(MXPT, MXCUR), YPT(MXPT*MXCUR)
      I = 0
      DO 10 ICUR = 1, NCUR
         DO 20 IX   = 1, NUMPT
           I = I + 1
           YPT(I) = YP(IX, ICUR)
20       CONTINUE
10    CONTINUE
      RETURN
C                       ****************************
      END
C
      SUBROUTINE POLIN2(X1A,X2A,YA,M,N,X1,X2,Y,DY)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=20,MMAX=20)
      DIMENSION X1A(M),X2A(N),YA(M,N),YNTMP(NMAX),YMTMP(MMAX)

      DO 12 J=1,M
        DO 11 K=1,N
          YNTMP(K)=YA(J,K)
11      CONTINUE
        CALL POLINT(X2A,YNTMP,N,X2,YMTMP(J),DY)
12    CONTINUE
      CALL POLINT(X1A,YMTMP,M,X1,Y,DY)
      RETURN
      END
C                        ****************************

      SUBROUTINE POLINT (XA,YA,N,X,Y,DY)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.) PAUSE
c          THEN
c             print *,X
c             PAUSE
c             endif
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

C     ****************************

      SUBROUTINE RATINTO(XA,YA,N,X,Y,DY)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=10,TINY=1.E-25)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      HH=ABS(X-XA(1))
      DO 11 I=1,N
        H=ABS(X-XA(I))
        IF (H.EQ.0.)THEN
          Y=YA(I)
          DY=0.0
          RETURN
        ELSE IF (H.LT.HH) THEN
          NS=I
          HH=H
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)+TINY
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          W=C(I+1)-D(I)
          H=XA(I+M)-X
          T=(XA(I)-X)*D(I)/H
          DD=T-C(I+1)
          IF(DD.EQ.0.)PAUSE
          DD=W/DD
          D(I)=C(I+1)*DD
          C(I)=T*DD
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
C                        ****************************


C                        ****************************

      SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=10,TINY=1.E-25,MXX=1050)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
	
      NS=1
      HH=ABS(X-XA(1))
      DO 11 I=1,N
        H=ABS(X-XA(I))
        IF (H.EQ.0.)THEN
          Y=YA(I)
          DY=0.0
          RETURN
        ELSE IF (H.LT.HH) THEN
          NS=I
          HH=H
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)+TINY
C	print *, 99, YA(I), I
11    CONTINUE
      Y=YA(NS)
C	print *, 100, Y, NS
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          W=C(I+1)-D(I)
          H=XA(I+M)-X
          T=(XA(I)-X)*D(I)/H
          DD=T-C(I+1)
          IF(DD.EQ.0.) GOTO 12
          DD=W/DD
          D(I)=C(I+1)*DD
          C(I)=T*DD
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
C	print *, 101, Y
13    CONTINUE
      RETURN
	
      END
C                        ****************************
      SUBROUTINE RBL (NU1, NU2, ILIN)

C                 Remove Trailing blanks    = 0   and put one blank at the end
C                 Remove Blank Lines:  ILIN = 1   removes all blank lines;
C                                             2   removes multiple blank lines
C                 Truncate to 80 record len = 3
C                 Also removes trailing blanks from all lines.

C                 NU1 : unit # of input  file to be processed.
C                 NU2 : unit # of output file
      CHARACTER LINE*132

      IF (ILIN .EQ. 3) THEN
        JLIN = -1
      ELSE
        JLIN = ILIN
      ENDIF

      Rewind (NU1)
      Rewind (NU2)

      DO 3 I = 1, 10000

      DO 5 IBLK = 0, 500
        Read (NU1, '(A)', END=10) LINE
        CALL RTB (LINE, LEN)
        IF (JLIN .EQ. -1) LEN = MIN (LEN, 80)
        IF (LEN .NE. 0) GOTO 6
    5 CONTINUE

    6 IF (IBLK .GT. 0 .AND. JLIN .NE. 1) WRITE (NU2, *)

      IF (JLIN .EQ. 0) THEN
        WRITE (NU2, '(A)') LINE(1:LEN)//' '
      ELSE
        WRITE (NU2, '(A)') LINE(1:LEN)
      ENDIF

    3 CONTINUE

   10 CONTINUE

      RETURN
C                        ****************************
      END
      FUNCTION ROMINT(F, A, B, N, ESTER, IRET)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

C                                       Romberg Integration of function F from
C                                       A to B using 2 ** N points.
C
C                                       EstEr is the estimated error.
C                                       Iret  is a return code
C                                            0   ---         O.K.
C                                            1   ---   N   is too big
C                                            2   ---   Lower limit A is bigger
C                                                      then upper limit B
      PARAMETER (MXN = 10)
      external f
C
      DIMENSION R(2, MXN+1)
C
      IRET = 0
      IF (N .GT. MXN) THEN
         IRET = 1
         RETURN
      ENDIF
C
      IF     (A .GT. B) THEN
         IRET = 2
         RETURN
      ELSEIF (A .EQ. B) THEN
         ROMINT = 0
         RETURN
      ENDIF
C
      H = B - A
      R(1, 1) = (F(A) + F(B)) * H / 2.D0
C
      DO 10 I = 2, N+1
C
         RT = 0
         XT = A - H / 2.D0
C
         DO 15 K = 1, 2**(I-2)
            RT = RT + F(XT + H * K)
   15    CONTINUE
C
         R(2, 1) = (R(1, 1) + H * RT) / 2.D0
C
         DO 16 J = 2, I
            AJ = 4.D0 ** (J - 1)
            R(2, J) = (AJ * R(2, J-1) - R(1, J-1)) / (AJ - 1.D0)
   16    CONTINUE
C
         H = H / 2.D0
         DO 17 J = 1, I
            R(1, J) = R(2, J)
   17    CONTINUE
C
   10 CONTINUE
C
      ROMINT = R(1, N+1)
      ESTER  = R(1, N+1) - R(1, N)
C
      RETURN
C                        ****************************
      END
C
      SUBROUTINE RTB (ch, lench)
C                           Set LENCH = length of CH, not counting trailing
C                           blanks and nulls
      CHARACTER*(*) ch
      integer lenmax,lench
C               slightly change excluding trailing Ctrl key
      lenmax=LEN(ch)
      do lench=lenmax,1,-1
      if(ch(lench:lench).gt.' ') return
      enddo
      lench=0
      return
C               *************************
      END

      SUBROUTINE SGL2NT (IACT, F1, F2, F3, DX, FINT, ESTER)

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E20 /
C                                                         Use quadratic formula
      TEM = DX * (4.*F1 + 3.*F2 + 2.*F3) / 9.
C                 Error est based on Diff between quadratic and linear integrals
      ER  = DX * (4.*F1 - 6.*F2 + 2.*F3) / 9.

C                          Invoke adaptive singular parametrization if IACT = 2
C                      Algorithm is based on the formula F(x) = AA + BB * x **CC
C                 where AA, BB & CC are determined from F(Dx/4), F(Dx/2) & F(Dx)

      IF (IACT .EQ. 2) THEN
          T1 = F2 - F1
          T2 = F3 - F2
          IF (T1*T2 .LE. 0.) GOTO 7
          T3  = T2 - T1
          IF (ABS(T3)*HUGE .LT. T1**2) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -D1)  GOTO 7
          BB  = T1**2 / T3
          AA  = (F1*F3 - F2**2) / T3
C                                          Estimated integral based on A+Bx**C
          TMP = DX * (AA + BB* 4.**CC / (CC + 1.))
C                                       Error estimate based on the difference
          ER = TEM - TMP
C                                              Use the improved integral value
          TEM= TMP
      ENDIF

    7 FINT = TEM
      ESTER= ER
      RETURN
C                        ****************************
      END

      SUBROUTINE SGLINT (IACT, F1, F2, F3, DX, FINT, ESTER)

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E20 /
C                                                         Use quadratic formula
      TEM = DX * (4.*F1 + 3.*F2 + 2.*F3) / 9.
C                 Error est based on Diff between quadratic and linear integrals
      ER  = DX * (4.*F1 - 6.*F2 + 2.*F3) / 9.

C                          Invoke adaptive singular parametrization if IACT = 2
C                      Algorithm is based on the formula F(x) = AA + BB * x **CC
C                 where AA, BB & CC are determined from F(Dx/4), F(Dx/2) & F(Dx)

      IF (IACT .EQ. 2) THEN
          T1 = F2 - F1
          T2 = F3 - F2
          IF (T1*T2 .LE. 0.) GOTO 7
          T3  = T2 - T1
          IF (ABS(T3)*HUGE .LT. T1**2) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -D1)  GOTO 7
          BB  = T1**2 / T3
          AA  = (F1*F3 - F2**2) / T3
C                                          Estimated integral based on A+Bx**C
          TMP = DX * (AA + BB* 4.**CC / (CC + 1.))
C                                       Error estimate based on the difference
          ER = TEM - TMP
C                                              Use the improved integral value
          TEM= TMP
      ENDIF

    7 FINT = TEM
      ESTER= ER
      RETURN
C                        ****************************
      END

      FUNCTION SIMP (NX, DX, F)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DIMENSION F(NX)
C
      IF (NX .EQ. 1) RETURN
      IF (NX .LT. 0 .OR. NX .GT. 10000) GOTO 99
      IF (NX .GT. 4) GOTO 50
      GOTO (20, 30, 40), NX-1
C
   20 SIMP = (F(1) + F(2)) * DX / 2D0
      RETURN
   30 SIMP = (F(1) + 4D0 * F(2) + F(3)) * DX / 3D0
      RETURN
   40 SIMP = (( F(1) + 4D0 * F(2) +     F(3)) / 3D0
     >       +(-F(2) + 8D0 * F(3) + 5D0 * F(4)) / 12D0 ) * DX
      RETURN
C
   50 SE = F(2)
      SO = 0
      NM1 = NX - 1
      DO 60 I = 4, NM1, 2
      IM1 = I - 1
      SE = SE + F(I)
      SO = SO + F(IM1)
   60 CONTINUE
      MS = MOD (NX, 2)
C
      IF (MS .EQ. 1) THEN
         SIMP = (F(1) + 4D0 * SE + 2D0 * SO + F(NX)) * DX / 3D0
      ELSE
         SIMP = (F(1) + 4D0 * SE + 2D0 * SO + F(NM1)) * DX / 3D0
     >         + (-F(NM1-1) + 8D0 * F(NM1) + 5D0 * F(NX)) * DX / 12D0
      END IF
C
      RETURN
C
   99 WRITE (NOUT, 999) NX
  999 FORMAT (/ 5X, 'NX = ', I6,
     >              'out of range in D-SIMP INTEGRATION ROUTINE')
      STOP
C                        ****************************
      END

C
C     subroutine to smooth evolved distribution
C

      subroutine smooth (nset,iknl,mxx,del,iv,xar,qu)

      implicit double precision (a-h,o-z)

      dimension xar(0:mxx),qu(mxx)

       DIMENSION XTEMP(10), TEMP(10), TEMP1(10)

       external ratint

C
C     compute slope just before and after x = del
C

      sl01 = (qu(iv-2)-qu(iv-1))/(xar(iv-2)-xar(iv-1))
      sl11 = (qu(iv+1)-qu(iv+2))/(xar(iv+1)-xar(iv+2))
C
C     set # of point used in extrapolation
c
      ktemp = 7
      n1 = 7

C
C     comp and comp1 give tolerance range in which the distribution can change 
C     from the two bins iv+1 and iv-1 to the iv bin
C     different values used depending wther polarized (iknl=1) or unpolarized 
c     (iknl=2) evolution is used!
C

      if (del.lt.0.01001) then

            if (iknl.gt.0) then

            comp = 1.1
            comp1 = 0.9

            else
            
            comp = 1.02
            comp1 = 0.98
            
            endif

            else

            if (iknl.gt.0) then

            comp = 1.05
            comp1 = 0.95

            else
            
            comp = 1.02
            comp1 = 0.98
            
            endif   

          endif

C
C     first check point. if function is 0 at x = del then extrapolate from
C     DGLAP region.
C
      if (qu(iv).eq.0D0) then

         do 355 ntemp = 1,ktemp

            temp(ntemp) = qu(iv+ntemp)
            xtemp(ntemp) = 1D0*(ntemp+1)
 355        continue

            call ratint(xtemp,temp,ktemp,1D0,res,err)

            qu(iv) = res

            endif
C
C     smoothe distribution if distribution not equal 0 zt x = del!
C

      if (abs(qu(iv)/qu(iv-1)).gt.comp
     >.or.abs(qu(iv)/qu(iv-1)).lt.comp1
     >) then    

 405     do 455 ntemp = 1,ktemp

            temp1(ntemp) = qu(iv-ntemp)
            xtemp(ntemp) = 1D0*(ntemp)
 455        continue

            call ratint(xtemp,temp1,ktemp,0D0,res,err)

            qu(iv) = res
C
C     checke whether first extrapolation produced somthing sensible, if not
C     redo it with different # of points in extrapolation
C
      if (abs(qu(iv)/qu(iv-1)).gt.comp
     >.or.abs(qu(iv)/qu(iv-1)).lt.comp1
     >) then

C     reduce # of extrapolation points

         ktemp = ktemp - 1

         goto 405

         endif

C     restore # of extrapolation points to old value

         ktemp = n1

        endif

C
C     final check whether slopes of curves and relative values at the bins
C     iv-1,iv+1 and iv are ok with respect to one another i.e. 
C     no abnormalities. Idea: use second order polynomial with two
C     support points, once in the DGLAP and once in the ERBL region
C     determine the value of the dsitribution at x = del both times and 
C     take the average.
C

        if ((sl01.lt.0D0.and.sl11.lt.0D0.and.
     >  qu(iv).gt.qu(iv-1).and.qu(iv).gt.
     >  qu(iv+1)).or.(sl01.gt.0D0.and.sl11.gt.0D0.and.
     >  qu(iv).lt.qu(iv-1).and.qu(iv).lt.
     >  qu(iv+1)).or.(sl01.lt.0D0.and.sl11.lt.0D0.and.
     >  qu(iv).lt.qu(iv-1).and.qu(iv).lt.
     >  qu(iv+1)).or.(sl01.gt.0D0.and.sl11.gt.0D0.and.
     >  qu(iv).gt.qu(iv-1).and.qu(iv).gt.
     >  qu(iv+1))) then

      DET =1D0/((xar(IV+2) - xar(IV-1))*(xar(IV+2) - xar(IV+1))) 
      DET1 =1D0/((xar(IV+1) - xar(IV-1))*(xar(IV+2) - xar(IV+1))) 
      DET2 =1D0/((xar(IV+1) - xar(IV-1))*(xar(IV+2) - xar(IV-1))) 
C 
      FF1 = QU(IV+2)*DET -QU(IV+1)*DET1 +QU(IV-1)*DET2 
C 
      FF2 = - (xar(IV-1) + xar(IV+1))*QU(IV+2)*DET  
     >          + (xar(IV+2) + xar(IV-1))*QU(IV+1)*DET1  
     >          - (xar(IV+1) + xar(IV+2))*QU(IV-1)*DET2 
C 
      FF3 = + xar(IV+1)*xar(IV-1)*QU(IV+2)*DET  
     >          - xar(IV-1)*xar(IV+2)*QU(IV+1)*DET1  
     >          + xar(IV+2)*xar(IV+1)*QU(IV-1)*DET2

      DET00 =1D0/((xar(IV-2) - xar(IV+1))*(xar(IV-2) - xar(IV-1))) 
      DET11 =1D0/((xar(IV-1) - xar(IV+1))*(xar(IV-2) - xar(IV-1))) 
      DET22 =1D0/((xar(IV-1) - xar(IV+1))*(xar(IV-2) - xar(IV+1))) 
C 
      F1 = QU(IV-2)*DET00 -QU(IV-1)*DET11 +QU(IV+1)*DET22 
C 
      F2 = - (xar(IV+1) + xar(IV-1))*QU(IV-2)*DET00 
     >          + (xar(IV-2) + xar(IV+1))*QU(IV-1)*DET11  
     >          - (xar(IV-1) + xar(IV-2))*QU(IV+1)*DET22 
C 
      F3 = + xar(IV-1)*xar(IV+1)*QU(IV-2)*DET00 
     >          - xar(IV+1)*xar(IV-2)*QU(IV-1)*DET11  
     >          + xar(IV-2)*xar(IV-1)*QU(IV+1)*DET22


      tempq1 = del**2*FF1+del*FF2+FF3

      tempq2 = del**2*F1+del*F2+F3

      qu(iv) = (tempq1+tempq2)/2D0

            endif

       if (iknl.gt.0.and.nset.eq.0) then

       qu(1) = - qu(iv)
       qu(2) = - qu(iv-1)

       elseif (iknl.gt.0.and.nset.eq.1) then

       qu(1) = qu(iv)
       qu(2) = qu(iv-1)

       elseif (iknl.lt.0.and.nset.eq.0) then

       qu(1) = qu(iv)
       qu(2) = qu(iv-1)

       elseif (iknl.lt.0.and.nset.eq.1) then

       qu(1) = - qu(iv)
       qu(2) = - qu(iv-1)

       endif

      return

      end

C
C     subroutine to smooth evolved distribution in ERBL region
C

      subroutine smooth1 (iknl,mxx,del,iv,xar,qu)

      implicit double precision (a-h,o-z)

      dimension xar(0:mxx),qu(mxx)

       DIMENSION XTEMP(10), TEMP(10), TEMP1(10)

       external ratint

C
C     compute slope in the four bins before x = del
C

      sl01 = (qu(iv-1)-qu(iv))/(xar(iv-1)-xar(iv))
      sl11 = (qu(iv-3)-qu(iv-2))/(xar(iv-3)-xar(iv-2))
C
C     set # of point used in extrapolation
c
      ktemp = 7
      n1 = 7

C
C     first check point. if sign changes. If so extrapolate! 
C     
C
 777  if (qu(iv-2).gt.0D0.and.qu(iv-1).lt.0D0.or.
     > qu(iv-2).lt.0D0.and.qu(iv-1).gt.0D0) then

         do 356 ntr = 1,0,-1

         do 355 ntemp = 1,ktemp

            temp(ntemp) = qu(iv-ntemp-ntr)
            xtemp(ntemp) = 1D0*(ntemp+1)
 355        continue

            call ratint(xtemp,temp,ktemp,1D0,res,err)

            qu(iv-ntr) = res

 356        continue

            ktemp = ktemp-1

            goto 777

            endif

C
C     reset number of extrapolation points
C

            ktemp = 7

C
C     smoothe distribution if slopes of the last four bins are different!
C

 888  if  (sl01.gt.0D0.and.sl11.lt.0D0.or.
     > sl01.lt.0D0.and.sl11.gt.0D0)then    

         do 456 itr = 1,0,-1

        do 455 ntemp = 1,ktemp

            temp1(ntemp) = qu(iv-ntemp-itr)
            xtemp(ntemp) = 1D0*(ntemp)
 455        continue

            call ratint(xtemp,temp1,ktemp,0D0,res,err)

            qu(iv-itr) = res

 456        continue

            sl01 = (qu(iv-1)-qu(iv))/(xar(iv-1)-xar(iv))
            sl11 = (qu(iv-3)-qu(iv-2))/(xar(iv-3)-xar(iv-2))

            ktemp = ktemp -1

            goto 888

            endif

      return

      end



C                      ================================
C                              LIBDMTH.FOR
C                        ----------------------------

      FUNCTION SMPNOL (NX, DX, FN, ERR)
C                                            DP Left-Open Simpson Integration
C     Inputs:  Nx, Dx, Fn(1:Nx)                             [F(1) is not used]
C     Output:  Err                                          (Error estimate)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FN(NX)

      MS = MOD(NX, 2)
      IF (NX .LE. 1 .OR. NX .GT. 1000) THEN
         PRINT *, 'NX =', NX, ' OUT OF RANGE IN SMPNOL!'
         STOP
      ELSEIF (NX .EQ. 2) THEN
         TEM = DX * FN(2)
      ELSEIF (NX .EQ. 3) THEN
         TEM = DX * FN(2) * 2.
      ELSE
         IF (MS .EQ. 0) THEN
            TEM = DX * (23.* FN(2) - 16.* FN(3) + 5.* FN(4)) / 12.
            TMP = DX * (3.* FN(2) - FN(3)) / 2.
            ERR = ABS(TEM - TMP)
            TEM = TEM + SMPSNA (NX-1, DX, FN(2), ER1)
            ERR = ABS(ER1) + ERR
         ELSE
            TEM = DX * (8.* FN(2) - 4.* FN(3) + 8.* FN(4)) / 3.
            TMP = DX * (3.* FN(2) + 2.* FN(3) + 3.* FN(4)) / 2.
            ERR = ABS(TEM - TMP)
            TEM = TEM + SMPSNA (NX-4, DX, FN(5), ER1)
c		print *, FN(5)
c		read (*,*)
            ERR = ABS(ER1) + ERR
         ENDIF
      ENDIF

      SMPNOL = TEM
      RETURN
C                        ****************************
      END
C
      FUNCTION SMPNOR (NX, DX, FN, ERR)
C                                              DP Right-Open Simpson Integration
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FN(NX)

      MS = MOD(NX, 2)
      IF (NX .LE. 1 .OR. NX .GT. 1000) THEN
         PRINT *, 'NX =', NX, ' OUT OF RANGE IN SMPNOR!'
         STOP
      ELSEIF (NX .EQ. 2) THEN
         TEM = DX * FN(Nx-1)
      ELSEIF (NX .EQ. 3) THEN
         TEM = DX * FN(NX-1) * 2.
      ELSE
        IF (MS .EQ. 0) THEN
         TEM = DX * (23.* FN(NX-1) - 16.* FN(NX-2) + 5.* FN(NX-3)) / 12.
         TMP = DX * (3.* FN(NX-1) - FN(NX-2)) / 2.
         ERR = ABS(TEM - TMP)
         TEM = TEM + SMPSNA (NX-1, DX, FN(1), ER1)
         ERR = ER1 + ERR
        ELSE
         TEM = DX * (8.* FN(NX-1) - 4.* FN(NX-2) + 8.* FN(NX-3)) / 3.
         TMP = DX * (3.* FN(NX-1) + 2.* FN(NX-2) + 3.* FN(NX-3)) / 2.
         ERR = ABS(TEM - TMP)
         TEM = TEM + SMPSNA (NX-4, DX, FN(1), ER1)
         ERR = ER1 + ERR
        ENDIF
      ENDIF

      SMPNOR = TEM
      RETURN
C                        ****************************
      END

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
        CALL WARNR(IW2,NWRT,'DX cannot be < 0. in SMPSNA', 'DX', DX,
     >               D0, D1, 0)
        SMPSNA = 0.
        RETURN
      ENDIF

      IF (NX .LE. 0 .OR. NX .GT. MAXX) THEN
        CALL WARNI(IW1, NWRT, 'NX out of range in SMPSNA', 'NX', NX,
     >               1, MAXX, 1)
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

      FUNCTION SMPSNF (FN, A, B, NX, ERR, IER)
C
C                       Does integral of F(X)dx from A TO B by the SIMPSON METHO
C
C                       Double precision version of SMPSN
C
C                       Input:          External function:      FN
C                                       Lower limit      :      A
C                                       Upper limit      :      B
C                                       Number of points :      Nx
C
C                       Uses (Nx-1) evenly spaced intervals.
C
C                       Output:         error estimate:         ERR
C                                       error code    :         IER
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      PARAMETER (MXPT = 1000)
C
      DIMENSION X(MXPT)
      external fn

C
      IF (NX .EQ. 1) RETURN
      IF (NX .LT. 0 .OR. NX .GT. MXPT) GOTO 99
C
      DX = (B - A) / (NX-1)
      IF (DX .LE. 0) THEN
        WRITE (NOUT, *) 'DX .LE. 0 in SMPSNF, DX =', DX
        SMPSNF = 0
        RETURN
      ENDIF
C
      DO 10 I = 1, NX
      X(I) = (A*(NX-I) + B*(I-1)) / (NX-1)
   10 CONTINUE
C
      IF (NX .GT. 4) GOTO 50
C
      GOTO (20, 30, 40), NX-1
   20 SMPSNF = (FN(X(1)) + FN(X(2))) * DX / 2D0
      RETURN
   30 SMPSNF = (FN(X(1)) + 4D0 * FN(X(2)) + FN(X(3))) * DX / 3D0
      RETURN
   40 SMPSNF = (( FN(X(1)) + 4D0 * FN(X(2)) +     FN(X(3))) / 3D0
     > + (-FN(X(2)) + 8D0 * FN(X(3)) + 5D0 * FN(X(4))) / 12D0 ) * DX
      RETURN
C
   50 SE = FN(X(2))
      SO = 0
      NM1 = NX - 1
      DO 60 I = 4, NM1, 2
      IM1 = I - 1
      SE = SE + FN(X(I))
      SO = SO + FN(X(IM1))
   60 CONTINUE
      MS = MOD (NX, 2)
      IF (MS .EQ. 1) THEN
        SMPSNF = (FN(X(1)) + 4D0 * SE + 2D0 * SO + FN(X(NX))) * DX / 3D0
        TRPZ = (FN(X(1)) + 2D0 * (SE + SO) + FN(X(NX))) * DX / 2D0
      ELSE
        SMPSNF =(FN(X(1)) + 4D0 * SE + 2D0 * SO + FN(X(NM1))) * DX / 3D0
     > +(-FN(X(NM1-1)) + 8D0 * FN(X(NM1)) + 5D0 * FN(X(NX))) * DX / 12D0
        TRPZ = (FN(X(1)) + 2D0 * (SE + SO + FN(X(NM1))) + FN(X(NX)))
     >          * DX / 2D0
      ENDIF
C
      ERR = SMPSNF - TRPZ
C
      RETURN
C
   99 WRITE (NOUT, 999) NX
  999 FORMAT (/ 5X, 'NX = ', I6,
     >              'out of range in SIMPSON INTEGRATION ROUTINE')
      STOP
C                        ****************************
      END
C
      SUBROUTINE SORTD (NCUR, NUMPT, XPT, YPT, SX, SY, IACT)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION XPT(NUMPT), YPT(NUMPT, NCUR)
      DIMENSION SX (NUMPT, NCUR), SY (NUMPT, NCUR)

      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / RANGED / XMIN, XMAX, DX, YMIN, YMAX, DY

      DATA BIGNUM / 1.0 E30 /

      IF (IACT .NE. 0 .AND. IACT .NE. 1)
     >   WRITE (NOUT, *) 'Illegel value of Iact in SORTD, Iact =', IACT

      IF (IACT .EQ. 0) THEN
         DO 10 IC = 1, NCUR
            SY(1, IC) = YPT(1, IC)
            SX(1, IC) = XPT(1)
            DO 20 IP = 2, NUMPT
               DO 30 IP1 = 1, IP - 1
                  IF (YPT(IP, IC) .GE. SY(IP1, IC)) THEN
                     YTEM = SY(IP1, IC)
                     XTEM = SX(IP1, IC)
                     SY(IP1, IC) = YPT(IP, IC)
                     SX(IP1, IC) = XPT(IP)
                     DO 40 IP2 = IP1 + 1, IP
                        XX = SX(IP2, IC)
                        YY = SY(IP2, IC)
                        SX(IP2, IC) = XTEM
                        SY(IP2, IC) = YTEM
                        XTEM = XX
                        YTEM = YY
   40                CONTINUE
                     GOTO 20
                  END IF
   30          CONTINUE
               SX(IP, IC) = XPT(IP)
               SY(IP, IC) = YPT(IP, IC)
   20       CONTINUE
   10    CONTINUE
      ENDIF
C                               Find range in x (that is, find xmin and xmax)
1     XMIN =  BIGNUM
      XMAX = -BIGNUM
      DO 31 I = 1,NUMPT
         XMIN = MIN(XMIN, XPT(I))
         XMAX = MAX(XMAX, XPT(I))
31    CONTINUE
c                               Find range in y (ymin and ymax)
      YMIN =  BIGNUM
      YMAX = -BIGNUM
      DO 21 J = 1, NCUR
         YMIN = MIN (YMIN, SY(NUMPT, J))
         YMAX = MAX (YMAX, SY(1, J))
21    CONTINUE

      RETURN
C                       ****************************
      END
C
      SUBROUTINE SORTI (NCUR, NUMPT, MPT, LPT, MST, LST)
C                                Integer version of SORT1 in the PLOTT package.
C                      Sorts the Arrays MPT and LPT in decreasing order of LPT.
C                                                           Output MST and LST
      DIMENSION MPT(NUMPT), LPT(NUMPT, NCUR)
      DIMENSION MST (NUMPT, NCUR), LST (NUMPT, NCUR)

      COMMON / IOUNIT / NIN, NOUT, NWRT

      DATA BIGNUM / 999999 /

C      IF (IACT .NE. 0 .AND. IACT .NE. 1)
C     >   WRITE (NOUT, *) 'Illegel value of Iact in SORT1, Iact =', IACT

      DO 10 IC = 1, NCUR
      LST(1, IC) = LPT(1, IC)
      MST(1, IC) = MPT(1)
                DO 20 IP = 2, NUMPT
                        DO 30 IP1 = 1, IP - 1
                        IF (LPT(IP, IC) .GE. LST(IP1, IC)) THEN
                                YTEM = LST(IP1, IC)
                                XTEM = MST(IP1, IC)
                                LST(IP1, IC) = LPT(IP, IC)
                                MST(IP1, IC) = MPT(IP)
                                DO 40 IP2 = IP1 + 1, IP
                                        XX = MST(IP2, IC)
                                        YY = LST(IP2, IC)
                                        MST(IP2, IC) = XTEM
                                        LST(IP2, IC) = YTEM
                                        XTEM = XX
                                        YTEM = YY
   40                           CONTINUE
                                GOTO 20
                        END IF
   30           CONTINUE
                        MST(IP, IC) = MPT(IP)
                        LST(IP, IC) = LPT(IP, IC)
   20   CONTINUE
   10 CONTINUE
C                           Find range in x (that is, find MMin and MMax)
1     MMIN =  BIGNUM
      MMAX = -BIGNUM
      DO 31 I = 1,NUMPT
                MMIN = MIN(MMIN, MPT(I))
                MMAX = MAX(MMAX, MPT(I))
31    CONTINUE
c                               Find range in y (LMin and LMax)
      LMIN =  BIGNUM
      LMAX = -BIGNUM
      DO 21 J = 1, NCUR
        LMIN = MIN (LMIN, LST(NUMPT, J))
        LMAX = MAX (LMAX, LST(1, J))
21    CONTINUE

      RETURN
C                       **********************
      END
C                      ================================

C     
      SUBROUTINE SPLINE (XA,F,N,IY,IX1,DEL,FF)

C
C     subroutine to calculate the second order derivative of the function in
C     array ya at the points in the array xa starting at point ix. needed
C     for subroutine SPLINT.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (MXX = 1050)

      DIMENSION XA(0:MXX), F(MXX), FF(3,MXX)


      DO 10 I=IY,N-1

         
      IF(XA(N-1).GT.DEL)THEN

          IF (I.GT.IY) THEN
		
      DET =1D0/((XA(I+1) - XA(I-1))*(XA(I+1) - XA(I)))
      DET1 =1D0/((XA(I) - XA(I-1))*(XA(I+1) - XA(I)))
      DET2 =1D0/((XA(I) - XA(I-1))*(XA(I+1) - XA(I-1)))
C
      FF(1,I) = F(I+1-IX1)*DET -F(I-IX1)*DET1 +F(I-1-IX1)*DET2
C
      FF(2,I) = - (XA(I-1) + XA(I))*F(I+1-IX1)*DET 
     >          + (XA(I+1) + XA(I-1))*F(I-IX1)*DET1 
     >          - (XA(I) + XA(I+1))*F(I-1-IX1)*DET2
C
      FF(3,I) = + XA(I)*XA(I-1)*F(I+1-IX1)*DET 
     >          - XA(I-1)*XA(I+1)*F(I-IX1)*DET1 
     >          + XA(I+1)*XA(I)*F(I-1-IX1)*DET2
     
      

        ELSEIF (I.EQ.IY.AND.IX1.EQ.0)THEN

      DET =1D0/((XA(I+2) - XA(I))*(XA(I+2) - XA(I+1)))
      DET1 =1D0/((XA(I+1) - XA(I))*(XA(I+2) - XA(I+1)))
      DET2 =1D0/((XA(I+1) - XA(I))*(XA(I+2) - XA(I)))
C
      FF(1,I) = F(I+2-IX1)*DET -F(I+1-IX1)*DET1 +F(I-IX1)*DET2
C
      FF(2,I) = - (XA(I) + XA(I+1))*F(I+2-IX1)*DET 
     >          + (XA(I+2) + XA(I))*F(I+1-IX1)*DET1 
     >          - (XA(I+1) + XA(I+2))*F(I-IX1)*DET2
C
      FF(3,I) = + XA(I+1)*XA(I)*F(I+2-IX1)*DET 
     >          - XA(I)*XA(I+2)*F(I+1-IX1)*DET1 
     >          + XA(I+2)*XA(I+1)*F(I-IX1)*DET2
    
      ELSEIF(I.EQ.Y .AND. XA(IX1+1).GT.DEL) THEN

      DET =1D0/((XA(I+2) - XA(I))*(XA(I+2) - XA(I+1)))
      DET1 =1D0/((XA(I+1) - XA(I))*(XA(I+2) - XA(I+1)))
      DET2 =1D0/((XA(I+1) - XA(I))*(XA(I+2) - XA(I)))
C
      FF(1,I) = F(I+2-IX1)*DET -F(I+1-IX1)*DET1 +F(I-IX1)*DET2
C
      FF(2,I) = - (XA(I) + XA(I+1))*F(I+2-IX1)*DET 
     >          + (XA(I+2) + XA(I))*F(I+1-IX1)*DET1 
     >          - (XA(I+1) + XA(I+2))*F(I-IX1)*DET2
C
      FF(3,I) = + XA(I+1)*XA(I)*F(I+2-IX1)*DET 
     >          - XA(I)*XA(I+2)*F(I+1-IX1)*DET1 
     >          + XA(I+2)*XA(I+1)*F(I-IX1)*DET2   
      
        ELSE

      DET =1D0/((XA(I+2) - DEL)*(XA(I+2) - XA(I+1)))
      DET1 =1D0/((XA(I+1) - DEL)*(XA(I+2) - XA(I+1)))
      DET2 =1D0/((XA(I+1) - DEL)*(XA(I+2) - DEL))
C
      FF(1,I) = F(I+2-IX1)*DET -F(I+1-IX1)*DET1 +F(I-IX1)*DET2
C
      FF(2,I) = - (DEL + XA(I+1))*F(I+2-IX1)*DET 
     >          + (XA(I+2) + DEL)*F(I+1-IX1)*DET1 
     >          - (XA(I+1) + XA(I+2))*F(I-IX1)*DET2
C
      FF(3,I) = + XA(I+1)*DEL*F(I+2-IX1)*DET 
     >          - DEL*XA(I+2)*F(I+1-IX1)*DET1 
     >          + XA(I+2)*XA(I+1)*F(I-IX1)*DET2     

        ENDIF
        
        ENDIF

      IF (XA(IY).LT.DEL.AND.XA(N-1).LT.DEL) THEN

           IF (I.GT.IY.AND.XA(N).LT.DEL) THEN
		
      DET =1D0/((XA(I+1) - XA(I-1))*(XA(I+1) - XA(I)))
      DET1 =1D0/((XA(I) - XA(I-1))*(XA(I+1) - XA(I)))
      DET2 =1D0/((XA(I) - XA(I-1))*(XA(I+1) - XA(I-1)))
C
      FF(1,I) = F(I+1-IX1)*DET -F(I-IX1)*DET1 +F(I-1-IX1)*DET2
C
      FF(2,I) = - (XA(I-1) + XA(I))*F(I+1-IX1)*DET 
     >          + (XA(I+1) + XA(I-1))*F(I-IX1)*DET1 
     >          - (XA(I) + XA(I+1))*F(I-1-IX1)*DET2
C
      FF(3,I) = + XA(I)*XA(I-1)*F(I+1-IX1)*DET 
     >          - XA(I-1)*XA(I+1)*F(I-IX1)*DET1 
     >          + XA(I+1)*XA(I)*F(I-1-IX1)*DET2

      ELSEIF (I.GT.IY.AND.XA(N).GT.DEL) THEN

      DET =1D0/((DEL - XA(I-1))*(DEL - XA(I)))
      DET1 =1D0/((XA(I) - XA(I-1))*(DEL - XA(I)))
      DET2 =1D0/((XA(I) - XA(I-1))*(DEL - XA(I-1)))
C
      FF(1,I) = F(I+1-IX1)*DET -F(I-IX1)*DET1 +F(I-1-IX1)*DET2
C
      FF(2,I) = - (XA(I-1) + XA(I))*F(I+1-IX1)*DET 
     >          + (DEL + XA(I-1))*F(I-IX1)*DET1 
     >          - (XA(I) + DEL)*F(I-1-IX1)*DET2
C
      FF(3,I) = + XA(I)*XA(I-1)*F(I+1-IX1)*DET 
     >          - XA(I-1)*DEL*F(I-IX1)*DET1 
     >          + DEL*XA(I)*F(I-1-IX1)*DET2
         
           
        ELSEIF (I.EQ.IY) THEN

      DET =1D0/((XA(I+2) - XA(I))*(XA(I+2) - XA(I+1)))
      DET1 =1D0/((XA(I+1) - XA(I))*(XA(I+2) - XA(I+1)))
      DET2 =1D0/((XA(I+1) - XA(I))*(XA(I+2) - XA(I)))
C
      FF(1,I) = F(I+2-IX1)*DET -F(I+1-IX1)*DET1 +F(I-IX1)*DET2
C
      FF(2,I) = - (XA(I) + XA(I+1))*F(I+2-IX1)*DET 
     >          + (XA(I+2) + XA(I))*F(I+1-IX1)*DET1 
     >          - (XA(I+1) + XA(I+2))*F(I-IX1)*DET2
C
      FF(3,I) = + XA(I+1)*XA(I)*F(I+2-IX1)*DET 
     >          - XA(I)*XA(I+2)*F(I+1-IX1)*DET1 
     >          + XA(I+2)*XA(I+1)*F(I-IX1)*DET2

        ENDIF

        ENDIF
     
 10     CONTINUE
     
            RETURN
C
            END
C

      SUBROUTINE SPLINED (XA,F,N,IY,IX,IX1,DEL,FF)

C
C     In this routine the integrand f(y) will be approximated
C     in the bin x_i, x_i+1 by making a polynomial approximation:
C                    f(y)= a*y^2 + b*y +c
C     where the coefficients a,b,c are determined by knowing f(y) in
C     the points x_i-1, x_i, x_i+1. This relationship yields a 3x3 matrix
C     which has then to be inverted to find a,b,c! 
C     The coeffiecints for the bins are
C     stored in the array FF as FF(i) = a, FF(i+1) = b, FF(i+2) = c
C     corresponding to the bin x_i, x_i+1. 
C

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (MXX = 1050)

      DIMENSION XA(0:MXX), F(0:MXX,0:MXX), FF(3,MXX)


C
C     Loop to calculate a,b,c.
C
     

      IF (IX1.EQ.1) THEN

      DO 11 I = N,IY-1,-1

         IF (I.GT.IY) THEN
		
      DET =1D0/((XA(I+1) - XA(I-1))*(XA(I+1) - XA(I)))
      DET1 =1D0/((XA(I) - XA(I-1))*(XA(I+1) - XA(I)))
      DET2 =1D0/((XA(I) - XA(I-1))*(XA(I+1) - XA(I-1)))
C
      FF(1,I) = F(IX,I+1)*DET -F(IX,I)*DET1 +F(IX,I-1)*DET2
C
      FF(2,I) = - (XA(I-1) + XA(I))*F(IX,I+1)*DET 
     >          + (XA(I+1) + XA(I-1))*F(IX,I)*DET1 
     >          - (XA(I) + XA(I+1))*F(IX,I-1)*DET2
C
      FF(3,I) = + XA(I)*XA(I-1)*F(IX,I+1)*DET 
     >          - XA(I-1)*XA(I+1)*F(IX,I)*DET1 
     >          + XA(I+1)*XA(I)*F(IX,I-1)*DET2

     

         ELSEIF (I.EQ.IY) THEN

      DET =1D0/((XA(I+1) - DEL)*(XA(I+1) - XA(I)))
      DET1 =1D0/((XA(I) - DEL)*(XA(I+1) - XA(I)))
      DET2 =1D0/((XA(I) - DEL)*(XA(I+1) - DEL))
C
      FF(1,I) = F(IX,I+1)*DET -F(IX,I)*DET1 +F(IX,I-1)*DET2
C
      FF(2,I) = - (DEL + XA(I))*F(IX,I+1)*DET 
     >          + (XA(I+1) + DEL)*F(IX,I)*DET1 
     >          - (XA(I) + XA(I+1))*F(IX,I-1)*DET2
C
      FF(3,I) = + XA(I)*DEL*F(IX,I+1)*DET 
     >          - DEL*XA(I+1)*F(IX,I)*DET1 
     >          + XA(I+1)*XA(I)*F(IX,I-1)*DET2

      

            ELSEIF (I.EQ.IY-1) THEN

      DET =1D0/((XA(I+2) - DEL)*(XA(I+2) - XA(I+1)))
      DET1 =1D0/((XA(I+1) - DEL)*(XA(I+2) - XA(I+1)))
      DET2 =1D0/((XA(I+1) - DEL)*(XA(I+2) - DEL))
C
      FF(1,I) = F(IX,I+2)*DET -F(IX,I+1)*DET1 +F(IX,I)*DET2
C
      FF(2,I) = - (DEL + XA(I+1))*F(IX,I+2)*DET 
     >          + (XA(I+2) + DEL)*F(IX,I+1)*DET1 
     >          - (XA(I+1) + XA(I+2))*F(IX,I)*DET2
C
      FF(3,I) = + XA(I+1)*DEL*F(IX,I+2)*DET 
     >          - DEL*XA(I+2)*F(IX,I+1)*DET1 
     >          + XA(I+2)*XA(I+1)*F(IX,I)*DET2

      
            ENDIF

 11         CONTINUE

            ELSE

       DO 10 I = IY, N-1         

      IF (I.GT.IY) THEN
		
      DET =1D0/((XA(I+1) - XA(I-1))*(XA(I+1) - XA(I)))
      DET1 =1D0/((XA(I) - XA(I-1))*(XA(I+1) - XA(I)))
      DET2 =1D0/((XA(I) - XA(I-1))*(XA(I+1) - XA(I-1)))
C
      FF(1,I) = F(I+1-IX1,IX)*DET -F(I-IX1,IX)*DET1 +F(I-1-IX1,IX)*DET2
C
      FF(2,I) = - (XA(I-1) + XA(I))*F(I+1-IX1,IX)*DET 
     >          + (XA(I+1) + XA(I-1))*F(I-IX1,IX)*DET1 
     >          - (XA(I) + XA(I+1))*F(I-1-IX1,IX)*DET2
C
      FF(3,I) = + XA(I)*XA(I-1)*F(I+1-IX1,IX)*DET 
     >          - XA(I-1)*XA(I+1)*F(I-IX1,IX)*DET1 
     >          + XA(I+1)*XA(I)*F(I-1-IX1,IX)*DET2
     
      

        ELSEIF (I.EQ.IY) THEN

      DET =1D0/((XA(I+2) - XA(I))*(XA(I+2) - XA(I+1)))
      DET1 =1D0/((XA(I+1) - XA(I))*(XA(I+2) - XA(I+1)))
      DET2 =1D0/((XA(I+1) - XA(I))*(XA(I+2) - XA(I)))
C
      FF(1,I) = F(I+2-IX1,IX)*DET -F(I+1-IX1,IX)*DET1 +F(I-IX1,IX)*DET2
C
      FF(2,I) = - (XA(I) + XA(I+1))*F(I+2-IX1,IX)*DET 
     >          + (XA(I+2) + XA(I))*F(I+1-IX1,IX)*DET1 
     >          - (XA(I+1) + XA(I+2))*F(I-IX1,IX)*DET2
C
      FF(3,I) = + XA(I+1)*XA(I)*F(I+2-IX1,IX)*DET 
     >          - XA(I)*XA(I+2)*F(I+1-IX1,IX)*DET1 
     >          + XA(I+2)*XA(I+1)*F(I-IX1,IX)*DET2

      
      
        ENDIF

 10         CONTINUE

            ENDIF

            RETURN
C
            END
C

      SUBROUTINE SPLINED1 (XA,F,N,IY,IX,IX1,DEL,FF)

C
C     subroutine to calculate the second order derivative of the function in
C     array ya at the points in the array xa starting at point ix. needed
C     for subroutine SPLINT.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (MXX = 1050)

      DIMENSION XA(0:MXX), F(0:MXX,0:MXX), FF(3,MXX)

      DO 10 I=IY,N-1

      IF (I.GT.IY .AND. XA(I+1).LE.DEL) THEN
		
      DET =1D0/((XA(I+1) - XA(I-1))*(XA(I+1) - XA(I)))
      DET1 =1D0/((XA(I) - XA(I-1))*(XA(I+1) - XA(I)))
      DET2 =1D0/((XA(I) - XA(I-1))*(XA(I+1) - XA(I-1)))
C
      FF(1,I) = F(IX,I+1-IX1)*DET -F(IX,I-IX1)*DET1 +F(IX,I-1-IX1)*DET2
C
      FF(2,I) = - (XA(I-1) + XA(I))*F(IX,I+1-IX1)*DET 
     >          + (XA(I+1) + XA(I-1))*F(IX,I-IX1)*DET1 
     >          - (XA(I) + XA(I+1))*F(IX,I-1-IX1)*DET2
C
      FF(3,I) = + XA(I)*XA(I-1)*F(IX,I+1-IX1)*DET 
     >          - XA(I-1)*XA(I+1)*F(IX,I-IX1)*DET1 
     >          + XA(I+1)*XA(I)*F(IX,I-1-IX1)*DET2
     
      A = (XA(I+1)+XA(I))/2.
      Y = FF(1,I)*A**2+FF(2,I)*A+FF(3,I)

        ELSEIF (I.EQ.IY) THEN

      DET =1D0/((XA(I+2) - XA(I))*(XA(I+2) - XA(I+1)))
      DET1 =1D0/((XA(I+1) - XA(I))*(XA(I+2) - XA(I+1)))
      DET2 =1D0/((XA(I+1) - XA(I))*(XA(I+2) - XA(I)))
C
      FF(1,I) = F(IX,I+2-IX1)*DET -F(IX,I+1-IX1)*DET1 +F(IX,I-IX1)*DET2
C
      FF(2,I) = - (XA(I) + XA(I+1))*F(IX,I+2-IX1)*DET 
     >          + (XA(I+2) + XA(I))*F(IX,I+1-IX1)*DET1 
     >          - (XA(I+1) + XA(I+2))*F(IX,I-IX1)*DET2
C
      FF(3,I) = + XA(I+1)*XA(I)*F(IX,I+2-IX1)*DET 
     >          - XA(I)*XA(I+2)*F(IX,I+1-IX1)*DET1 
     >          + XA(I+2)*XA(I+1)*F(IX,I-IX1)*DET2

            
      A = (XA(I+1)+XA(I))/2.
      Y = FF(1,I)*A**2+FF(2,I)*A+FF(3,I)


        ENDIF

 10      CONTINUE
        
            RETURN
C
            END
C
      SUBROUTINE SPLINED2 (XA,F,N,IY,IX,IX1,DEL,FF)

C
C     subroutine to calculate the second order derivative of the function in
C     array ya at the points in the array xa starting at point ix. needed
C     for subroutine SPLINT.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (MXX = 1050)

      DIMENSION XA(0:MXX), F(0:MXX,0:MXX), FF(3,MXX)

      DO 10 I=IY,N-1

         IF (I.GT.IY.AND.XA(I+1).GT.DEL) THEN

      DET =1D0/((DEL - XA(I-1))*(DEL - XA(I)))
      DET1 =1D0/((XA(I) - XA(I-1))*(DEL - XA(I)))
      DET2 =1D0/((XA(I) - XA(I-1))*(DEL - XA(I-1)))
C
      FF(1,I) = F(IX,I+1-IX1)*DET -F(IX,I-IX1)*DET1 +F(IX,I-1-IX1)*DET2
C
      FF(2,I) = - (XA(I-1) + XA(I))*F(IX,I+1-IX1)*DET 
     >          + (DEL + XA(I-1))*F(IX,I-IX1)*DET1 
     >          - (XA(I) + DEL)*F(IX,I-1-IX1)*DET2
C
      FF(3,I) = + XA(I)*XA(I-1)*F(IX,I+1-IX1)*DET 
     >          - XA(I-1)*DEL*F(IX,I-IX1)*DET1 
     >          + DEL*XA(I)*F(IX,I-1-IX1)*DET2


      ELSEIF (I.GT.IY .AND. XA(I+1).LE.DEL) THEN
		
      DET =1D0/((XA(I+1) - XA(I-1))*(XA(I+1) - XA(I)))
      DET1 =1D0/((XA(I) - XA(I-1))*(XA(I+1) - XA(I)))
      DET2 =1D0/((XA(I) - XA(I-1))*(XA(I+1) - XA(I-1)))
C
      FF(1,I) = F(IX,I+1-IX1)*DET -F(IX,I-IX1)*DET1 +F(IX,I-1-IX1)*DET2
C
      FF(2,I) = - (XA(I-1) + XA(I))*F(IX,I+1-IX1)*DET 
     >          + (XA(I+1) + XA(I-1))*F(IX,I-IX1)*DET1 
     >          - (XA(I) + XA(I+1))*F(IX,I-1-IX1)*DET2
C
      FF(3,I) = + XA(I)*XA(I-1)*F(IX,I+1-IX1)*DET 
     >          - XA(I-1)*XA(I+1)*F(IX,I-IX1)*DET1 
     >          + XA(I+1)*XA(I)*F(IX,I-1-IX1)*DET2

        ELSEIF (I.EQ.IY) THEN

      DET =1D0/((XA(I+2) - XA(I))*(XA(I+2) - XA(I+1)))
      DET1 =1D0/((XA(I+1) - XA(I))*(XA(I+2) - XA(I+1)))
      DET2 =1D0/((XA(I+1) - XA(I))*(XA(I+2) - XA(I)))
C
      FF(1,I) = F(IX,I+2-IX1)*DET -F(IX,I+1-IX1)*DET1 +F(IX,I-IX1)*DET2
C
      FF(2,I) = - (XA(I) + XA(I+1))*F(IX,I+2-IX1)*DET 
     >          + (XA(I+2) + XA(I))*F(IX,I+1-IX1)*DET1 
     >          - (XA(I+1) + XA(I+2))*F(IX,I-IX1)*DET2
C
      FF(3,I) = + XA(I+1)*XA(I)*F(IX,I+2-IX1)*DET 
     >          - XA(I)*XA(I+2)*F(IX,I+1-IX1)*DET1 
     >          + XA(I+2)*XA(I+1)*F(IX,I-IX1)*DET2

        ENDIF

 10      CONTINUE
        
            RETURN
C
            END
C

      SUBROUTINE SPLINT (YA,IX,X,Y)

C
C     subroutine to interpolate the integrand in the intervall ix..ix+1
C     using the polynomial approximation. uses result from SPLINE namely ya.
C

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (MXX = 1050)

      DIMENSION YA(3,MXX)

      Y = YA(1,IX)*X**2 + YA(2,IX)*X + YA(3,IX)

      RETURN
C
      END
C


      FUNCTION TNGLE(X,Y,Z)
C                                       "Triangle Function" of the three sides.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DATA RDOFF / 1E-11 /

      AMX = MAX (X*X, Y*Y, Z*Z)

      TMP = X*X + Y*Y + Z*Z - 2.* (X*Y + Y*Z + Z*X)
      ATMP= ABS(TMP)

      IF     (ATMP .LT. AMX * RDOFF) THEN
        TMP = ATMP
      ELSEIF (TMP .LT. 0.) THEN
        PRINT '(A, 4(1PE12.3))', 'X,Y,Z, TMP =', X, Y, Z, TMP
        STOP 'Negative argument in TNGLE function; check for errors!'
      ENDIF

      TNGLE = SQRT (TMP)

      RETURN
C           ****************************
      END

      SUBROUTINE TOT2LZ
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

      SUBROUTINE TOTALZ
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END
C
C                            2nd copy of ADZINT to be used for double integrals
C List of GLOBAL Symbols

C     FUNCTION   ADZ2NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZ2PL (F, I, IER)
C     SUBROUTINE ADZ2AL (F,I)
C     SUBROUTINE SGL2NT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOT2LZ
C     FUNCTION   INT2SZ ()
C
C     COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
C    > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
C    > ICTA, ICTB, NUMINT, IB
C                        ****************************

       SUBROUTINE TRMSTR(STRING,ILEN)

C   Removes leading spaces and returns true length of a character string

       CHARACTER STRING*(*),SPACE*1

       DATA SPACE/' '/

       ILEN=0

       IF (STRING.EQ.SPACE) RETURN
C                           Remove leading spaces
1      IF (STRING(1:1).NE.SPACE) GOTO 2
          STRING=STRING(2:)
          GOTO 1
2      CONTINUE
C                           Count up trailing spaces
       DO 3 I=LEN(STRING),1,-1
          IF (STRING(I:I).NE.SPACE) THEN
             ILEN=I
             RETURN
          END IF
3      CONTINUE
C               *************************
       END

      FUNCTION TRNGLE (X,Y,Z, IRT)
C                                       "Triangle Function" of the three sides.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DATA RDOFF / 1E-11 /

      IRT = 0

      AMX = MAX (X*X, Y*Y, Z*Z)

      TMP = X*X + Y*Y + Z*Z - 2.* (X*Y + Y*Z + Z*X)
      ATMP= ABS(TMP)

      IF     (ATMP .LT. AMX * RDOFF) THEN
        TMP = ATMP
        IRT = 1
      ELSEIF (TMP .LT. 0.) THEN
        PRINT '(A, 4(1PE12.3))', 'X,Y,Z, TMP =', X, Y, Z, TMP
        STOP 'Negative argument in TRNGLE function; check for errors!'
      ENDIF

      TRNGLE = SQRT (TMP)

      RETURN
C           ****************************
      END

        SUBROUTINE UC(A)
C               Converts A to all upper case.
C               System dependent: assumes ICHAR(uc letter) - ICHAR(lc letter)
C               is constant.
        CHARACTER A*(*), C*(1)
        INTEGER I

        ENTRY UPCASE(A)
        DO 1 I=1, LEN(A)
                C = A(I:I)
C
C                       Use ASCII ordering for detecting lc:
                IF ( LGE(C, 'a') .AND. LLE(C, 'z') )THEN
                        A(I:I) = CHAR(ICHAR(C)+ICHAR('A')-ICHAR('a'))
                        ENDIF
 1              CONTINUE
        RETURN
C               *************************
        END
C
      SUBROUTINE UpC (A, La, UpA)

C  This is a variation of the two old routines UC(A) and UpCase(A). Here
C  the converted value is return to the new variable UpA, rather than
C  the input dummy variable A, to avoid bombing the program with error
C  "access violation, reason mask=04" when the routine is given a CONSTANT
C  actual argument (rather than a character variable argument).

C  The inconvenience of this new version is that the calling program must
C  declare the length of UpA explicitly; and it better be .Ge. Len(A).

C  We trim any excess trailing characters in UpA and replace with spaces.
C  To be exact, the returned value should be considered as UpA(1:La).

      CHARACTER A*(*), UpA*(*), C*(1)
      INTEGER I, La, Ld

      La = Len(A)
      Lb = Len(UpA)

      If (Lb .Lt. La) Stop 'UpCase conversion length mismatch!'

      Ld = ICHAR('A')-ICHAR('a')

      DO 1 I = 1, Lb

        If (I .Le. La) Then
         c = A(I:I)
         IF ( LGE(C, 'a') .AND. LLE(C, 'z') ) THEN

           UpA (I:I) = CHAR(Ichar(c) + ld)
         Else
           UpA (I:I) = C

         ENDIF
        Else
         UpA (I:I) = ' '
        Endif

 1    CONTINUE

      RETURN
C               *************************
      END

      SUBROUTINE WARNI (IWRN, NWRT, MSG, NMVAR, IVAB,
     >                  IMIN, IMAX, IACT)
C                                            Integer version of WarnR
      CHARACTER*(*) MSG, NMVAR

      Data Nmax / 1000 /

      IW = IWRN
      IV = IVAB

      IF  (IW .EQ. 0) THEN
         PRINT '(1X,A/1X, 2A,I10 /A,I4)', MSG, NMVAR, ' = ', IV,
     >         ' For all warning messages, check file unit #', NWRT
         IF (IACT .EQ. 1) THEN
         PRINT       '(A/2I10)', ' The limits are: ', IMIN, IMAX
         WRITE (NWRT,'(A/2I10)') ' The limits are: ', IMIN, IMAX
         ENDIF
      ENDIF

      If (Iw .LT. Nmax) Then
         WRITE (NWRT,'(1X,A/1X,2A, I10)') MSG, NMVAR, ' = ', IV
      Elseif (Iw .Eq. Nmax) Then
         Print '(/A/)', '!!! Severe Warning, Too many errors !!!'
         Print '(/A/)', '    !!! Check The Error File !!!'
         Write (Nwrt, '(//A//)')
     >     'Too many warnings, Message suppressed !!'
      Endif

      IWRN = IW + 1

      RETURN
C               *************************
      END

      SUBROUTINE WARNR (IWRN, NWRT, MSG, NMVAR, VARIAB,
     >                  VMIN, VMAX, IACT)

C   Subroutine to handle warning messages.  Writes the (warning) message
C   and prints out the name and value of an offending variable to SYS$OUT
C   the first time, and to output file unit # NWRT in subsequent times.
C
C   The switch IACT decides whether the limits (VMIN, VMAX) are active or
C   not.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      CHARACTER*(*) MSG, NMVAR
      Data Nmax / 1000 /

      IW = IWRN
      VR = VARIAB

      IF  (IW .EQ. 0) THEN
         PRINT '(1X, A/1X,2A,1PD16.7/A,I4)', MSG, NMVAR, ' = ', VR,
     >         ' For all warning messages, check file unit #', NWRT
         IF (IACT .EQ. 1) THEN
         PRINT       '(A/2(1PE15.4))', ' The limits are: ', VMIN, VMAX
         WRITE (NWRT,'(A/2(1PE15.4))') ' The limits are: ', VMIN, VMAX
         ENDIF
      ENDIF

      If (Iw .LT. Nmax) Then
         WRITE (NWRT,'(I5, 2A/1X,2A,1PD16.7)') IW, '   ', MSG,
     >                  NMVAR, ' = ', VR
      Elseif (Iw .Eq. Nmax) Then
         Print '(/A/)', '!!! Severe Warning, Too many errors !!!'
         Print '(/A/)', '    !!! Check The Error File !!!'
         Write (Nwrt, '(//A//)')
     >     '!! Too many warnings, Message suppressed from now on !!'
      Endif

      IWRN = IW + 1

      RETURN
C           ****************************
      END

      FUNCTION ZBRNT(FUNC, X1, X2, TOL, IRT)

C                          Return code  IRT = 1 : limits do not bracket a root;
C                                             2 : function call exceeds maximum
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (ITMAX = 1000, EPS = 3.E-12)
      external func

      IRT = 0
      TOL = ABS(TOL)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.)  THEN
        PRINT *, 'Root must be bracketed for ZBRNT.'
        IRT = 1
        RETURN
      ENDIF
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRNT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
      PRINT *, 'ZBRNT exceeding maximum iterations.'
      IRT = 2
      ZBRNT=B
      RETURN
      END
C                        ****************************












