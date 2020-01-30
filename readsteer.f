*CMZ :          05/10/2004  13.48.05  by  H1 Saclay
*-- Author :    Unknown   17/12/2003



      subroutine readsteer

      implicit none

      include 'dvcs.common'


      integer StdCin,StdCout


      StdCin  = 5
      StdCout = 6

      WRITE(StdCout,*) 'looking for input cards file... '
      CALL FFINIT(1000)
      OPEN(UNIT=StdCin,FILE='dvcs.steer',STATUS='unknown')
      CALL FFSET ( 'LINP' , StdCin )
      CALL FFSET ( 'LOUT' , StdCout)
      CALL FFSET ( 'SIZE' , 10)

C --- generate in exp(var)
         CALL FFKEY('EXP',stEXP,1,'INTEGER')

C --- calculate asymetries
CCC         CALL FFKEY('ASYM',stASYM,1,'INTEGER')

C --- Generator SEED
         CALL FFKEY('SEED',stISEED,1,'INTEGER')

C --- Fixed target or collider
         CALL FFKEY('FIXED',stFIXED,1,'LOGICAL')

C --- Beam Energies
         CALL FFKEY('ELEP',stELEP,1,'REAL')
         CALL FFKEY('ETARG',stETARG,1,'REAL')

C --- RAD
         CALL FFKEY('IRAD',stIRAD,1,'INTEGER')

C --- Elastic or proton dissociation
         CALL FFKEY('IELAS',stIELAS,1,'INTEGER')

C ---  Decays of resonances (when IELAS=0)
         CALL FFKEY('IRFRA',stIRFRA,1,'INTEGER')

C --- Treatment of dissociated proton in the continuum :
         CALL FFKEY('PROSPLIT',stPROSPLIT,1,'REAL')

C --- MY**2 dependence (when IELAS=0)
         CALL FFKEY('EPSM',stEPSM,1,'REAL')

C --- TINTIN
         CALL FFKEY('TINTIN',stTINTIN,1,'LOGICAL')
CCC         CALL FFKEY('SIGN',stSIGN,1,'REAL')
         CALL FFKEY('BTIN',stBTIN,1,'REAL')
         CALL FFKEY('RTIN',stRTIN,1,'REAL')
         CALL FFKEY('F2QCD',stF2QCD,1,'LOGICAL')
         CALL FFKEY('DIPOLE',stDIPOLE,1,'LOGICAL')


C --- Lepton charge
         CALL FFKEY('LCHAR',stLCHAR,1,'REAL')

C --- Polarisations
         CALL FFKEY('LPOL',stLPOL,1,'REAL')
         CALL FFKEY('TPOL',stTPOL,1,'REAL')

C --- Z/A of target
         CALL FFKEY('ZTAR',stZTAR,1,'REAL')
         CALL FFKEY('ATAR',stATAR,1,'REAL')

C --- SPIN
         CALL FFKEY('SPIN',stSPIN,1,'REAL')

C --- Integration / Generation / Both
         CALL FFKEY('IGEN',stIGEN,1,'INTEGER')

C --- Number of events to generate
         CALL FFKEY('NGEN',stNGEN,1,'INTEGER')
         CALL FFKEY('NPRINT',stNPRINT,1,'INTEGER')

C --- Parameters for BASES
         CALL FFKEY('NCALL',stNCALL,1,'INTEGER')
         CALL FFKEY('ITMX1',stITMX1,1,'INTEGER')
         CALL FFKEY('ITMX2',stITMX2,1,'INTEGER')

C --- Debug flag
         CALL FFKEY('IDEBUG',stIDEBUG,1,'INTEGER')

C --- Number of x & q2 points in the amplitudes grid
         CALL FFKEY('NXGRID',stNX,1,'INTEGER')
         CALL FFKEY('NQGRID',stNQ,1,'INTEGER')
         CALL FFKEY('NTGRID',stNT,1,'INTEGER')

C --- Process to generate
         CALL FFKEY('IPRO',stIPRO,1,'INTEGER')

C --- LO or NLO
         CALL FFKEY('IORD',stIORD,1,'INTEGER')

C --- Kinematic domain
         CALL FFKEY('XMIN',stXMIN,1,'REAL')
         CALL FFKEY('XMAX',stXMAX,1,'REAL')
         CALL FFKEY('QMIN',stQMIN,1,'REAL')
         CALL FFKEY('QMAX',stQMAX,1,'REAL')
         CALL FFKEY('TMIN',stTMIN,1,'REAL')
         CALL FFKEY('TMAX',stTMAX,1,'REAL')

C --- Bounds on MY**2 when dissociated proton
         CALL FFKEY('MYMIN',stMYMIN,1,'REAL')
         CALL FFKEY('MYMAX',stMYMAX,1,'REAL')

C --- Integration over angles
         CALL FFKEY('NSET',stNSET,1,'INTEGER')
         CALL FFKEY('PHI',stPHI,1,'REAL')
         CALL FFKEY('PPHI',stPPHI,1,'REAL')
         CALL FFKEY('THETA',stTHETA,1,'REAL')

C --- With or without Twist 3
         CALL FFKEY('TWIST3',stTWIST3,1,'LOGICAL')

C --- Parameters for t-dependence
         CALL FFKEY('ITFORM',stITFORM,1,'INTEGER')
         CALL FFKEY('BQCST',stBQ,1,'REAL')
         CALL FFKEY('BQSLOPE',stBQ_SLOPE,1,'REAL')
         CALL FFKEY('BGCST',stBG,1,'REAL')
         CALL FFKEY('BGSLOPE',stBG_SLOPE,1,'REAL')
         CALL FFKEY('Q0SQ',stQ02,1,'REAL')
         CALL FFKEY('X0DEF',stX0,1,'REAL')

C --- Kinematic cuts
         CALL FFKEY('THLMIN',stTHLMIN,1,'REAL')
         CALL FFKEY('THLMAX',stTHLMAX,1,'REAL')
         CALL FFKEY('ELMIN',stELMIN,1,'REAL')
         CALL FFKEY('THGMIN',stTHGMIN,1,'REAL')
         CALL FFKEY('THGMAX',stTHGMAX,1,'REAL')
CCC         CALL FFKEY('EGMIN',stEGMIN,1,'REAL')
         CALL FFKEY('EIMIN',stEIMIN,1,'REAL')

         CALL FFKEY('YMIN',stYMIN,1,'REAL') ! s.fazio October 2010
         CALL FFKEY('YMAX',stYMAX,1,'REAL')
         CALL FFKEY('M12MIN',stM12MIN,1,'REAL')
         CALL FFKEY('ELMB',stELMB,1,'REAL')
         CALL FFKEY('EGMB',stEGMB,1,'REAL')

      CALL FFGO

      write(6,*) 'Valeur de IDEBUG ',stIDEBUG

      return
      end

