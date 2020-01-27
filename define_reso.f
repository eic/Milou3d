*CMZ :          29/05/2004  19.18.48  by  H1 Saclay
*-- Author :    H1 Saclay   29/05/2004


* --------------------------------------------------------

      SUBROUTINE DEFINE_RESO
*
* Define the Nstar resonances in the Pythia commons
* (actually the particle's names)
* cf Benno's routine DECNST :
*                 name     PDG code      name     PDG code
*                 N(1440)+    12212  N(1440)0    12112
*                 N(1520)+     2124  N(1520)0     1214
*                 N(1535)+    22212  N(1535)0    22112
*                 N(1650)+    32212  N(1650)0    32112
*                 N(1675)+     2216  N(1675)0     2116
*                 N(1680)+    12216  N(1680)0    12116
*                 N(1700)+    22124  N(1700)0    21214
*                 N(1710)+    42212  N(1710)0    42112
*                 N(1720)+    32124  N(1720)0    31214
*
* Use the compressed KCcodes 41 - (41+18)
*
* Used only in the p-dissociated case.
*
* --------------------------------------------------------

      include 'dvcs.common'


      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON /LUDAT3/ MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON /LUDAT4/ CHAF(500)
      CHARACTER*8     CHAF

      CHARACTER*8 RNAMES(18)
      INTEGER     RCODES(18)

      DATA (RNAMES(I),I=   1, 18)/
     + 'N(1440)+','N(1520)+','N(1535)+','N(1650)+',
     + 'N(1675)+','N(1680)+','N(1700)+','N(1710)+',
     + 'N(1720)+',
     + 'N(1440)0','N(1520)0','N(1535)0','N(1650)0',
     + 'N(1675)0','N(1680)0','N(1700)0','N(1710)0',
     + 'N(1720)0'
     + /


      if (stIELAS.eq.1) goto 999


      do i=1,18
       KC = 41 + (i-1)
c       KCHG(KC,4) = RCODES(i)
       CHAF(KC)   = RNAMES(i)
       kF = KC
      enddo


      goto 999

      CHAF(LUCOMP(2124))  = 'N(1520)+'
      CHAF(LUCOMP(22212)) = 'N(1535)+'
      CHAF(LUCOMP(32212)) = 'N(1650)+'
      CHAF(LUCOMP(2216))  = 'N(1675)+'
      CHAF(LUCOMP(12216)) = 'N(1680)+'
      CHAF(LUCOMP(22124)) = 'N(1700)+'
      CHAF(LUCOMP(42212)) = 'N(1710)+'
      CHAF(LUCOMP(32124)) = 'N(1720)+'

      CHAF(LUCOMP(12112)) = 'N(1440)0'
      CHAF(LUCOMP(1214))  = 'N(1520)0'
      CHAF(LUCOMP(22112)) = 'N(1535)0'
      CHAF(LUCOMP(32112)) = 'N(1650)0'
      CHAF(LUCOMP(2116))  = 'N(1675)0'
      CHAF(LUCOMP(12116)) = 'N(1680)0'
      CHAF(LUCOMP(21214)) = 'N(1700)0'
      CHAF(LUCOMP(42112)) = 'N(1710)0'
      CHAF(LUCOMP(31214)) = 'N(1720)0'


999   CONTINUE
      RETURN
      END

