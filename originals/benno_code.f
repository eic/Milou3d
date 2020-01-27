*CMZ :          29/05/2004  19.16.09  by  H1 Saclay
*-- Author :    H1 Saclay   29/05/2004


      integer function benno_code(ikc)

*
* Translate a Nstar resonance KC code (41, 42, ..)
* into Benno's code, for decay using DECNST
*

      INTEGER     RCODES(18)

      DATA (RCODES(I),I=   1, 18)/
     + 12212, 2124, 22212, 32212, 2216, 12216, 22124, 42212, 32124,
     + 12112, 1214, 22112, 32112, 2116, 12116, 21214, 42112, 31214
     + /


      ireso = ikc - 41 + 1
      if (ireso.le.0.or.ireso.gt.18)
     =  write(6,*) ' --- Problem in BENNO_CODE '

      benno_code = RCODES(ireso)


      return
      end
