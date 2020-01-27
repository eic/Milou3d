*CMZ :          29/05/2004  19.08.58  by  H1 Saclay
*-- Author :    H1 Saclay   29/05/2004


      integer function KC_FREE(ICODE_RES)

*
* Translate Benno's PDGcode for Nstar resonances
* into a free compressed pythia code KC
* (free : 41 - 80)

      INTEGER     RCODES(18)

      DATA (RCODES(I),I=   1, 18)/
     + 12212, 2124, 22212, 32212, 2216, 12216, 22124, 42212, 32124,
     + 12112, 1214, 22112, 32112, 2116, 12116, 21214, 42112, 31214
     + /


      kc_free = -1

      do 10 i=1,18
       j = rcodes(i)
       if (ICODE_RES.eq.j) then
        kc_free = 41 + (i-1)
        goto 11
       endif
10    continue

11    continue

      if (kc_free.eq.-1) write(6,*) '---- Problem in KC_FREE !'

      return
      end
