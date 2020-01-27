*CMZ :          28/05/2004  09.03.53  by  H1 Saclay
*-- Author :    Unknown   28/11/2003


      program gdvcs

* Main routine
*


      include 'dvcs.common'


      common/dvcsNMAX/nmax_dvcs


      external ludata,pydata 

      DATA IEV /0/
      SAVE IEV



         CALL DVCS_INI(mfail)
         if (mfail.eq.1) goto 999



      kevt = 0
      nevmax = nmax_dvcs


C--------------------GENERATE NEXT EVENT------------------------------


   10    CONTINUE

         IEV    = IEV + 1

C---generate an event
         call gener_evt


         if (iev.le.stNPRINT) then
          call lulist(1)
         endif

         if (iev.le.stNGEN) goto 10




C-------------------END OF EVENT GENERATION --------------------------

        call dvcs_end



999   continue
      stop
      end

