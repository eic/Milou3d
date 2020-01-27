*CMZ :          28/05/2004  09.03.53  by  H1 Saclay
*-- Author :    Unknown   28/11/2003


      program gdvcs

* Main routine
*


      include 'dvcs.common'

* ------------- PYTHIA common -------------
cv      include 'pythia.inc'
      integer      N,  K(4000,5)
      real       P(4000,5), V(4000,5)
      COMMON /LUJETS/ N,K,P,V

      common/dvcsNMAX/nmax_dvcs

      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)

      external ludata,pydata 

      DATA IEV /0/
      SAVE IEV


c read cards and do bases integration
         CALL DVCS_INI(mfail)
         if (mfail.eq.1) goto 999

         call Print_asc_head(35)

      kevt = 0
      nevmax = nmax_dvcs


C--------------------GENERATE NEXT EVENT------------------------------

       open(33,file='./ascii.txt',status='unknown')
       mstu(11)=33


   10    CONTINUE

         IEV    = IEV + 1

C---generate an event
         call gener_evt


         if (iev.le.stNPRINT) then
          call lulist(1)
         endif

         call PRINT_asc(35,iev,asc_flag)

         if (iev.le.stNGEN) goto 10


C-------------------END OF EVENT GENERATION --------------------------

        call dvcs_end

999   continue
      stop
      end

