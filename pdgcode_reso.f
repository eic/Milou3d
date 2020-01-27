*CMZ :          28/05/2004  08.42.49  by  H1 Saclay
*-- Author :    H1 Saclay   28/05/2004


      SUBROUTINE PDGCODE_RESO(ITYPE)

      real*8 res
      real*8 x_main,q_main,phi_main,t_main,ym_main

* --- To retrieve MY (dissociated p)
      common/dvcs_OUT/res
      common/dvcs_VAR/x_main,q_main,phi_main,t_main,ym_main


      ITYPE = -1
      DMXP = sngl(ym_main)

      ibeamp = 2212     !  adapter pour cible quelconque

        IF (DMXP .LT. 1.48) THEN
          ITYPE  = SIGN (12212, IBEAMP)
        ELSE IF (DMXP .LT. 1.60) THEN
          ITYPE  = SIGN (2124, IBEAMP)
        ELSE IF (DMXP .LT. 1.90) THEN
          R = H1RN()
          IF (R .LT. 0.5) THEN
            ITYPE  = SIGN (12216, IBEAMP)
          ELSE IF (R .LT. 0.83) THEN
            ITYPE  = SIGN (22124, IBEAMP)
          ELSE
            ITYPE  = SIGN (42212, IBEAMP)
          END IF
        ELSE
          ITYPE  = SIGN (2210, IBEAMP)
        END IF

      if (itype.eq.-1) then
       write(6,*)' --- Problem in PDGCODE_RESO !!!'
      endif

      RETURN
      END
