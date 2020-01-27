      INTEGER  FUNCTION MINDAT(LUN)
      implicit none
*     
      INTEGER LUN
      
      include ?
*     
      MINDAT = 0
*
*      include 'cuts.f'
      if (esclre2.lt.4.0)  return
      if (esclre1.lt.10.0) return
      if (esclrr1.gt.3.2)  return
      if (esclrr2.gt.3.2)  return
      if (ESCLRRA1.lt.9.5) return
      if (ESCLRRA2.lt.9.5) return
      if (elar60.gt.1.0)   return
      if (esclre1+esclre2.lt.20) return
      if (esclesum-esclre1-esclre2.gt.0.5) return
*     
      MINDAT = 1
*
      write(lun,*) ESCLRX1-VTXX,ESCLRY1-VTXY,ESCLRZ1-VTXZ-1.0,
     &             ESCLRX2-VTXX,ESCLRY2-VTXY,ESCLRZ2-VTXZ-1.0,
     &             VTSX,VTSY
      RETURN
      END
      

      
