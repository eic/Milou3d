	program convol
	implicit none
	real*8 lamlo,lamnlo
	integer gpd, pol,nx,ndel,nq,nfl,tw,ord,k

C
C       setting initial twist
C	
	tw = 2

C
C       nx = # of points in X, nq = # of points in Q^2, 
C       ndel = # of points in zeta (x_bj)
C
	
	nx = 499
	ndel = 58
	nq = 15

C
C       do both twist-2 (k=2) and twist-3 (k=3),
C

	do 1 k = 2,3

	   tw = k
C
C       which GPD: H,\tilde H (gpd=1) or E, \tilde E (gpd=2)
C

	   do 2 gpd = 1,2
C
C       order: LO (ord=1), NLO (ord=2)
C

	      do 3 ord = 1,2

		 if (ord.eq.2.and.tw.eq.3) goto 3
C
C unpolarized (pol=1) and/or polarized (pol=2)
C		 

		 do 4 pol = 1,2
	  
C
C       Lambda_QCD for LO/NLO and unpolarized or polarized
C


	if (pol.eq.1) then
	lamlo = 0.220D0
	lamnlo = 0.323D0
	else
	  lamlo = 0.174D0
	  lamnlo = 0.246D0
	endif 

	if (tw.eq.2.and.pol.eq.1) then

	   nfl = 4

	   else

	   nfl = 3

	endif
 
	call convolution(tw,gpd,pol,ord,lamlo,lamnlo,nfl,nx,ndel,nq)

 4	enddo
 3	enddo
 2	enddo
 1	enddo

	end
