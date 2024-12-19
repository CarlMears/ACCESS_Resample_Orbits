	SUBROUTINE moon_vector(days, moonVec,moonDis)	 !returns vector in J2000 system, dis is km
	IMPLICIT NONE

	REAL(8), INTENT(IN) :: days
	REAL(8), INTENT(OUT):: moonVec(3),moonDis  
	REAL(8) T,L_0, l, l_p, F, D, L_M, B_M, dL, twiceD, twicel, l_2D, l_p_2d,l_plus_2D, twicel_twiceD
	REAL(8)  cosL_M, sinL_M, cosB_M, sinB_M, eps,coseps, sineps 


      T=days/36525.D0

 	eps=23.4393D0 - 0.01300D0*T
	coseps=cosd(eps)
	sineps=sind(eps)

	L_0 = 218.31617d0 + 481266.48368d0*T	 !degrees
	l   = 134.96292d0 + 477198.86753d0*T
	l_p = 357.52543d0 +  35999.04944d0*T
	F   =  93.27283d0 + 483202.01873d0*T
	D   = 297.85027d0 + 445267.11135d0*T

	twiceD = 2*D
	twicel = 2*l
	l_2D = l - twiceD
	l_p_2d = l_p - twiceD
	l_plus_2D = l + twiceD
	twicel_twiceD = twicel - twiceD
	
	dL = 22640*dsind(l)        + 769*dsind(twicel)   - 4586*dsind(l_2D)     + 2370*dsind(twiceD)
     &      -668*dsind(l_p)      - 412*dsind(2*F)       - 212*dsind(twicel_twiceD) 
     &      -206*dsind(l+l_p_2D) + 192*dsind(l_plus_2D) - 165*dsind(l_p_2D) 
     &      +148*dsind(l-l_p)    - 125*dsind(D)         - 110*dsind(l+l_p) 
     &       -55*dsind(2*F-twiceD)   !arcseconds
		
	L_M = L_0 +  dL/3600.d0

	B_M = 18520*dsind(F + (dL +  412*dsind(2*F) +     541*dsind(l_p))/3600.d0)
     &       -526*dsind(F-twiceD) + 44*dsind(l_2D+F )  - 31*dsind(F-l_plus_2D)
     &        -25*dsind(F-twicel) - 23*dsind(l_p_2D+F) + 21*dsind(F-l) 
     &       + 11*dsind(F-l_p-twiceD)  !arcseconds

	B_M = B_M/3600.d0	 !degrees

	moonDis= 385000 - 20905*dcosd(l)         - 3699*dcosd(l_2D)         - 2956*dcosd(twiceD)
     &                   -570*dcosd(twicel)    +  246*dcosd(twicel_twiceD) - 205*dcosd(l_p_2D) 
     &                   -171*dcosd(l_plus_2D) -  152*dcosd(l + l_p_2D)  	 !km
	
	cosL_M = dcosd(L_M)
	sinL_M = dsind(L_M)
	cosB_M = dcosd(B_M)
	sinB_M = dsind(B_M)
		
	moonVec(1) = cosL_M * cosB_M
	moonVec(2) = coseps * sinL_M * cosB_M - sineps * sinB_M
	moonVec(3) = sineps * sinL_M * cosB_M + coseps * sinB_M     
	
	END SUBROUTINE moon_vector		
