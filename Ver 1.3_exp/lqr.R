## LQR implementation

lqr = function(sys){
	A = sys$A
	B = sys$B
	R = sys$R
	Q = sys$Q
	
	nx = dim(Q)[1]
	nu = dim(R)[1]
	N = dim(Q)[3]

	P = array(0,c(nx,nx,N))
	L = array(0,c(nu,nx,N-1))
	P[,,N] = Q[,,N]
	
	for(k in (N-1):1){
		Bt.P = t(B[,,k]) %*% P[,,k+1]
		L[,,k] = pinv(R[,,k] + (Bt.P %*% B[,,k])) %*% Bt.P %*% A[,,k]
		
		A_BL = A[,,k] - (B[,,k] %*% L[,,k])
		P[,,k] = Q[,,k] + (t(L[,,k]) %*% R[,,k] %*% L[,,k]) + (t(A_BL) %*% P[,,k+1] %*% A_BL)
	}
	## u = - L * x
	L
}
