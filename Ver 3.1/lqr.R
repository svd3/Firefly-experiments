## LQR implementation

lqr = function(sys){
	A = sys$A
	B = sys$B
	R = sys$R
	Q = sys$Q
	C0 = sys$C0 # process noise
	D0 = sys$D0 # measurement noise
	H = sys$H # measurement matrix
	
	nx = dim(Q)[1]
	ny = dim(H)[1]
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
	
	P = array(0,c(nx,nx,N))
	K = array(0,c(nx,ny,N-1))
	P[,,1] = sys$S1
	## find kalman gain forward pass
	for(k in 1:(N-1)){
	   K[,,k] = A[,,k]%*%P[,,k]%*%t(H[,,k])%*%pinv(H[,,k]%*%P[,,k]%*%t(H[,,k]) + D0[,,k]%*%t(D0[,,k]))
	   P[,,k+1] = (A[,,k] - K[,,k]%*%H[,,k])%*%P[,,k]%*%t(A[,,k]) + C0[,,k]%*%t(C0[,,k])
	}
	
	## u = - L * x
	list(L =L, K =K)
}
