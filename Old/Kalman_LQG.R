## hi KaLman LQG
library(pracma)
tr = function(A){
	sum(diag(A))
}

kalman_lqg = function(systm){
    #print("kalman_lqg :: ")
    A = systm$A
    B = systm$B
    C = systm$C
    C0 = systm$C0
    H = systm$H
    D = systm$D
    D0 = systm$D0
    E0 = systm$E0
    Q = systm$Q
    R = systm$R
    X1 = systm$X1
    S1 = systm$S1
    ans = kalman(A,B,C,C0, H,D,D0, E0, Q,R, X1,S1, NSim =0,Init=1,Niter=0)
    ans
}

kalman = function(As,Bs,Cs,C0s, Hs,Ds,D0s, E0s, Q,Rs, X1,S1, NSim,Init,Niter){
    #################
    #print("kalman :: ")# Initialization
    
    # numerical parameters
    MaxIter = 500;
    Eps = 1e-15;
    
    # determine sizes
    szX = size(As,1);
    szU = size(Bs,2);
    szY = size(Hs,1);
    szC = size(Cs,3);
    szC0 = size(C0s,2);
    szD = size(Ds,3);
    szD0 = size(D0s,2);
    szE0 = size(E0s,2);
    N = size(Q,3);
    
    # initialize missing optional parameters
    # initialize policy and filter
    K = array(0, c(szX,szY,N-1));
    L = array(0, c(szU,szX,N-1));
    
    #######################################################################
    # run iterative algorithm - until convergence or MaxIter
    Cost = array(0, MaxIter)
    prevCost = 10; # some random big value for first loop
    for(iter in 1:MaxIter){
    	 # initialize covariances
    	 SiE = S1;
         #SiX = X1%*%t(X1);
         temp_X = as.vector(t(X1))
         SiX = outer(X1,X1)
    	 SiXE = array(0, c(szX,szX));
    	 
    	 # forward pass - recompute Kalman filter  
    	 for(k in 1:(N-1)){
             # adapt this loop for time-varying systems
             A = As[,,k]
             B = Bs[,,k]
             C0 = C0s[,,k]
             C = Cs[,,,k]
             H = Hs[,,k]
             D0 = D0s[,,k]
             D = Ds[,,,k]
             E0 = E0s[,,k]
             # compute Kalman gain
             temp = SiE + SiX + SiXE + t(SiXE);
             if(size(D,2)==1){
                 DSiD = diag(diag(temp)*D^2);
             } else{
                 DSiD = array(0, c(szY,szY));
                 for(i in 1:szD){
                     DSiD = DSiD + D[,,i]%*%temp%*%t(D[,,i]);
                 }
             }
             K[,,k] = A%*%SiE%*%t(H)%*%ginv(H%*%SiE%*%t(H)+D0%*%t(D0)+DSiD);
             
             # compute new SiE
             E0.tE0 = as.numeric(E0%*%t(E0))
             newE =  E0.tE0 + C0%*%t(C0) + (A-K[,,k]%*%H)%*%SiE%*%t(A);
             LSiL = L[,,k]%*%SiX%*%t(L[,,k]);
             if(size(C,2)==1){
                 newE = newE + B%*%diag(diag(LSiL)*C^2)%*%t(B);
             } else{
                 for(i in 1:szC){
                     newE = newE + B%*%C[,,i]%*%LSiL%*%t(C[,,i])%*%t(B);
                 }
             }
             
             # update SiX, SiE, SiXE
             E0.tE0 = as.numeric(E0%*%t(E0))
             SiX = E0.tE0 + K[,,k]%*%H%*%SiE%*%t(A) + (A-B%*%L[,,k])%*%SiX%*%t(A-B%*%L[,,k]) + (A-B%*%L[,,k])%*%SiXE%*%t(H)%*%t(K[,,k]) + K[,,k]%*%H%*%t(SiXE)%*%t(A-B%*%L[,,k]);
             SiE = newE;
             SiXE = (A-B%*%L[,,k])%*%SiXE%*%t(A-K[,,k]%*%H) - E0.tE0;
         }
    	 
    	 # first pass initialization
    	 if(iter==1){
    	 	if(Init ==0){	# open loop
    	 		K = array(0, c(szX,szY,N-1));
    	 	} else if(Init ==2){	#random
                K = array(rnorm(szX*szY*(N-1)), c(szX,szY,N-1));
    	 	}
    	 }
    	 
    	 # initialize optimal cost-to-go function
    	 Sx = Q[,,N];
    	 Se = array(0, c(szX,szX));
    	 Cost[iter] = 0;
         
    	 # backward pass - recompute control policy
    	 for(k in (N-1):1){
             A = As[,,k]
             B = Bs[,,k]
             C0 = C0s[,,k]
             C = Cs[,,,k]
             H = Hs[,,k]
             D0 = D0s[,,k]
             D = Ds[,,,k]
             E0 = E0s[,,k]
             R = Rs[,,k]
             
             # update Cost
             E0.tE0 = as.numeric(E0%*%t(E0))
             Cost[iter] = Cost[iter] + tr(Sx%*%C0%*%t(C0)) + tr(Se%*%(K[,,k]%*%D0%*%t(D0)%*%t(K[,,k]) + E0.tE0 + C0%*%t(C0)));
             
             # Controller
             temp = R + t(B)%*%Sx%*%B;
             BSxeB = t(B)%*%(Sx+Se)%*%B;
             if(size(C,2)==1){
                 temp = temp + diag(diag(BSxeB)*C^2);
             } else{
                 for(i in 1:size(C,3)){
                     temp = temp + t(C[,,i])%*%BSxeB%*%C[,,i];
                 }
             }
             L[,,k] = ginv(temp)%*%t(B)%*%Sx%*%A;
             
             # compute new Se
             newE = t(A)%*%Sx%*%B%*%L[,,k] + t(A-K[,,k]%*%H)%*%Se%*%(A-K[,,k]%*%H);
             
             # update Sx and Se
             Sx = Q[,,k] + t(A)%*%Sx%*%(A-B%*%L[,,k]);
             KSeK = t(K[,,k])%*%Se%*%K[,,k];
             if(size(D,2)==1){
                 Sx = Sx + diag(diag(KSeK)*D^2);
             } else{
                 for(i in 1:szD){
                     Sx = Sx + t(D[,,i])%*%KSeK%*%D[,,i];
                 }
             }
             Se = newE;
         }
    	 
    	 # adjust cost
    	 Cost[iter] = Cost[iter] + t(X1)%*%Sx%*%X1 + tr((Se+Sx)%*%S1);
    	 
    	 # check convergence of Cost
         # check here
    	 if((Niter>1 & iter>=Niter) | (Niter==1 & iter>1 & abs(prevCost-Cost[iter])<Eps)
    	 		| (Niter==1 & iter>20 & sum(diff(dist(iter-10:iter))>0)>3)){
    	 			break;
    	 }
    	 prevCost = Cost[iter]
    }
    
    # compute average trajectory
    Xa = array(0, c(szX,N))
    Xa[,1] = X1
    for(k in 1:(N-1)){
        u = -L[,,k]%*%Xa[,k]
        Xa[,k+1] = A%*%Xa[,k] + B%*%u
    }
    
    # simulate noisy trajectories
    Xsim = 0; CostSim = 0
    #print(size(L))
    list(K=K, L=L, Cost=Cost, Xa=Xa, Xsim = Xsim, CostSim=CostSim, iters=iter)
}
