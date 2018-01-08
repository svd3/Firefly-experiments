## kalman lqg my implementation
library(pracma)
tr = function(A){
    sum(diag(A))
}

kalman_lqg = function(sys, maxIter = 100, eps = 1e-15){
    Init=1
    #print("my kalman_lqg :: ")
    As = sys$A
    Bs = sys$B
    Cs = sys$C
    C0s = sys$C0
    Hs = sys$H
    Ds = sys$D
    D0s = sys$D0
    E0s = sys$E0
    Q = sys$Q
    Rs = sys$R
    X1 = sys$X1
    S1 = sys$S1
    
    #maxIter = 500;
    #eps = 1e-15;
    
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
    
    K = array(0, c(szX,szY,N-1));
    L = array(0, c(szU,szX,N-1));
    
    cost = array(0, maxIter)
    prevCost = 10; # some random big value for first loop
    
    for(iter in 1:maxIter){
        #print(iter)
        # initialize covariances
        SiE = S1;
        SiX = X1%*%t(X1);
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
            E0 = E0s#[,,k] # what all can be time varying ??
            
            # compute Kalman gain
            temp = SiE + SiX + SiXE + t(SiXE);
            
            DSiD = array(0, c(szY,szY));
            for(i in 1:szD){
                DSiD = DSiD + D[,,i]%*%temp%*%t(D[,,i]);
                #DSiD = DSiD + D%*%temp%*%t(D);
            }
            
            K[,,k] = A%*%SiE%*%t(H)%*%ginv(H%*%SiE%*%t(H)+D0%*%t(D0)+DSiD);
            
            # compute new SiE
            newE =  E0%*%t(E0) + C0%*%t(C0) + (A-K[,,k]%*%H)%*%SiE%*%t(A);
            LSiL = L[,,k]%*%SiX%*%t(L[,,k]);
            for(i in 1:szC){
                #newE = newE + B%*%C[,,i]%*%LSiL%*%t(C[,,i])%*%t(B);
                newE = newE + C[,,i]%*%LSiL%*%t(C[,,i]);
            }
            
            # update SiX, SiE, SiXE
            SiX = E0%*%t(E0) + K[,,k]%*%H%*%SiE%*%t(A) + (A-B%*%L[,,k])%*%SiX%*%t(A-B%*%L[,,k]) + (A-B%*%L[,,k])%*%SiXE%*%t(H)%*%t(K[,,k]) + K[,,k]%*%H%*%t(SiXE)%*%t(A-B%*%L[,,k]);
            SiE = newE;
            SiXE = (A-B%*%L[,,k])%*%SiXE%*%t(A-K[,,k]%*%H) - E0%*%t(E0);
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
        cost[iter] = 0;
        
        # backward pass - recompute control policy
        for(k in (N-1):1){
            A = As[,,k]
            B = Bs[,,k]
            C0 = C0s[,,k]
            C = Cs[,,,k]
            H = Hs[,,k]
            D0 = D0s[,,k]
            D = Ds[,,,k]
            E0 = E0s#[,,k]
            R = Rs[,,k]
            
            # update Cost
            cost[iter] = cost[iter] + tr(Sx%*%C0%*%t(C0)) + tr(Se%*%(K[,,k]%*%D0%*%t(D0)%*%t(K[,,k]) + E0%*%t(E0) + C0%*%t(C0)));
            
            # Controller
            temp = R + t(B)%*%Sx%*%B;
            #BSxeB = t(B)%*%(Sx+Se)%*%B;
            BSxeB = (Sx+Se);
            for(i in 1:size(C,3)){
                temp = temp + t(C[,,i])%*%BSxeB%*%C[,,i];
            }
            L[,,k] = ginv(temp)%*%t(B)%*%Sx%*%A;
            
            # compute new Se
            newE = t(A)%*%Sx%*%B%*%L[,,k] + t(A-K[,,k]%*%H)%*%Se%*%(A-K[,,k]%*%H);
            
            # update Sx and Se
            Sx = Q[,,k] + t(A)%*%Sx%*%(A-B%*%L[,,k]);
            KSeK = t(K[,,k])%*%Se%*%K[,,k];
            for(i in 1:szD){
                Sx = Sx + t(D[,,i])%*%KSeK%*%D[,,i];
            }
            Se = newE;
        }
        
        # adjust cost
        cost[iter] = cost[iter] + t(X1)%*%Sx%*%X1 + tr((Se+Sx)%*%S1);
        
        # check convergence of Cost
        # check here
        if(iter>1 & abs(prevCost-cost[iter])<eps){
            break;
        }
        prevCost = cost[iter]
    }
    list(K=K, L=L, cost=cost, iters=iter)
}