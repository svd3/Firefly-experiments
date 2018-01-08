library(pracma)
tr = function(A){
    sum(diag(A))
}

lqg_all = function(sys, maxIter = 100, eps = 1e-15){
    #dynamic matrices
    As = sys$A
    Bs = sys$B
    Hs = sys$H
    Es = sys$E
    
    #noise matrices
    Cxs = sys$Cx
    Cus = sys$Cu
    C0s = sys$C0
    Dxs = sys$Dx
    Dus = sys$Du
    D0s = sys$D0
    X1 = sys$X1
    S1 = sys$S1
    
    #cost matrices
    Qs = sys$Q
    Rs = sys$R
    qs = sys$q
    rs = sys$r
    Ps = sys$P
    
    nx = size(As,1)
    nu = size(Bs,2)
    ny = size(Hs,1)
    nw = size(C0s,2)
    nv = size(D0s,2)
    N = size(Qs,3)
    
    L = array(0, c(nu,nx,N-1))
    l = array(0, c(nu,N-1))
    K = array(0, c(nx,ny,N-1))
    cost = 0
    for(iter in 1:maxIter){
        ######### Kalman Filter
        
        SiX = X1%*%t(X1)
        SiE = S1
        SiXE = array(0, c(nx,nx))
        mx = X1; me = array(0,nx)
        for(k in 1:(N-1)){
            A = As[,,k]
            B = Bs[,,k]
            H = Hs[,,k]
            E = Es[,,k]
            
            C0 = C0s[,,k]
            Cx = Cxs[,,,k]
            Cu = Cus[,,,k]
            D0 = D0s[,,k]
            Dx = Dxs[,,,k]
            Du = Dus[,,,k]
            
            Pk = array(0, c(ny,ny))
            Mk = array(0, c(nx,nx))
            
            mxe = mx + me  ## mx + me needed below!!
            for(i in 1:nw){
                Mk = Mk + C0[,i]%*%t(C0[,i])
                temp = Cx[,,i]%*%mxe%*%t(C0[,i]) + 
                    Cu[,,i]%*%(l[,k]+L[,,k]%*%mx)%*%t(C0[,i]) + 
                    Cu[,,i]%*%(l[,k]%*%t(mxe) + L[,,k]%*%(SiX + SiXE))%*%t(Cx[,,i])
                Mk = Mk + temp + t(temp)
                Mk = Mk + Cx[,,i]%*%(SiX + SiE + SiXE + t(SiXE))%*%t(Cx[,,i]) +
                    Cu[,,i]%*%(l[,k]%*%t(l[,k] + L[,,k]%*%mx) +
                                   L[,,k]%*%(mx%*%t(l[,k]) + SiX%*%t(L[,,k])))%*%t(Cu[,,i])
            }
            for(i in 1:nv){
                Pk = Pk + D0[,i]%*%t(D0[,i])
                temp = Dx[,,i]%*%mxe%*%t(D0[,i]) + 
                    Du[,,i]%*%(l[,k]+L[,,k]%*%mx)%*%t(D0[,i]) + 
                    Du[,,i]%*%(l[,k]%*%t(mxe) + L[,,k]%*%(SiX + SiXE))%*%t(Dx[,,i])
                Pk = Pk + temp + t(temp)
                Pk = Pk + Dx[,,i]%*%(SiX + SiE + SiXE + t(SiXE))%*%t(Dx[,,i]) +
                    Du[,,i]%*%(l[,k]%*%t(l[,k] + L[,,k]%*%mx) +
                                   L[,,k]%*%(mx%*%t(l[,k]) + SiX%*%t(L[,,k])))%*%t(Du[,,i])
            }
            
            K[,,k] = A%*%SiE%*%t(H)%*%pinv(H%*%SiE%*%t(H) + Pk)
            
            KH = K[,,k]%*%H; A_BL = A + B%*%L[,,k]; Bl = B%*%l[,k]
            temp = A_BL%*%SiXE%*%t(KH) + B%*%l[,k]%*%t(A_BL%*%mx + KH%*%me)
            SiX = A_BL%*%SiX%*%t(A_BL) + KH%*%SiE%*%t(A) +
                temp + t(temp) + Bl%*%t(Bl)
            
            SiE = (A-KH)%*%SiE%*%t(A) + Mk
            SiXE = A_BL%*%SiXE%*%t(A-KH) + Bl%*%t((A-KH)%*%me)
            
            mx = A_BL%*%mx + KH%*%me + Bl
            me = (A-KH)%*%me
        }
        
        ############## Control Gain
        Q = Qs
        q = qs
        Sx = Q[,,N]
        Sxh = Sxxh = array(0,c(nx,nx)) #xh = x^ or xhat
        sx = qs[,N]
        sxh = array(0,nx)
        
        J = array(0, c(nu,nu,N-1))
        g = array(0, c(nu,N-1)) # think g = rk
        G = Gx = Gxh = array(0, c(nu,nx,N-1)) ## dim of dl/du.dx = P
    
        for(k in (N-1):1){
            A = As[,,k]
            B = Bs[,,k]
            H = Hs[,,k]
            E = Es[,,k]
            
            C0 = C0s[,,k]
            Cx = Cxs[,,,k]
            Cu = Cus[,,,k]
            D0 = D0s[,,k]
            Dx = Dxs[,,,k]
            Du = Dus[,,,k]
            
            R = Rs[,,k]
            #Q = Qs ## check
            #q = qs
            r = rs[,k]
            P = Ps[,,k]
            
            KH = K[,,k]%*%H  ## important
            
            J[,,k] = R + t(B)%*%(Sx + Sxh + 2*Sxxh)%*%B
            g[,k] = r + t(B)%*%(sx + sxh)
            Gx[,,k] = P + t(B)%*%(Sx + Sxxh)%*%A + t(B)%*%(Sxh + Sxxh)%*%KH
            Gxh[,,k] = t(B)%*%(Sxh + Sxxh)%*%(A - KH)    
            for(i in 1:nw){
                CuSx = t(Cu[,,i])%*%Sx
                J[,,k] = J[,,k] + CuSx%*%Cu[,,i]
                g[,k] = g[,k] + CuSx%*%C0[,i]
                Gx[,,k] = Gx[,,k] + CuSx%*%Cx[,,i]
            }
            for(i in 1:nv){
                DKSK = t(K[,,k]%*%Du[,,i])%*%Sxh%*%K[,,k]
                J[,,k] = J[,,k] + DKSK%*%Du[,,i]
                g[,k] = g[,k] + DKSK%*%D0[,i]
                Gx[,,k] = Gx[,,k] + DKSK%*%Dx[,,i]
            }
            G[,,k] = Gx[,,k] + Gxh[,,k]
            
            ## set J to positive semidefinite??
            invJ = pinv(J[,,k]); L[,,k] = -invJ%*%G[,,k]
            l[,k] = -invJ%*%g[,k] ## ??????????????????????????????
            #e = eigen(J[,,k])
            #print(e$values)
            
            newSx = Q[,,k] + t(A)%*%Sx%*%A + t(KH)%*%Sxh%*%KH + 2*t(A)%*%Sxxh%*%KH
            
            temp = t(L[,,k])%*%Gxh[,,k] ## need to find L and l first
            newSxh = t(A-KH)%*%Sxh%*%(A-KH) + t(L[,,k])%*%J[,,k]%*%L[,,k] + temp + t(temp)
            
            newSxxh = (t(KH)%*%Sxh + t(A)%*%Sxxh)%*%(A-KH) + t(Gx[,,k])%*%L[,,k]
            
            new_sx = q[,k] + t(A)%*%sx + t(KH)%*%sxh + t(Gx[,,k])%*%l[,k]
            
            new_sxh = t(A-KH)%*%sxh + t(L[,,k])%*%J[,,k]%*%l[,k] + 
                t(L[,,k])%*%g[,k] + t(Gxh[,,k])%*%l[,k]
            
            for(i in 1:nw){
                CxSx = t(Cx[,,i])%*%Sx
                newSx = newSx + CxSx%*%Cx[,,i]
                new_sx = new_sx + CxSx%*%C0[,i]
            }
            for(i in 1:nv){
                DKSK = t(K[,,k]%*%Dx[,,i])%*%Sxh%*%K[,,k]
                newSx = newSx + DKSK%*%Dx[,,i]
                new_sx = new_sx + DKSK%*%D0[,i]
            }
            
            Sx= newSx
            Sxh = newSxh
            Sxxh = newSxxh
            sx = new_sx
            sxh = new_sxh
        }
        
        cost[iter+1] = tr(Sx + Sxh + 2*Sxxh)
        
        if(iter>1 & abs(cost[iter+1]-cost[iter])<eps){
            break;
        }
    }
    list(K=K, L=L, l=l, cost=cost, iters=iter)
}