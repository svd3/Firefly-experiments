## Notations:
## x = actual trajectory;   u = actual control
## x_ = nominal trajectory; u_ = nominal control
## dx = deviation from nominal i.e. "dx = x - x_"
## du = deviation from nominal i.e. "du = u - u_"
## dy = deviation from nominal i.e. "dy = y - y_"
## dt = small time step


## make dimensions global

## Libraries:
library(pracma)
library(MASS)


linearize = function(x_, u_, f, g, h, l, Fn, G){
    ## load dimensions
    nx = dim(x_)[1]
    nu = dim(u_)[1]
    N = dim(x_)[2]
    ## add ny
    ## ny as of now global
    
    ## dx(k+1) = A(k)*dx(k) + B(k)*du(k)
    ## dy(k) = C(k)*dx(k) + D(k)*du(k)
    ## A(k) = I + {df/dx @x_(k)}*dt  note: df/dx (partial derivative) not deviation!!
    ## B(k) = {df/du @u_(k)}*dt
    
    # state augmentation
    dxa = c(rep(0,(nx+nu)),1) ## (x, u, 1)
    na = nx + nu + 1
    sys = list()
    sys$A = array(0,c(na,na,N-1))
    sys$B = array(0,c(na,nu,N-1))
    sys$X1 = c(dxa)
    sys$S1 = 0*diag(na)
    sys$S1[na,na] = 0
    sys$H = array(0,c(ny,na,N-1))
    sys$Q = array(0,c(na,na,N))
    sys$R = array(0,c(nu,nu,N-1))
    
    # system noises
    sys$C0 = array(0,c(na,nw,N-1));  C0 = array(0,c(na,nw))
    sys$D0 = array(0,c(ny,nv,N-1)); D0 = array(0,c(ny,nv))
    sys$C = array(0,c(na,nu,nw,N-1)); C = array(0,c(na,nu,nw))
    sys$D = array(0,c(ny,na,nv,N-1)); D = array(0,c(ny,na,nv))
    sys$E0 = array(0,c(na,nw))
    ## na is the new nx
        
    ## linearize (A_k, B_k, C_k, D_k)
    for(k in 1:(N-1)){
        f_x = jacobian(f, x0 = x_[,k], u = u_[,k]);
        A = diag(nx) + f_x*dt
        
        f_u = jacobian(f, x0 = u_[,k], x = x_[,k]);
        B = f_u*dt
        
        H = jacobian(g, x0 = x_[,k], u = u_[,k]); #dg/dx
        H = diag(nx)
        Q = dt*hessian(l, x0 = x_[,k], u = u_[,k]) ## d^2l/dx^2
        q = dt*jacobian(l, x0 = x_[,k], u = u_[,k])
        s = dt*l(x_[,k],u_[,k])
        
        if(k==1){
            r = array(0,nu)
            R = array(0,c(nu,nu))
        } else{
            r = dt*jacobian(l, x0 = u_[,k-1], x = x_[,k-1])
            R = dt*hessian(l, x0 = u_[,k-1], x = x_[,k-1]) ## d^2l/du^2
        }
        
        #C0 = array(0, c(na,na)); C0[1:nx,1:nx] = 0.01*diag(nx)
        
        # noises
        C0[1:nx,] = sqrt(dt)*Fn(x_[,k], u_[,k])
        Fn_u = jacobian(Fn, x0 = u_[,k], x = x_[,k])
        for(i in 1:nw){
            m1 = (i-1)*nx + 1; m2 = i*nx
            C[1:nx,,i] = Fn_u[m1:m2,]
        }
        C = sqrt(dt)*C
        
        D0 = G(x_[,k], u_[,k])/sqrt(dt) 
        G_x = jacobian(G, x0 = x_[,k], u = u_[,k])
        for(j in 1:nv){
            m1 = (j-1)*ny + 1; m2 = j*ny
            D[,1:nx,j] = G_x[m1:m2,]
        }
        D = D/sqrt(dt)
        
        sys$A[1:nx,1:nx,k] = A; sys$A[na,na,k] = 1
        sys$B[1:nx,1:nu,k] = B; sys$B[(nx+1):(nx+nu),1:nu,k] = diag(nu)
        sys$Q[1:nx,1:nx,k] = Q;
        sys$Q[na,1:nx,k] = q/2; sys$Q[1:nx,na,k] = q/2
        sys$Q[(nx+1):(nx+nu),(nx+1):(nx+nu),k] = R;
        sys$Q[na,(nx+1):(nx+nu),k] = r/2; sys$Q[(nx+1):(nx+nu),na,k] = r/2
        sys$Q[na,na,k] = s
        sys$H[,1:nx,k] = H
        sys$C0[,,k] = C0
        sys$D0[,,k] = D0
        sys$C[,,,k] = C
        sys$D[,,,k] = D
    }
    k = N
    Q = hessian(h, x0 = x_[,N]);
    q = jacobian(h, x0 = x_[,N])
    R = dt*hessian(l, x0 = u_[,N-1], x = x_[,N-1]);
    r = dt*jacobian(l, x0 = u_[,N-1], x = x_[,N-1])
    s = h(x_[,N])
    
    sys$Q[1:nx,1:nx,k] = Q;
    sys$Q[na,1:nx,k] = q/2; sys$Q[1:nx,na,k] = q/2
    sys$Q[(nx+1):(nx+nu),(nx+1):(nx+nu),k] = R;
    sys$Q[na,(nx+1):(nx+nu),k] = r/2; sys$Q[(nx+1):(nx+nu),na,k] = r/2
    sys$Q[na,na,k] = s
    #sys$Q[1:nx,1:nx,N] = Q; sys$Q[nx+1,1:nx,N] = q/2; sys$Q[1:nx,nx+1,N] = q/2; sys$Q[nx+1,nx+1,N] = h(x_[,N])
    #sys$Q[,,N] = 2*sys$Q[,,N]
    
    # noises
    #sys$C0 = C0 
    #sys$C = 0
    #sys$D0 = 0.01*diag(ny)
    #sys$D = 0
    #sys$E0 = 0
    sys
}
