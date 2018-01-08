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


linearize2 = function(x_, u_, f, g, h, l, Fn, G){
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

    sys = list()
    #dynamics matrices
    sys$A = array(0,c(nx,nx,N-1))
    sys$B = array(0,c(nx,nu,N-1))
    sys$H = array(0,c(ny,nx,N-1))
    sys$E = array(0,c(ny,nu,N-1))
    
    sys$X1 = rep(0,nx)
    sys$S1 = 0*diag(nx)
    
    #cost matrices
    sys$Q = array(0,c(nx,nx,N))
    sys$R = array(0,c(nu,nu,N-1))
    sys$q = array(0,c(nx,N))
    sys$r = array(0,c(nu,N-1))
    sys$P= array(0,c(nu,nx,N-1))
    
    #noise matrices
    sys$C0 = array(0,c(nx,nw,N-1));  C0 = array(0,c(nx,nw))
    sys$Cx = array(0,c(nx,nx,nw,N-1));  Cx = array(0,c(nx,nx,nw))
    sys$Cu = array(0,c(nx,nu,nw,N-1));  Cu = array(0,c(nx,nu,nw))
    
    sys$D0 = array(0,c(ny,nv,N-1)); D0 = array(0,c(ny,nv))
    sys$Dx = array(0,c(ny,nx,nv,N-1));  Dx = array(0,c(ny,nx,nv))
    sys$Du = array(0,c(ny,nu,nv,N-1));  Du = array(0,c(ny,nu,nv))

        
    ## linearize (A_k, B_k, C_k, D_k)
    for(k in 1:(N-1)){
        f_x = jacobian(f, x0 = x_[,k], u = u_[,k]); #df/dx
        A = diag(nx) + f_x*dt
        
        f_u = jacobian(f, x0 = u_[,k], x = x_[,k]); #df/du
        B = f_u*dt
        
        H = jacobian(g, x0 = x_[,k], u = u_[,k]); #dg/dx
        #H = diag(nx)
        E = jacobian(g, x0 = u_[,k], x = x_[,k]); #dg/du
        
        Q = dt*hessian(l, x0 = x_[,k], u = u_[,k]) ## d^2l/dx^2
        q = dt*jacobian(l, x0 = x_[,k], u = u_[,k])
        s = dt*l(x_[,k],u_[,k])
        r = dt*jacobian(l, x0 = u_[,k], x = x_[,k])
        R = dt*hessian(l, x0 = u_[,k], x = x_[,k]) ## d^2l/du^2
        P = array(0,c(nu,nx)) # zero for now, we'll update later
        
        #C0 = array(0, c(na,na)); C0[1:nx,1:nx] = 0.01*diag(nx)
        
        # noises
        C0 = sqrt(dt)*Fn(x_[,k], u_[,k]) ## Fn is nx x nw
        Fn_u = jacobian(Fn, x0 = u_[,k], x = x_[,k])
        Fn_x = jacobian(Fn, x0 = x_[,k], u = u_[,k])
        for(i in 1:nw){
            m1 = (i-1)*nx + 1; m2 = i*nx
            Cu[,,i] = Fn_u[m1:m2,]
            Cx[,,i] = Fn_x[m1:m2,]
        }
        Cu = sqrt(dt)*Cu
        Cx = sqrt(dt)*Cx
        
        D0 = G(x_[,k], u_[,k])/sqrt(dt) 
        G_x = jacobian(G, x0 = x_[,k], u = u_[,k])
        G_u = jacobian(G, x0 = u_[,k], x = x_[,k])
        for(j in 1:nv){
            m1 = (j-1)*ny + 1; m2 = j*ny
            Dx[,,j] = G_x[m1:m2,]
            Du[,,j] = G_u[m1:m2,]
        }
        Dx = Dx/sqrt(dt)
        Du = Du/sqrt(dt)
        
        sys$A[,,k] = A;
        sys$B[,,k] = B;
        sys$H[,,k] = H;
        sys$E[,,k] = E;
        
        sys$Q[,,k] = Q;
        sys$R[,,k] = R;
        #sys$P[,,k] = P;
        sys$q[,k] = q;
        sys$r[,k] = r;

        sys$C0[,,k] = C0
        sys$Cx[,,,k] = Cx
        sys$Cu[,,,k] = Cu
        sys$D0[,,k] = D0
        sys$Dx[,,,k] = Dx
        sys$Du[,,,k] = Du
    }
    k = N
    Q = hessian(h, x0 = x_[,N]);
    q = jacobian(h, x0 = x_[,N])
    sys$Q[,,N] = Q;
    sys$q[,N] = q;
    sys
}
