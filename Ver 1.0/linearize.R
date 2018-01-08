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


linearize = function(x_,u_,f,g,h,l){
    ## load dimensions
    
    ## dx(k+1) = A(k)*dx(k) + B(k)*du(k)
    ## dy(k) = C(k)*dx(k) + D(k)*du(k)
    ## A(k) = I + {df/dx @x_(k)}*dt  note: df/dx (partial derivative) not deviation!!
    ## B(k) = {df/du @u_(k)}*dt
    
    # state augmentation
    dxa = c(rep(0,nx),1)
    sys = list()
    sys$A = array(0,c((nx+1),(nx+1),N-1))
    sys$B = array(0,c((nx+1),(nu+1),N-1))
    sys$X1 = c(dxa)
    sys$S1 = diag(nx+1)
    sys$S1[nx+1,nx+1] = 0
    sys$H = array(0,c(ny,(nx+1),N-1))
    sys$Q = array(0,c((nx+1),(nx+1),N))
    sys$R = array(0,c((nu+1),(nu+1),N-1))
    
    ## linearize (A_k, B_k, C_k, D_k)
    for(k in 1:(N-1)){
        f_x = jacobian(f, x0 = x_[,k], u = u_[,k]); dim(f_x) = c(nx,nx)
        A = diag(nx) + f_x*dt
        
        f_u = jacobian(f, x0 = u_[,k], x = x_[,k]); dim(f_u) = c(nx,nu)
        B = f_u*dt
        
        H = jacobian(g, x0 = x_[,k], u = u_[,k]); dim(H) = c(ny,nx) #dg/dx
        H = diag(nx)
        Q = dt*hessian(l, x0 = x_[,k], u = u_[,k]) ## d^2l/dx^2
        R = dt*hessian(l, x0 = u_[,k], x = x_[,k]) ## d^2l/du^2
        q = dt*jacobian(l, x0 = x_[,k], u = u_[,k])
        r = dt*jacobian(l, x0 = u_[,k], x = x_[,k])
        s = dt*l(x_[,k],u_[,k])
        
        sys$A[1:nx,1:nx,k] = A; sys$A[nx+1,nx+1,k] = 1
        sys$B[1:nx,1:nu,k] = B; sys$B[nx+1,nu+1,k] = 0
        sys$Q[1:nx,1:nx,k] = Q; sys$Q[nx+1,1:nx,k] = q/2; sys$Q[1:nx,nx+1,k] = q/2; sys$Q[nx+1,nx+1,k] = s
        sys$R[1:nu,1:nu,k] = R; sys$R[nu+1,1:nu,k] = r/2; sys$R[1:nu,nu+1,k] = r/2; sys$R[nu+1,nu+1,k] = 0
        sys$H[,1:nx,k] = H
    }
    Q = hessian(h, x0 = x_[,N]); q = jacobian(h, x0 = x_[,N])
    sys$Q[1:nx,1:nx,N] = Q; sys$Q[nx+1,1:nx,N] = q/2; sys$Q[1:nx,nx+1,N] = q/2; sys$Q[nx+1,nx+1,N] = h(x_[,N])
    sys$Q[,,N] = 2*sys$Q[,,N]
    # noises
    sys$C0 = 0
    sys$C = 0
    sys$D0 = 0
    sys$D = 0
    sys$E0 = 0
    sys
}
