## Iterative lqg - Updated 11-04-2017 (Working) 
## Version 2.1

source("Utils/newton.R")
source("Utils/to_World.R")
source("Utils/optimal_tf.R")
source("Kalman_LQG/kalman_lqg_2.R")
source("Linearize/linearize.R")
source("sys_dynamics.R")
source("simulate.R")

## add state dependent noise!!

ilqg = function(x0, N, maxIter=100, eps= 1e-6, logg = F,...){
    
    cost_vec = 0
    ## initialize a trajectory
    soln = init_trajectory(x0, N,...)
    x_ = soln$x_ ; u_ = soln$u_; y_ = soln$y_
    
    init_x = x_
    init_u = u_
    init_cost = compute_cost(x_,u_)  ## initial cost
    cost_vec[1] = init_cost
    
    if(logg) { 
        print(paste("initial cost = ", init_cost))
    }
    
    iter = 1
    converged = F
    
    while(!converged){
        update_nominal = T
        prev_cost = cost_vec[iter]
        
        ## Linearize
        sys = linearize(x_,u_,f,g,h,l,Fn,G)  ## linearize around initial nominal trajectory
        
        ## Compute gains
        gains = kalman_lqg(sys) ## lqg_all(sys)
        K = gains$K; L = gains$L; #lx = gains$l
        
        ## always check if Lu is 0 i.e. L[,(nx+1):(nx+nu),k] = 0 for all k
        
        ## line search to find control deviations (du)
        line_search = T
        alp = 1
        while(line_search){
            dx = dx_hat = array(0,c(nx,N))
            du = array(0,c(nu,N-1))
            dy = array(0, c(ny,N))
            u = array(0,c(nu,N-1))
            x = array(0,c(nx,N))
            y = array(0,c(ny,N))
            x[,1] = x0
            
            for(k in 1:(N-1)){
                Lx = -L[,1:nx,k] ; lx = -L[,nx+nu+1,k]
                Kx = K[1:nx,,k]
                
                du[,k] = Lx%*%dx_hat[,k] + alp*lx
                u[,k] = u_[,k] + du[,k]
                
                x[,k+1] = x[,k] + f(x[,k],u[,k])*dt #+ Fn(x[,k],u[,k])%*%rnorm(nw)
                y[,k+1] = g(x[,k+1], u[,k]) #+ G(x[,k+1],u[,k])%*%rnorm(nv)
                dy[,k+1] = y[,k+1] - y_[,k+1]
                dx[,k+1] = x[,k+1] - x_[,k+1]
                
                A = sys$A[1:nx,1:nx,k]; B = sys$B[1:nx,1:nu,k]; H = sys$H[,1:nx,k]
                
                dx_hat[,k+1] = A %*% dx_hat[,k] + B %*% du[,k] + 
                    Kx %*% (dy[,k+1] - H %*% dx_hat[,k]) #new
                
            }
            cost = compute_cost(x,u)
            if(cost>prev_cost){
                #line_search = T
                alp = 0.9*alp
                if(alp<1e-6){
                    alp = 0
                    line_search=F
                    update_nominal = F
                }
            } else{  ## cost <= prev_cost
                line_search = F
                update_nominal  =T
            }
        }
        
        ## set new nominal trajectory to the updated one before entering new iteration or exiting iterations
        if(update_nominal){
            x_ = x
            u_ = u
            y_ = y
        } 
        cost = compute_cost(x_,u_)
        
        if(logg){
            print(paste("iter", iter, ": cost =",cost))
        }
        ## convergence criteria: u_ doesn't change much or max iters
        ## norm(u-u_) < eps
        #eps = 1e-6
        if(iter >= maxIter | abs(cost-prev_cost)<eps){
            converged = T
        }
        iter = iter + 1
        cost_vec[iter] = cost
    }
    ## return
    lu = array(0, c(nu,N-1))
    Lu = array(0, c(nu,nx,N-1))
    Ku = array(0, c(nx,ny,N-1))
    cx = array(0, c(nx,N-1))
    cy = array(0, c(ny,N-1))
    A = array(0, c(nx,nx,N-1))
    B = array(0, c(nx,nu,N-1))
    H = array(0, c(ny,nx,N-1))
    for(k in 1:(N-1)){
        Ku[,,k] = K[1:nx,,k]
        Lu[,,k] = -L[,1:nx,k] ; lx = -L[,nx+nu+1,k]
        lu[,k] = u_[,k] - Lu[,,k]%*%x_[,k] + lx 
        A[,,k] = sys$A[1:nx,1:nx,k]; B[,,k] = sys$B[1:nx,1:nu,k]; H[,,k] = sys$H[,1:nx,k]
        cx[,k] = x_[,k+1] - A[,,k]%*%x_[,k] - B[,,k]%*%u_[,k]
        cy[,k] = y_[,k] - H[,,k]%*%x_[,k]
    }
    list(K = Ku, L = Lu, l = lu, cx=cx, cy = cy, A = A, B=B, H=H, x = x_, y = y_, u = u_, cost = cost, dt = dt)
}