## Iterative lqg - Updated 11-04-2017 (Working) 
## Version 2.1

source("Utils/newton.R")
source("Utils/to_World.R")
source("Kalman_LQG/lqr.R")
source("Kalman_LQG/new_kalman.R")
source("Linearize/linearize2.R")
source("sys_dynamics.R")

## add state dependent noise!!

ilqg = function(x0, N, maxIter=100, eps  = 1e-6, noise = F, logg = F, plt = F,...){

    cost_vec = 0
    ## initialize a trajectory
    soln = init_trajectory(x0, N,...) # , mode= "lin_x")
    x_ = soln$x_ ; u_ = soln$u_

    init_x = x_
    init_u = u_
    init_cost = compute_cost(x_,u_)  ## initial cost
    cost_vec[1] = init_cost

    if(logg) { 
        print(paste("initial cost = ", init_cost))
    }

    iter = 1
    #maxIter = 100
    converged = F

    while(!converged){
        update_nominal = T
        #prev_cost = compute_cost(x_,u_)
        prev_cost = cost_vec[iter]
        sys = linearize2(x_,u_,f,g,h,l,Fn,G)  ## linearize around initial nominal trajectory
        if(noise){
            #gains = kalman_lqg(sys)
            gains = lqg_all(sys) ## find gains
            K = gains$K; L = gains$L; lx = gains$l
        } else{
            ## can use this for simple additive noise (much much faster)
            gains = lqr(sys)
            K = gains$K; L = -gains$L; lx = 0 ## negative to match lqg_all convention
            ## now we have " u = Lx "
        }
        #print(dim(L))
        ## always check if Lu is 0 i.e. L[,(nx+1):(nx+nu),k] = 0 for all k

        ## line search to find control deviations (du)
        ## line search minimizes cost w.r.t previous iteration
        line_search = T
        alp = 1
        while(line_search){
            dx = dx_hat = array(0,c(nx,N))
            du = array(0,c(nu,N-1))
            dy = array(0, c(ny,N))
            u = array(0,c(nu,N-1))
            x = array(0,c(nx,N))
            x[,1] = x0
            for(k in 1:(N-1)){
                # dx_hat[1] = dx[1] (first element)
                Ck = array(0,c(nx,nw))
                Dk = array(0,c(ny,nv))
                for(i in 1:nv){
                    Dk[,i] = sys$D0[,i,k] + sys$Dx[,,i,k]%*%dx[,k] + sys$Du[,,i,k]%*%du[,k]
                }
                for(i in 1:nw){
                    Ck[,i] = sys$C0[,i,k] + sys$Cx[,,i,k]%*%dx[,k] + sys$Cu[,,i,k]%*%du[,k]
                }
                
                du[,k] = L[,,k]%*%dx_hat[,k] + alp*lx[,k]
                #du[,k] = -L[,1:nx,k]%*% dx_hat[,k] - alp*L[,nx+nu+1,k] #new
                dy[,k] = sys$H[,,k]%*%dx[,k] +sys$E[,,k]%*%du[,k] + Dk%*%rnorm(nv)

                dx_hat[,k+1] = sys$A[,,k]%*%dx_hat[,k] + sys$B[,,k]%*%du[,k] + 
                    K[,,k]%*%(dy[,k] - (sys$H[,,k]%*%dx_hat[,k] +sys$E[,,k]%*%du[,k])) #new

                dx[,k+1] = sys$A[,,k]%*%dx[,k] + sys$B[,,k]%*%du[,k] + Ck%*%rnorm(nw)

                u[,k] = u_[,k] + du[,k]
                x[,k+1] = x[,k] + f(x[,k],u[,k])*dt
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
            } else{
                line_search = F
                update_nominal  =T
            }
            #print(alp)
        }

        ## set new nominal trajectory to the updated one before entering new iteration or exiting iterations
        if(update_nominal){
            x_ = x
            u_ = u
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

    ## plotting
    if(plt){
        rx = range(x_[1,]); ry = range(x_[1,]); rx=ry = c(-1,3)
        plot(t(x_), type='l', col=1, xlim = rx, ylim = ry, asp=1, xlab = NA, ylab = NA)
        #par(new=T)
        #plot(t(init_x), type='l', col=2, xlim = rx, ylim = ry,asp=1, xlab = NA, ylab = NA)
        title(xlab = "x", ylab = "y")
        grid()
    }

    ##return
    list(x = x_, u = u_, cost = cost, cost_vec = cost_vec, init_x = init_x, init_u = init_u, init_cost = init_cost, dt = dt)
}
