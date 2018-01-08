## Iterative lqg - Updated 03-04-2017 (Working)
source("newton.R")
source("kalman_lqg_2.R")
source("sys_dynamics.R")
source("linearize.R")
source("lqr.R")
source("to_World.R")

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
        sys = linearize(x_,u_,f,g,h,l,Fn,G)  ## linearize around initial nominal trajectory
        if(noise){
            gains = kalman_lqg(sys) ## find gains
            K = gains$K; L = gains$L
        } else{
            ## can use this for simple additive noise (much much faster)
            gains = lqr(sys)
            K = gains$K; L = gains$L
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
                #Lx = L[,1:nx,k];  l = L[,nx+nu+1,k]
                ## apply line search
                y_noise = 0
                for(i in 1:nv){
                    y_noise = y_noise + sys$D[,1:nx,i,k]%*%dx[,k]*rnorm(1)
                }
                x_noise = 0
                for(i in 1:nw){
                    x_noise = x_noise + sys$C[1:nx,,i,k]%*%du[,k]*rnorm(1)
                }
                
                dy[,k] = sys$H[,1:nx,k]%*%dx[,k] + sys$D0[,,k]%*%rnorm(ny) + y_noise #new

                du[,k] = -L[,1:nx,k]%*% dx_hat[,k] - alp*L[,nx+nu+1,k] #new

                dx_hat[,k+1] = sys$A[1:nx,1:nx,k]%*%dx_hat[,k] + sys$B[1:nx,1:nu,k]%*%du[,k] + 
                K[1:nx,,k]%*%(dy[,k] - sys$H[,1:nx,k]%*%dx_hat[,k]) #new

                dx[,k+1] = sys$A[1:nx,1:nx,k]%*%dx[,k] + sys$B[1:nx,1:nu,k]%*%du[,k] + 
                sys$C0[1:nx,1:nx,k]%*%rnorm(nx) + x_noise

                u[,k] = u_[,k] + du[1:nu,k]
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
