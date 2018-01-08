## Iterative lqg - Updated 03-03-2017 (Working)
setwd("Documents/Acads/FireflyProject-master/Rcodes/Ver 1.0")
source("newton.R")
source("kalman_lqg.R")
source("sys_dynamics.R")
source("linearize.R")
source("lqr.R")

run = function(x0, N, maxIter=100, plt = T){
	
	all_x = list()
    cost_vec = 0
    
	## initialize a trajectory
	u_ = array(0,c(nu,N-1))
	x_ = array(0,c(nx,N))
	for(i in 1:nx){
	    x_[i,] = seq(x0[i],xf[i], length = N)
	}
	for(k in 1:(N-1)){
	    xdot = (x_[,k+1]-x_[,k])/dt
	    u_[,k] = root2(f,xdot,x0=c(1,1), x=x_[,k])
	}
	
	## for checking if controls above create the trajectory
    "if(F){
        x_[,1] = x0
        for(k in 1:(N-1)){
            xdot = (x_[,k+1]-x_[,k])/dt
            u_[,k] = root2(f,xdot,x0=c(1,1), x=x_[,k])
            x_[,(k+1)] = x_[,k] + f(x_[,k],u_[,k])*dt
        }
    }"
	
	init_x = t(x_)
    init_u = t(u_)
    all_x[[1]] = t(x_)
	init_cost = compute_cost(x_,u_)  ## initial cost
	print(paste("initial cost = ", init_cost))
    cost_vec[1] = init_cost
    #sys = linearize(x_,u_,f,g,h,l)  # linearize around initial nominal trajectory
    #init_sys = sys ## why is this needed?
	
    iter = 1
	#maxIter = 100
	converged = F
    
	while(!converged){
        #prev_cost = compute_cost(x_,u_)
	    prev_cost = cost_vec[iter]
        sys = linearize(x_,u_,f,g,h,l)  ## linearize around initial nominal trajectory
	    #gains = kalman_lqg(sys) ## find gains
	    #K = gains$K; L = gains$L
	    L = lqr(sys)
	
	    ## line search to find control deviations (du)
        ## line search minimizes cost w.r.t previous iteration
	    line_search = T
	    alp = 1
	    while(line_search){
	        dx = array(0,c(nx,N))
	        du = array(0,c(nu,N-1))
	        u = array(0,c(nu,N-1))
	        x = array(0,c(nx,N))
	        x[,1] = x0
	        for(k in 1:(N-1)){
	            ## apply line search
	            du[,k] = -L[1:nu,1:nx,k]%*% dx[,k] - alp*L[1:nu,nx+1,k]
	            dx[,k+1] = sys$A[1:nx,1:nx,k]%*%dx[,k] + sys$B[1:nx,1:nu,k]%*%du[,k]
	            u[,k] = u_[,k] + du[1:nu,k]
	            x[,k+1] = x[,k] + f(x[,k],u[,k])*dt
	        }
	        cost = compute_cost(x,u)
	        if(cost>prev_cost){
	            #line_search = T
	            alp = 0.9*alp
	        } else{
	            line_search = F
	        }
	    }
	    
        #cost = compute_cost(x,u)
	    print(paste("iter", iter, ": cost =",cost))
	    
	    ## convergence criteria: u_ doesn't change much or max iters
	    ## norm(u-u_) < eps
	    if(iter >= maxIter | abs(cost-prev_cost)<1e-6){
	        converged = T
	    }
	    
	    ## set new nominal trajectory to the updated one before entering new iteration or exiting iterations
	    x_ = x
	    u_ = u
	    iter = iter + 1
        all_x[[iter]] = t(x_)
        cost_vec[iter] = cost
	}
    x_ = t(x_)
    u_ = t(u_)
    
    ## plotting
    if(plt){
        rx = range(x_[,1]); ry = range(x_[,2])
        plot(x_, type='l', col=1, xlim = rx, ylim = ry, asp=1, xlab = NA, ylab = NA)
        par(new=T)
        plot(init_x, type='l', col=2, xlim = rx, ylim = ry,asp=1, xlab = NA, ylab = NA)
        title(xlab = "x", ylab = "y")
        grid()
    }
    
    ##return
	list(all_x = all_x, x = x_, u = u_, cost = cost_vec)
}

