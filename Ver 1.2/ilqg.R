## Iterative lqg - Updated 06-03-2017 (Working)
source("newton.R")
source("kalman_lqg.R")
source("sys_dynamics.R")
source("linearize.R")
source("lqr.R")
source("to_World.R")


ilqg = function(x0, N, maxIter=100, eps  = 1e-6, noise = F, logg = F, plt = F,...){
	
    cost_vec = 0
	## initialize a trajectory
    soln = init_trajectory(x0, N,...) # , mode= "lin_x")
    x_ = soln$x_ ; u_ = soln$u_
	
	## for checking if controls above create the trajectory
    "if(F){
        x_[,1] = x0
        for(k in 1:(N-1)){
            xdot = (x_[,k+1]-x_[,k])/dt
            u_[,k] = root2(f,xdot,x0=c(1,1), x=x_[,k])
            x_[,(k+1)] = x_[,k] + f(x_[,k],u_[,k])*dt
        }
    }"
	
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
        #prev_cost = compute_cost(x_,u_)
	    prev_cost = cost_vec[iter]
        sys = linearize(x_,u_,f,g,h,l)  ## linearize around initial nominal trajectory
	    if(noise){
	    	gains = kalman_lqg(sys) ## find gains
	    	K = gains$K; L = gains$L
	    } else{
	    	L = lqr(sys)
	    }
        ## always check if Lu is 0 i.e. L[,(nx+1):(nx+nu),k] = 0 for all k
        
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
	            #Lx = L[,1:nx,k];  l = L[,nx+nu+1,k]
	            ## apply line search
	            du[,k] = -L[,1:nx,k]%*% dx[,k] - alp*L[,nx+nu+1,k]
	            dx[,k+1] = sys$A[1:nx,1:nx,k]%*%dx[,k] + sys$B[1:nx,1:nu,k]%*%du[,k]
	            u[,k] = u_[,k] + du[1:nu,k]
	            x[,k+1] = x[,k] + f(x[,k],u[,k])*dt
	        }
	        cost = compute_cost(x,u)
	        if(cost>prev_cost){
	            #line_search = T
	            alp = 0.9*alp
	            if(alp<1e-13){
	                alp = 0
	            }
	        } else{
	            line_search = F
	        }
	        #print(abs(cost-prev_cost))
	    }
	    if(logg){
		    print(paste("iter", iter, ": cost =",cost))
		}
	    
	    ## convergence criteria: u_ doesn't change much or max iters
	    ## norm(u-u_) < eps
	    #eps = 1e-6
	    if(iter >= maxIter | abs(cost-prev_cost)<eps){
	        converged = T
	    }
	    
	    ## set new nominal trajectory to the updated one before entering new iteration or exiting iterations
	    x_ = x
	    u_ = u
	    iter = iter + 1
        cost_vec[iter] = cost
	}
    
    ## plotting
    if(plt){
        rx = range(x_[1,]); ry = range(x_[1,])
        plot(t(x_), type='l', col=1, xlim = rx, ylim = ry, asp=1, xlab = NA, ylab = NA)
        par(new=T)
        plot(t(init_x), type='l', col=2, xlim = rx, ylim = ry,asp=1, xlab = NA, ylab = NA)
        title(xlab = "x", ylab = "y")
        grid()
    }
    
    ##return
	list(x = x_, u = u_, cost = cost, cost_vec = cost_vec, init_x = init_x, init_u = init_u, init_cost = init_cost, dt = dt)
}

init_trajectory = function(x0, N, mode = "control", rand =F){
    if(mode == "lin_x"){
        ## initialize a trajectory
        u_ = array(0,c(nu,N-1))
        x_ = array(0,c(nx,N))
        for(i in 1:nx){
            x_[i,] = seq(x0[i],xf[i], length = N)
        }
        for(k in 1:(N-1)){
            xdot = (x_[,k+1]-x_[,k])/dt
            u_[,k] = root2(f,xdot,x0=rep(1,nx), x=x_[,k])
        }
    }
    if(mode == "control"){
        ## initialize a trajectory
        if(rand){
            temp = rnorm(nu*(N-1),sd=0.1)
        } else{
            temp = 0.1
        }
        u_ = array(temp,c(nu,N-1))
        x_ = array(0,c(nx,N))
        x_[,1] = x0
        for(k in 1:(N-1)){
            x_[,k+1] = x_[,k] + f(x_[,k], u_[,k])*dt
        }
    }
    list(x_ = x_, u_ = u_)
}
