source("iLQG.R")
## state defination : x = (x, y, o (angle in rad))	position and orientation 
## control variables : u = (v, w) velocity and angular velocity

# define the system for iLQG

f = function(x, u) {
	## f defines the deterministic part of the system dynamics. ##
	#print("f :: ")
	return(c(-x[2]*u[2], x[1]*u[2]-u[1]))
}

Fn = function(x, u) {
	## F defines the stochastic part of the system dynamics. ##
	#print("F :: ")
	return(C0)
}

g = function(x, u) {
	## g defines the deterministic part of the system observables. ##
	#print("g :: ")
	return(u)
}

G = function(x, u) {
	## G defines the stochastic part of the system observables. ##
	#print("G :: ")
	return(D0)
}

l = function(x, u) {
	## l defines the system costs prior to the final state. ##
	#print("l :: cost ::")
	(t(x)%*%Q%*%x) + (t(u)%*%R%*%u)
}

h = function(x, u) {
	## h defines the system costs in the final state. ##
	#print("h :: final cost ::")
	return(0)
}

#dt = 0.1

gen_ilqg = function(f, Fn, g, G, h, l, x0, d, xf = 0){
    #print("gen_ilqg :: ")
    N = d$N
    ## nu = dim of control variables
    #dt = 1 ## time update
    nx = d$nx # length(x0) ## dim of state variables
    nu = d$nu; 
    u_n = array(0, c(nu, N-1)) ## nominal control
    u_p = u_n ## actual control
    x_n = array(0, c(nx, N)) ## nominal control
    x_p = x_n ## actual control
    x_hat = x_n ## estimate of x (deviation) from observable
    x_p[,1] = x_n[,1] = x0
    ny = d$ny #length(g(x_n[,1], u_n[,1])) ## dim of observables
    y_n = array(0, c(ny, N)) ## nominal obs
    y_p = y_n ## actual obs
    Lx = array(0, dim = c(nu, nx, N-1)) ## for control law (lx + Lx . X)
    lx = array(0, dim = c(nu, N-1))
    K = array(0, dim = c(nx, ny, N-1)) ## Kalman Gain
	
    nw = d$nw #0 ## dim(Fn(x_n[,0], u_n[,0]))[1]
    nv = d$nv #0 ## dim(G(x_n[,0], u_n[,0]))[1]
    d0 = d
    for(k in 1:(N-1)) {
        print(k)
        if(k==1 | norm2(x_p[,k]-x_n[,k]) > 0.1*norm2(x_n[,k]) ){
            d$N = N-k+1
            solution = iterative_lqg(f,Fn,g,G,h,l, x_p[,k], d, xf)
            #print(size(x_n))
            x_n[,k:N] = solution[[1]]
            u_n[,k:(N-1)] = solution[[2]]
            Lx[,,k:(N-1)] = solution[[3]]
            lx[,k:(N-1)] = solution[[4]]
            K[,,k:(N-1)] = solution[[5]]
            systm = linearize_and_quadratize(f,Fn,g,G,h,l, x_n, u_n,d0)
            # calculate the nominal observations
            for(j in 1:(N-1)){
            		y_n[ ,j] = g(x_n[,j], u_n[,j]) ## change this part
            }
        }
        u_p[,k] = (Lx[,,k]%*%x_hat[,k]) + lx[,k] + u_n[,k]  # calculate the control input
        x_p[,k+1] = x_p[,k] + f(x_p[,k], u_p[,k])*dt
                            + Fn(x_p[,k], u_p[,k])%*%rnorm(nw)  # calculate the next state
        if(k==0){
            y_p[,k] = y_n[,k]
        }
        y_p[,(k+1)] = g(x_p[,k],u_p[,k]) + G(x_p[,k], u_p[,k])%*%rnorm(nv)
        
        y = y_p[,k] - y_n[,k]
        #print(dim(systm$A))
        A = systm$A[1:nx,1:nx,k]
        B = systm$B[1:nx,,k]
        H = systm$H[,1:nx,k]
        u = u_p[,k] - u_n[,k]
        x_hat[,k+1] = A%*%x_hat[,k] + B%*%u + K[,,k]%*%(y - H%*%x_hat[,k])
    }
    list(x_p, u_p)
}
