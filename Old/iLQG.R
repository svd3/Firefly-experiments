#library(numDeriv) ## jacobian and hessian
library(pracma) ## jacobian
library(MASS)  ## ginv
source("Kalman_LQG.R")

#dt = 1.0
#size = function(x){ if(is.null(dim(x))){ length(x)} else{dim(x)}}
norm2 = function(x) sqrt(sum(x^2))

compute_control_trajectory = function(f, x, d){
	#print("compute_control_trajectory :: ")
	## Compute the control trajectory from the state trajectory and
    ## the function that describes the system dynamics.
    nu = d$nu
    N = d$N #size(x)[2]
    u = array(0, c(nu, N-1))
    ## Calculate the the control input estimate, u_hat, for the first time step.
    u_hat = rep(0, nu)
    #dt = 1.0  # until there's a reason to use something else
    dx = x[,2] - x[,1]
    dfdu = jacobian(f, x0 = u_hat, x = x[,1])
    u_hat = ginv(dfdu)%*%(dx/dt - f(x[,1], u_hat))
    for(k in 1:(N-1)){
    	dfdu = jacobian(f, x0 = u_hat, x = x[,k])
    	dx = x[,k+1] - x[,k]
    	du = ginv(dfdu)%*%(dx/dt - f(x[,k], u_hat))
    	u_hat = u_hat + du
    	u[,k] = u[,k] + u_hat
    }
    u
}


compute_state_trajectory = function(f, x0, u,d){
	#print("compute_state_trajectory :: ")
    N = d$N #size(u)[2] + 1
    nx = d$nx #size(x0)[1]
    x = array(0, c(nx, N))
    x[,0] = x0
    for(k in 1:(N-1)){
    	x[,k+1] = x[,k] + f(x[,k], u[,k])*dt
    }
    x
}

initial_trajectory = function(f, x0, xf, d){
	#print("initial_trajectory :: ")
	#nx = size(x0)[1]
    nx = d$nx; nu = d$nu; N = d$N
	x = array(0, c(nx, N))
	for(i in 1:nx){ x[i,] = seq(x0[i], xf[i], len = N) }
	#x[2,] = x0[2]*exp(-1*0:(N-1))
	#dx = (xf - x0)/(N-1)
	#for(i in 1:N){ x[,i] = x0 + (i-1)*dx }
	u = compute_control_trajectory(f, x, d)
	list(x=x, u=u)
}


## lets create a system dim variable contatining all dimensions !!
linearize_and_quadratize = function(f, Fn, g, G, h, l_cost, x, u, d){
	#print("linearize_and_quadratize :: ")
	nx = d$nx #size(x)[1]
	nu = d$nu #size(u)[1]
	nxa = nx + nu + 1
	szC0 = d$nw #size(Fn(x[,0], u[,0]))[2]
	ny = d$ny #size(g(x[,0], u[,0]))[1]
	szD0 = d$nv #size(G(x[,0], u[,0]))[2]
	N = d$N #size(x)[2]
	
	systm = list()
	
	# build the vector for the initial augmented state
	x0a = rep(0, nx)
	u0a = c(rep(0, nu), 1)
	
	systm$X1 = c(x0a, u0a)
	S1 = diag(nxa)
	S1[nxa, nxa] = 0
	systm$S1 = S1
	systm$A = array(0, c(nxa, nxa, N-1))
	systm$B = array(0, c(nxa, nu, N-1))
	systm$C0 = array(0, c(nxa, szC0, N-1))
	systm$C = array(0, c(nu, nu, szC0, N-1))
	systm$H = array(0, c(ny, nxa, N-1))
	systm$D0 = array(0, c(ny, szD0, N-1))
	systm$D = array(0, c(ny, nxa, szD0, N-1))
	systm$Q = array(0, c(nxa, nxa, N))
	systm$R = array(0, c(nu, nu, N-1))
	
	for(k in 1:(N-1)){
		dfdx = jacobian(f, x0 = x[,k], u = u[,k]); dim(dfdx) = c(nx,nx)
		A = dfdx*dt + diag(nx)
		dfdu = jacobian(f, x0 = u[,k], x = x[,k]); dim(dfdu) = c(nx,nu)
		B = dfdu*dt
		C0 = sqrt(dt)*Fn(x[,k], u[,k])
		dFdu = jacobian(Fn, x0 = u[,k], x = x[,k]); dim(dFdu) = c(nx, szC0, nu)
		dFdu = aperm(dFdu, c(1,3,2))
		C = array(0, c(nu, nu, szC0))
		for(j in 1:szC0){ C[,,j] = ginv(B) %*% dFdu[,,j] }
		C = sqrt(dt)*C
        H = dgdx = jacobian(g, x0 = x[,k], u = u[,k]); dim(dgdx) = c(nu,nx)
        dGdx = jacobian(G, x0 = x[,k], u = u[,k]); dim(dGdx) = c(ny,szD0,nx)
        dGdx = aperm(dGdx, c(1,3,2))
        D = dGdx/dt
        qs = dt*l(x[,k],u[,k])
        dldx = jacobian(l, x0 = x[,k], u = u[,k])
        q = dt*dldx
        d2ldx2 = hessian(l, x0 = x[,k], u = u[,k])
        Q = dt*d2ldx2
        if(k==1){
            r = array(0, nu)
            R = array(0, c(nu,nu))
        } else{
            dldu = jacobian(l, x0 = u[,k-1], x = x[,k-1])
            d2ldu2 = hessian(l, x0 = u[,k-1], x = x[,k-1])
            r = dt*dldu
            R = dt*d2ldu2
        }
        Aa = array(0, c(nxa,nxa))
        Aa[1:nx,1:nx] = A
        Aa[nxa,nxa] = 1
        systm$A[,,k] = Aa
        Ba = array(0,c(nxa,nu))
        Ba[1:nx,] = B
        Ba[(nx+1):(nx+nu),] = diag(nu)
        systm$B[,,k] = Ba
        C0a = array(0, c(nxa,szC0))
        C0a[1:nx,] = C0
        systm$C0[,,k] = C0a
        Ha = array(0, c(ny,nxa))
        Ha[,1:nx] = H
        systm$H[,,k] = Ha
        for(j in 1:size(D)[3]){
            Da = array(0, c(ny,nxa))
            Da[,1:nx] = D[,,j]
            systm$D[,,j,k] = Da
        }
        Qa = array(0, c(nxa,nxa))
        Qa[1:nx,1:nx] = Q
        Qa[1:nx,nx+nu] = Qa[nx+nu,1:nx] = q/2
        Qa[(nx+1):(nx+nu),(nx+1):(nx+nu)] = R
        Qa[(nx+1):(nx+nu),nx+nu] = Qa[nx+nu,(nx+1):(nx+nu)] = r/2
        Qa[nxa,nxa] = qs
        systm$Q[,,k] = Qa
        systm$R[,,k] = array(0, c(nu,nu))
	}
    qs = h(x[,N])
    q = dhdx = jacobian(h, x0 = x[,N])
    Q = d2hdx2 = hessian(h, x0 = x[,N])
    dldu = jacobian(l, x0 = u[,N-1], x = x[,N-1])
    r = dt*dldu
    d2ldu2 = hessian(l, x0 = u[,N-1], x = x[,N-1])
    R = dt*d2ldu2
    Qa = array(0, c(nxa,nxa))
    Qa[1:nx,1:nx] = Q
    Qa[1:nx,nx+nu] = Qa[nx+nu,1:nx] = q/2
    Qa[(nx+1):(nx+nu),(nx+1):(nx+nu)] = R
    Qa[(nx+1):(nx+nu),nx+nu] = Qa[nx+nu,(nx+1):(nx+nu)] = r/2
    Qa[nxa,nxa] = qs
	systm$Q[,,k] = Qa
    # iLQG does not accommodate noise added to the state estimate
    systm$E0 = array(0,c(1,1,N))
    systm
}


compute_cost = function(h, l, x, u, d){
	#print("compute_cost :: ")
	cost = 0
	N = d$N #size(x)[2]
	for(k in 1:(N-1)){
		cost = cost + l(x[,k], u[,k])*dt
	}
	cost = cost + h(x[,N])
	cost
}

update_trajectory = function(f, x_n, u_n, La){
	#print("update_trajectory :: ")
	N = size(La)[3] + 1
    nxa = size(La)[2]
	nu = size(u_n)[1]
	nx = size(x_n)[1]
	#print(nx)
	u_p = array(0, c(nu, N-1))
	x_p = array(0, c(nx, N))
    x_p[,1] = x_n[,1]
    l = array(0, c(nu, N-1))
    L = array(0, c(nu,nx,N-1))
    Lu = array(0, c(nu,nu,N-1))
    #print(size(L))
    #print(size(La))
    for(k in 1:(N-1)){
        x = x_p[,k] - x_n[,k]
        L[,,k] = La[,1:nx,k]
        l[,k] = La[,nxa,k]
        Lu[,,k] = La[,(nx+1):(nx+nu),k]
        if(length(which(Lu!=0))!=0){
            print("ERROR: file iLQG.R, fn update_trajectory:: Lu != 0")
        }
        u = -l[,k] - L[,,k]%*%x
        #print(u)
        u_p[,k] = u_n[,k] + u
        x_p[,k+1] = x_p[,k] + f(x_p[,k],u_p[,k])*dt
    }
    list(x_p=x_p, u_p=u_p, L=L, l=l)
}

iterative_lqg = function(f,Fn,g,G,h,l,x0,d,xf=0){
	#print("iterative_lqg :: ")
    N = d$N
    nx = d$nx #length(x0)
    if(xf==0){
        xf = array(0,nx)
    }
    nu = d$nu
    soln = initial_trajectory(f,x0,xf,d)
    x_n = soln$x; u_n = soln$u
    #print(u_n)
    cost = compute_cost(h,l,x_n,u_n,d)
    print(paste("init tarjectory cost:", cost))
    systm = linearize_and_quadratize(f, Fn, g, G, h, l, x_n, u_n, d)
    
    init_systm = systm
    
    soln = kalman_lqg(systm)
    #print(soln)
    K = soln$K; L = soln$L; Cost = soln$Cost;
    Xa = soln$Xa; XSim = soln$XSim; CostSim = soln$CostSim; iters = soln$iters
    print(L[,,1])
    hasnt_converged = T
    iter = 1
    while(hasnt_converged){
        if(FALSE){
            ## plotting if needed
        }
        trajectory = update_trajectory(f, x_n, u_n, L)
        x_n = trajectory[[1]]; u_n = trajectory[[2]];
        L_n = trajectory[[3]]; l_n = trajectory[[4]];
        #print(u_n)
        prev_cost = cost
        cost = compute_cost(h,l,x_n,u_n,d)
        print(paste0("ilqg iter:", iter, ":: cost:", cost))
        if(abs(cost-prev_cost) < 0.1){
            ## check condition
          hasnt_converged = FALSE
        } else{
            systm = linearize_and_quadratize(f, Fn, g, G, h, l, x_n, u_n,d)
            soln = kalman_lqg(systm)
            K = soln$K; L = soln$L; Cost = soln$Cost;
            Xa = soln$Xa; XSim = soln$XSim; CostSim = soln$CostSim; iters = soln$iters
        }
        #print(hasnt_converged)
        iter = iter +1
    }
    list(x_n, u_n, L_n, l_n, K[1:nx,,])
}
