simulation = function(x0, res){
    N = dim(res$u)[2] + 1
    x = xh = array(0,c(nx,N)) 
    y = yh = array(0,c(ny,N))
    u = array(0,c(nu,N))
    x[,1] = xh[,1] = x0
    l = res$l
    L = res$L
    K = res$K
    A = res$A; B = res$B; H = res$H; cx = res$cx; cy = res$cy
    for(k in 1:(N-1)){
        u[,k] = l[,k] + L[,,k]%*%xh[,k]
        x[,k+1] = x[,k] + f(x[,k],u[,k])*dt #+ Fn(x[,k],u[,k])%*%rnorm(nw)
        y[,k+1] = g(x[,k+1],u[,k]) #+ G(x[,k+1],u[,k])%*%rnorm(nv)
        xh[,k+1] = A[,,k]%*%xh[,k] + B[,,k]%*%u[,k] + cx[,k]
        yh[,k+1] = H[,,k]%*%xh[,k+1] + cy[,k]
        xh[,k+1] = xh[,k+1] + K[,,k]%*%(y[,k+1] - yh[,k+1])
        #xh[,k+1] = xh[,k] + f(xh[,k],u[,k])*dt + K[,,k]%*% (y[,k+1] - g(xh[,k], u[,k]))
    }
    return(list(x = x, y= y, u = u))
}

f_act = function(x, u) {
    ## f defines the deterministic part of the system dynamics. ##
    v = x[4:5]; #psi = x[6]; tau = 0.05 + exp(psi); #10*sigmoid(psi,10,1.4)
    tau = 0.1
    amax = 1+log(1+exp(2*(tau-Tf)))/2
    # tau = 0 -- velocity
    # tau = infi -- acceleration
    return(c(v[1]*cos(x[3]), v[1]*sin(x[3]), v[2], (-v + amax*u*gain)/tau))#, -t_psi*psi-0.0))
}