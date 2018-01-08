######### system initialization ########
#N <<- as.integer(readline(prompt = "Enter N:  "))
nu = 2
nw = nx = 5
nv = ny = 2
dt = 0.1
N = 50
Tf = N*dt
x0 = c(0,0, pi/2,0,0)
xt = c(2,2,0,0,0)
Qf = 1000*diag(nx) # check Qf !!
Qf[3,3] = 0; #Qf[6,6] = 0
Q = array(0,c(nx,nx))
R = diag(nu) #array(c(1,0,0,1), c(2,2))

####### System Dynamics ######
f = function(x, u) {
    ## f defines the deterministic part of the system dynamics. ##
    v = x[4:5]; #psi = x[6]; tau = 0.05 + exp(psi); #10*sigmoid(psi,10,1.4)
    tau = 0.1
    amax = 1+log(1+exp(2*(tau-Tf)))/2
    # tau = 0 -- velocity
    # tau = infi -- acceleration
    #bias = c(vb, vb, wb, vb, wb)
    bias = 1
    gain = c(2,0.9)
    return(bias*c(v[1]*cos(x[3]), v[1]*sin(x[3]), v[2], (-v + amax*u*gain)/tau))#, -t_psi*psi-0.0))
}

f_act = function(x, u) {
    ## f defines the deterministic part of the system dynamics. ##
    v = x[4:5]; #psi = x[6]; tau = 0.05 + exp(psi); #10*sigmoid(psi,10,1.4)
    tau = 0.1
    amax = 1+log(1+exp(2*(tau-Tf)))/2
    # tau = 0 -- velocity
    # tau = infi -- acceleration
    gain = c(2,0.9)
    return(c(v[1]*cos(x[3]), v[1]*sin(x[3]), v[2], (-v + amax*u*gain)/tau))#, -t_psi*psi-0.0))
}

Fn = function(x,u){
    B = matrix(1, ncol=nu, nrow=nx)
    ## noise in state dynamics
    ## dim = nx X nw , we use nw = nx
    C0 = 0.01*(diag(nw)) #+ 0.00*diag(B%*%u)
    #C0[6,6] = 0
    C0
}

g = function(x, u) {
    ## g defines the deterministic part of the system observables. ##
    #bias = c(vb, wb)
    bias = 1
    bias*x[4:5]
}

G = function(x,u){
    H = matrix(0, ncol=nx, nrow=ny)
    H[,4:5] = diag(2)
    ## noise in observation
    # dim = ny x nv , we use nv = ny 
    D0 = 0.01*(diag(ny)) #+ 0.00*diag(H%*%x)
    D0
}

####### Cost functions ########

l = function(x, u) {
    ## l defines the instantaneous costs prior to the final state. ##
    u2 = (t(u)%*%R%*%u); u_c = (u[1]/1)^2 + (u[2]/1)^2
    (t(x-xt)%*%Q%*%(x-xt)) + u2 + 100*(u_c)^8
}

h = function(x) {
    ## h defines the system costs in the final state. ##
    (t(x-xt)%*%Qf%*%(x-xt))
    #xf = c(0.5,1); xf2 = c(1.5,0)
    #-exp(- t(x-xf)%*%(x-xf)) #- exp(-t(x-xf2)%*%(x-xf2))
}

compute_cost = function(x_,u_){
    N = dim(x_)[2]
    cost = h(x_[,N])
    for(k in 1:(N-1)){
        cost = cost + l(x_[,k],u_[,k])*dt
    }
    cost
}

control_cost = function(u_){
    N = dim(u_)[2] + 1
	cost = 0
	for(k in 1:(N-1)){
        cost = cost + (t(u_[,k]%*%R%*%u_[,k]))*dt
    }
    cost
}

state_cost = function(x_){
    N = dim(x_)[2]
	cost = h(x_[,N])
	for(k in 1:(N-1)){
        cost = cost + (t(x_[,k]%*%Q%*%x_[,k]))*dt
    }
    cost
}
###########################

### initial trajectory ######
init_trajectory = function(x0, N,rand =F){

    ## initialize a trajectory
    if(rand){
        temp = rnorm(nu*(N-1),sd=0.1)
    } else{
        temp = 0
    }
    u_ = array(temp,c(nu,N-1))
    x_ = array(0,c(nx,N))
    x_[,1] = x0
    for(k in 1:(N-1)){
        x_[,k+1] = x_[,k] + f(x_[,k], u_[,k])*dt
    }
    y_ = array(0,c(ny,N))
    for(k in 1:(N-1)){
        y_[,k+1] = g(x_[,k+1], u_[,k])
    }
    list(x_ = x_, u_ = u_, y_ = y_)
}

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
        x[,k+1] = x[,k] + f(x[,k],u[,k])*dt + Fn(x[,k],u[,k])%*%rnorm(nw)
        y[,k+1] = g(x[,k+1],u[,k]) + G(x[,k+1],u[,k])%*%rnorm(nv)
        xh[,k+1] = A[,,k]%*%xh[,k] + B[,,k]%*%u[,k] + cx[,k]
        yh[,k+1] = H[,,k]%*%xh[,k+1] + cy[,k]
        xh[,k+1] = xh[,k+1] + K[,,k]%*%(y[,k+1] - yh[,k+1])
        #xh[,k+1] = xh[,k] + f(xh[,k],u[,k])*dt + K[,,k]%*% (y[,k+1] - g(xh[,k], u[,k]))
    }
    return(list(x = x, y= y, u = u))
}