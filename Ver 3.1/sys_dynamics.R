######### system initialization ########
nv = nu = 2
nw = nx = 5
ny = 2
N = 100
dt = 0.1
x0 = c(0,0, pi/2,0,0)
xt = c(2,2,0,0,0)
Qf = 1000*diag(5) # check Qf !!
Qf[3,3] = 0
Q = array(0,c(5,5))
R = diag(2) #array(c(1,0,0,1), c(2,2))


####### System Dynamics ######
f = function(x, u) {
    ## f defines the deterministic part of the system dynamics. ##
    v = x[4:5]
    # tau = 0 -- acceleration
    # tau = infi -- velocity
    return(c(v[1]*cos(x[3]), v[1]*sin(x[3]), v[2], (-tau*v + amax*u)))
}

Fn = function(x,u){
    ## noise in state dynamics
    ## dim = nx X nw , we use nw = nx
    C0 = 0.01*(diag(nw)) #+ 0.01*diag(u) ## (Bu %*% t(Bu))
    C0
}

g = function(x, u) {
    ## g defines the deterministic part of the system observables. ##
    u
}

G = function(x,u){
    ## noise in observation
    # dim = ny x nv , we use nv = ny 
    D0 = 0.01*(diag(ny)) + 0.01*diag(u)
    D0
}

####### Cost functions ########

l = function(x, u) {
    ## l defines the instantaneous costs prior to the final state. ##
    u2 = (t(u)%*%R%*%u); u_c = (u[1]/1)^2 + (u[2]/1)^2
    (t(x)%*%Q%*%x) + u2 + 100*(u_c)^8
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
init_trajectory = function(x0, N, mode = "control", rand =F){
    if(F){
        if(mode == "lin_x"){
        ## initialize a trajectory
        u_ = array(0,c(nu,N-1))
        x_ = array(0,c(nx,N))
        for(i in 1:nx){
            x_[i,] = seq(x0[i],xf[i], length = N)
        }
        for(k in 1:(N-1)){
            xdot = (x_[,k+1]-x_[,k])/dt
            u_[,k] = root2(f,xdot,x0=rep(1,nx),eps = 1e-2, x=x_[,k])
        }
    }}
    if(mode == "control"){
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
    }
    list(x_ = x_, u_ = u_)
}