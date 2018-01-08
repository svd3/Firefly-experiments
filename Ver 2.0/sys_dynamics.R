######### system initialization ########
nv = nu = 2
ny = nw = nx = 2
N = 100
dt = 0.1
x0 = c(2,2)
xf = c(0,0)
Qf = 1000*diag(2) 
Q = array(0,c(2,2))
R = array(c(10,0,0,2), c(2,2))


####### System Dynamics ######
f = function(x, u) {
    ## f defines the deterministic part of the system dynamics. ##
    
    return(c(-x[2]*u[2], x[1]*u[2]-u[1]))
}

Fn = function(x,u){
    ## noise in state dynamics
    ## dim = nx X nw , we use nw = nx
    C0 = 0.005*(diag(nw)) #+ 0.2*diag(u)# 0.1*u%*%t(u) ## ## (Bu %*% t(Bu))
    C0
}

g = function(x, u) {
    ## g defines the deterministic part of the system observables. ##
    #H = diag(nx)
    #return(H%*%x)
    #return(c(-x[2]*u[2], x[1]*u[2]-u[1]))
    return(c(-u[2]*(1 + (x[1]/x[2])), (x[1]*u[2]-u[1])/x[2]^2))
}

G = function(x,u){
    ## noise in observation
    # dim = ny x nv , we use nv = ny 
    D0 = 0.005*(diag(ny)) #+ 0.01*diag(u)#0.5*u%*%t(u) ## (Hx %*% t(Hx))
    D0
}

####### Cost functions ########

l = function(x, u) {
    ## l defines the instantaneous costs prior to the final state. ##
    (t(x)%*%Q%*%x) + sum(u^2) + 100*(sum(u^2))^8 #(t(u)%*%R%*%u)+ (t(u)%*%R%*%u)
}

h = function(x) {
    ## h defines the system costs in the final state. ##
    (t(x)%*%Qf%*%x)
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
    }
    if(mode == "control"){
        ## initialize a trajectory
        if(rand){
            temp = rnorm(nu*(N-1),sd=1)
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






### Defining plant (actual world simulation) !!!
## gets control input and gives observation output and maintains state
init_plant = function(init_state){
    sys_x <<- init_state
    lockBinding("sys_x", globalenv())
}
plant = function(u){
    ## initialize state before running
    ## maintains states on its own
    
    unlockBinding("sys_x", globalenv())
    nw = dim(Fn(sys_x,u))[2]
    nv = dim(G(sys_x,u))[2]
    
    ## updates
    sys_x <<- sys_x + f(sys_x, u)*dt + Fn(sys_x, u)%*%rnorm(nw)
    y = g(sys_x, u) + G(sys_x, u)%*%rnorm(nv)
    lockBinding("sys_x", globalenv())
    y
}
stop_plant = function(){
    unlockBinding("sys_x", globalenv())
}

