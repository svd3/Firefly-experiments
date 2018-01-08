
f = function(x, u) {
    ## f defines the deterministic part of the system dynamics. ##
    return(c(-x[2]*u[2], x[1]*u[2]-u[1]))
}

g = function(x, u) {
    ## g defines the deterministic part of the system observables. ##
    return(x)
}

l = function(x, u) {
    ## l defines the instantaneous costs prior to the final state. ##
    (t(x)%*%Q%*%x) + (t(u)%*%R%*%u)
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
            u_[,k] = root2(f,xdot,x0=c(1,1), x=x_[,k])
        }
    }
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


#plant = function(x,u){
 #   ## gives next output
  #  x_new = x + f(x,u)*dt
   # y_new = g(x,u)
    #return list(x_new, y_new)
#}

## system initialization
nu = 2
nx = 2
ny = 2
N = 100
dt = 0.1
x0 = c(2,2)
xf = c(0,0)
Qf = 1000*diag(2) 
Q = array(0,c(2,2))
R = array(c(10,0,0,2), c(2,2))
