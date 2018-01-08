## system initialization
nu = 2
nx = 3
ny = 3
N = 100
dt = 0.1
x0 = c(0,0,pi/2)
xf = c(2,2,0); xf2 = c(2,1,0)
Qf = 1000*diag(3)
Qf[3,3] = 0
Q = array(0,c(3,3))
R = array(c(10,0,0,2), c(2,2))


f = function(x, u) {
    ## f defines the deterministic part of the system dynamics. ##
    return(c(u[1]*cos(x[3]), u[1]*sin(x[3]), u[2]))
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
    min(t(x-xf)%*%Qf%*%(x-xf), t(x-xf2)%*%Qf%*%(x-xf2))
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
