## updated on 06-03-2017
## Working

#setwd("Documents/Acads/FireflyProject-master/Rcodes/Ver 1.1")
source("ilqg.R")
source("to_World.R")

## iterative ilqg
#x0 = c(2,2)
#N = 20
#res = gen_trajectory(x0,N)

#ilqg_par = list(maxIter = 50, eps = 1e-3, logg = T)
gen_trajectory = function(x0, N,...){
    soln = init_trajectory(x0,N)
    init_x = x = x_ = soln$x_ ;  init_u = u = u_ = soln$u_
    prev_cost = init_cost = compute_cost(x_,u_)
    print(paste("initial cost = ", init_cost))
    iter = 1
    for(n in N:2){
        soln = ilqg(x0, n, ...) # maxIter = 50, eps=1e-3, logg =T)
        x_[,(N-n+1):N] = soln$x
        u_[,(N-n+1):(N-1)] = soln$u
        cost = compute_cost(x_,u_)
        
        if(cost < prev_cost){
          ## take steps
          ## update actual trajectory
          x = x_
          u = u_
        }
        cost = compute_cost(x,u)
        print(paste("time_step", iter, ": cost =",cost))
        prev_cost = cost
        x0 = x[,(N-n+2)]
        iter = iter+1
    }
    #x_ = x
    #u_ = u
    list(x = x, u = u, cost = cost, init_x = init_x, init_u = init_u, init_cost = init_cost, dt = dt)
}
