setwd("Documents/Acads/FireflyProject-master/Rcodes/Ver 1.2")

source("ilqg.R")

# change if needed
Qf = 1000*diag(2) 
Q = array(0,c(2,2))
R = array(c(10,0,0,2), c(2,2))

## iterative ilqg
x0 = c(2,2)
dt = 0.1
N = 100

#res2 = gen_trajectory(x0,N,maxIter=10)
res = ilqg(x0, N, maxIter = 100, logg = T, mode = "lin_x", rand = T)
plot(t(res$x), col=2, type='l', xlim = c(-1,3), ylim = c(-1,3), asp=1, xlab = NA, ylab = NA); grid()
par(new=T)

W = to_World(res$u)
plot(t(W$x), col=3, type='l', xlim = c(-1,3), ylim = c(-1,3), asp=1, xlab = NA, ylab = NA); grid()
par(new=T)
vec_plot(W,xlim = c(-1,3), ylim = c(-1,3))


## Cost per time

cost_ps = 0
k = 1
kx = 10 ## important parameter relates to juice drop rewards
dt = 1
range_N = seq(5,30,by=1)
for(N in range_N){
    res = ilqg(x0, N, maxIter = 100, logg = T, mode = "lin_x")
    cost_ps[k] = (res$cost - kx)/N
    k = k+1
}

plot(range_N, cost_ps, type='o')

ku = exp(mean(log(cost_ps)+ 2*log(range_N))) #appx. and kx = 0
