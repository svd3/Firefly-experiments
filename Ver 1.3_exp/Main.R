setwd("Documents/Acads/FireflyProject-master/Rcodes/Ver 1.3_exp")

source("ilqg.R")

# change if needed
x0 = c(0,0,pi/2)
xf = c(2,2,0)
Qf = 1*diag(3)
Qf[3,3] = 0
Q = array(0,c(3,3))
R = array(c(10,0,0,2), c(2,2))

## iterative ilqg
#x0 = c(2,2)
dt = 0.1
N = 100

res = ilqg(x0, N, maxIter = 100, logg = T, mode = "control", rand = F)
vec_plot(res,xlim = c(-1,3), ylim = c(-1,3))

#plot(t(res$x), col=2, type='l', xlim = c(-1,3), ylim = c(-1,3), asp=1, xlab = NA, ylab = NA); grid()
#par(new=T)
#W = to_World(res$u)
#plot(t(W$x), col=3, type='l', xlim = c(-1,3), ylim = c(-1,3), asp=1, xlab = NA, ylab = NA); grid()
#par(new=T)



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
