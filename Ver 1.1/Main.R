setwd("Documents/Acads/FireflyProject-master/Rcodes/Ver 1.1")
source("gen_trajectory.R")
# change if needed
Qf = 1000*diag(2) 
Q = array(0,c(2,2))
R = array(c(10,0,0,2), c(2,2))


## iterative ilqg
x0 = c(2,2)
dt = 0.1
N = 500

#res2 = gen_trajectory(x0,N,maxIter=10)
res = ilqg(x0, N, maxIter = 100, logg = T)
plot(t(res$x), col=2, type='l', xlim = c(-1,3), ylim = c(-1,3), asp=1, xlab = NA, ylab = NA); grid()
par(new=T)

W = to_World(res)
plot(t(W$x), col=3, type='l', xlim = c(-1,3), ylim = c(-1,3), asp=1, xlab = NA, ylab = NA); grid()
par(new=T)
vec_plot(W,xlim = c(-1,3), ylim = c(-1,3))
