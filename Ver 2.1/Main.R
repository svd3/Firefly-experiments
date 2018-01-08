setwd("Documents/Acads/FireflyProject-master/Rcodes/Ver 2.1")

## We'll have to implement dy(k) = H(k)*dx(k) + noise
## dx^(k+1) = A(k)*dx^(k) + B(k)*u(k) + K(k)*(dy(k) - H(k)*x^(k))

source("ilqg.R")

## iterative ilqg
x0 = c(2,2)
dt = 1
N = 20

res = ilqg(x0, N, maxIter = 100, eps=1e-6, noise = T, logg = T, mode = "control", rand = F)
#plot(t(res$x), col=2, type='l', xlim = c(-1,3), ylim = c(-1,3), asp=1, xlab = NA, ylab = NA); grid()
par(new=T)
W = to_World(res$u)
vec_plot(W,xlim = c(-1,3), ylim = c(-1,3), col=3)


## Cost per time

cost_ps = 0
k = 1
kx = 0 ## important parameter relates to juice drop rewards
dt = 1
range_N = seq(5,30,by=1)
for(N in range_N){
    res = ilqg(x0, N, maxIter = 100,  eps=1e-2, noise = F, logg = T, mode = "control")
    cost_ps[k] = (res$cost - kx)/N
    k = k+1
}

plot(range_N, cost_ps, type='o')

ku = exp(mean(log(cost_ps)+ 2*log(range_N))) #appx. and kx = 0


###### R off diagonal elements matter
e = c(-0.99,seq(-0.9,0.9,0.1),0.99)
for(i in 1:length(e)){
    R[1,2] = R[2,1] = e[i]
    res = ilqg(x0, N, maxIter = 100, eps=1e-6, noise = T, logg = F, mode = "control", rand = F)
    print(i*100/length(e))
    W = to_World(res$u)
    vec_plot(W,xlim = c(-1,3), ylim = c(-1,3), col=i)
    if(i != length(e)){
        par(new=T)
    }
}
