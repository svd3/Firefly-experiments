#setwd("~/Documents/Acads/FireflyProject-master/Rcodes/Ver 3.0 test")

## We'll have to implement dy(k) = H(k)*dx(k) + noise
## dx^(k+1) = A(k)*dx^(k) + B(k)*u(k) + K(k)*(dy(k) - H(k)*x^(k))

source("ilqg.R")


## iterative ilqg
res = ilqg(x0, N, maxIter = 100, eps=1e-6, noise = T, logg = T, mode = "control", rand = F)
par(new=T)
vec_plot(res$x[1:2,], res$x[3,],xlim = c(-0.5,3), ylim = c(-0.5,3), col=4, asp=1)


plot(c(0,res$u[1,]), type='l', ylim = c(0,0.5), col=2, ylab = NA)
par(new=T)
plot(res$x[4,], type='l', ylim = c(0,0.5), col=1, ylab = NA)
par(new=T)
plot(res_a$x[4,], type='l', ylim = c(0,0.5), col=3, ylab = NA)
par(new=T)
plot(res_t$x[4,], type='l', ylim = c(0,0.5), col=3, ylab = NA)




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
