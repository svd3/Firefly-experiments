j= j+1#trial number
tar = data$subs[[1]]$trial[[j]]$fly/100
X = data$subs[[1]]$trial[[j]]$trj
X[,1:2] = X[,1:2]/100

source("ilqg.R")
N = round(nrow(X)*0.012/0.1)
xt = c(tar,0,0,0)

res = ilqg(x0, N, maxIter = 100, eps=1e-5,logg = T, rand = F)
sim = simulation(x0,res)

vec_plot(t(X[,1:2]), (pi/2-X[,3]), xlim = c(-6,6), ylim = c(0,6), col=1)
points(tar[1],tar[2],pch=17, col=1)
par(new=T)
vec_plot(sim$x[1:2,], sim$x[3,], xlim = c(-6,6), ylim = c(0,6), col=2)


eqf = function(p){
    R = array(c(p[1],p[2],p[2],p[3]), c(2,2))
    vb =p[4]; wb = p[5]
    det(R)
}
func = function(p){
    R = array(c(p[1],p[2],p[2],p[3]), c(2,2))
    vb =p[4]; wb = p[5]
    res = ilqg(x0, N, maxIter = 100, eps=1e-5,logg = F, rand = F)
    #cat('.'); cat("\014");
    sim = simulation(x0,res)
    sum((X2 - t(sim$x[1:3,]))^2)
}