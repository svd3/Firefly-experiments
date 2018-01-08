#setwd("~/Documents/Acads/FireflyProject-master/Rcodes/Ver 3.2.1")

source("ilqg.R")

res = ilqg(x0, N, maxIter = 100, eps=1e-5,logg = T, rand = F)
sim = simulation(x0,res)
par(new=T)
vec_plot(sim$x[1:2,], sim$x[3,],xlim = c(0,2),ylim = c(0,2),col=2,asp=1)
points(tar[1],tar[2],pch=17)

for(i in 17:24){
  X = data$subs[[1]]$trial[[i]]$trj
  tar = data$subs[[1]]$trial[[i]]$fly
  
  vec_plot(t(X[,1:2]),(pi/2-X[,3]), col=i, xlim = c(-500,500), ylim = c(0,600), asp=1)
  points(x=tar[1],y=tar[2],pch=17, col=i)
  par(new=T)
}


func = function(p){
    #R <<- diag(p[1:2])+ 0.01*diag(2)
    gain <<- p[1:2]
    wm <<- p[3]
    res = ilqg(x0, N, maxIter = 100, eps=1e-5,logg = F, rand = F)
    print(paste("cost = ", res$cost))
    #cat('.'); cat("\014");
    sim = simulation(x0,res)
    #vec_plot(sim$x[1:2,], sim$x[3,],xlim = c(0,2),ylim = c(0,2),col=2,asp=1)
    ans = mean((sim$x-X0)^2)
    print(paste("value = ", ans))
    ans
}

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


