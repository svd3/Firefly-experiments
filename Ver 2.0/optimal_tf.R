# find optimal time

###########  function
optimal_N = function(Nseq, kx = 50){
    i=1
    cost=  0
    res = list()
    for(N in Nseq){
        res[[i]] = ilqg(x0, N, maxIter = 100, eps=1e-2, noise = T, logg = F, mode = "control", rand = F)
        cost[i] = (res[[i]]$cost - kx)/N
        #print(i)
        i = i+1
    }
    imin = which.min(cost)
    print("fin.")
    list(Nopt = Nseq[imin], cost, res = res[[imin]])
}


## heirachical approach
Nseq = seq(10,100,10) ## 10 iterations
ans = optimal_N(Nseq, kx=50)

Nseq = seq(ans[[1]]-10, ans[[1]]+10, 5) #5 here
ans = optimal_N(Nseq, kx = 50)

Nseq = seq(ans[[1]]-5, ans[[1]]+5, 1) #10 here
ans = optimal_N(Nseq, kx = 50)

## total only 25 needed as opposed to 90!!  can be made smarter
########## 

 
tf = Nseq*dt
par(bg = grey(0.98))
par(mar = c(3.5,3.75,1,1))
plot(tf, cost,type='l', col=2, lwd=2, xlab = "", ylab = "", ylim = c(0,115))
axis(1, at = c(1,3,5,7,9))
title(xlab = TeX("$t_f$"), ylab =TeX("$J^*(t_f)$"), line=2.5)

grid()



par(bg = grey(0.98))
par(mar = c(3.5,3.75,1,1))
kx = c(50,60,75,90,100)
rate = list()
for(i in 1:5){
    rate[[i]] = (cost - kx[i])/tf
    plot(tf, rate[[i]], type='l', ylim = c(-30,30), yaxt="n", col=i, lwd = 1.25, xlab = "", ylab = "")
    #abline(v = tf[which.min(rate[[i]])], lty=2, col=i)
    t_opt = tf[which.min(rate[[i]])]
    J =min(rate[[i]])
    lines(x = rep(t_opt,2), y = c(J,0), col=i, lty = 2, lwd =1.25)
    points(t_opt,J, pch=19, col=i, cex = 0.7)
    #grid()
    if(i!=5){
        par(new=T)
    }
}
abline(h=0, lty=2)
axis(1, at = c(1,3,5,7,9))
axis(2, at = c(0,15,30, -15,-30))
title(xlab = TeX("$t_f$"), ylab =TeX("$J_r^*(t_f)$"), line=2.5)
legend(x=8,y=29, c(TeX("$k_x = 50$"),TeX("$k_x = 60$"),TeX("$k_x = 75$"),TeX("$k_x = 90$"),TeX("$k_x = 100$") ), col= c(1,2,3,4,5), lwd=rep(1.25,5), cex = 0.8)


x = seq(-2,2,0.01)
y = 1000*x^2

par(bg = grey(0.98))
par(mar = c(3.5,3.5,1,1))
plot(x,y, type='l', lwd = 1.25, col=2, ylim = c(-60,60), xlim = c(-1,1), xlab ="", ylab ="", yaxt ="n", xaxt="n")
abline(h=0, lty=2)
axis(1, at = c(-1,0,1))
axis(2,at = c(-50,-25,0,25,50))

arrows(0,0, 0,-50, 0.05, 20, col=3)
arrows(0,-50, 0,0, 0.05, 20, col=3)
arrows(0,5, 0.22,5, 0.05, 20, col=3)
arrows(0.22,5, 0,5, 0.05, 20, col=3)
text(0.05,-25,TeX("$k_x$"), cex = 0.9)
text(0.1,9,TeX("$\\Delta x$"), cex = 0.8)
title(xlab = TeX("$x$"), ylab =TeX("$h(x)$"), line=2)


Q = c(500,1000,1250,1500,2000)

par(bg = grey(0.98))
par(mar = c(3.5,3.5,1,1))
for(i in 1:5){
    y = Q[i]*(x^2 - 0.2^2)
    plot(x,y, type='l', lwd = 1.25, col=i, ylim = c(-80,60), xlim = c(-1,1), xlab ="", ylab ="", yaxt ="n", xaxt="n")
    if(i!=5){
        par(new=T)
    }
}
abline(h=0, lty=2)
axis(1, at = c(-1,0,1))
axis(2,at = c(-50,-25,0,25,50))
title(ylab = TeX("$h(x_{t_f})$"), xlab =TeX("$x_{t_f}$"), line=2)
legend(x=0.6,y=-40, c(TeX("$k_x = 20$"),TeX("$k_x = 40$"),TeX("$k_x = 50$"),TeX("$k_x = 60$"),TeX("$k_x = 80$") ), col= c(1,2,3,4,5), lwd=rep(1.25,5), cex = 0.8)

par(bg = grey(0.98))
par(mar = c(3.5,3.5,1,1))
plot(Vseq, u, type='o', col = 2, lwd =1.5, pch = 19, cex = 0.4, ylim = c(0.5,1), xlab = "", ylab ="", yaxt = "n", xaxt = "n")
abline(h=1, lty=2)
axis(1, c(1000,1e+5,2e+5,3e+5,4e+5,5e+5))
axis(2,at = c(0.6,0.8,1))
#title(ylab = "|| u ||", xlab = "Q or final cost", line=2.5)
#title(ylab = TeX("$\\lVert u \\rVert$"), xlab =TeX("$Q$"), line=2)
title(ylab = "|| u ||", xlab = TeX("$k_x$"), line=2.5)




par(mar = c(3.5,3.5,1,1))
r = runif(200,3,6); th = runif(200,50*pi/180,140*pi/180)
for(i in 1:200){ 
    x0 = c(r[i]*cos(th[i]), r[i]*sin(th[i]))
    #x0 = c(runif(1,-3,3), runif(1,1,3))
    N = 100
    res = ilqg(x0, N, maxIter = 100, eps=1e-2, noise = T, logg = T, mode = "control", rand = F)
    
    #res = ans$res
    W = to_World(res$u)
    #vec_plot(W,xlim = c(-1,3), ylim = c(-1,3), col=i)
    plot(t(W$x), col=grey(0.6), type='l', lwd=1.5, xlim = c(-6,6), ylim = c(0,6), xlab = "", ylab = "", yaxt = "n", xaxt="n");
    points(x0[1], x0[2], pch=19, col=grey(0.2), cex = 0.6)
    if(i!=200){
        par(new=T)
    }
}

