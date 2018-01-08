#velocity control # res_1
# ta -> infi => velocity control
tau_l = c(0,0.2,0.5,1,2,5,10)
resl = list()
for(i in 1:length(tau_l)){
    tau = tau_l[i]
    N = opt_time()[[1]]
    resl[[i]] = ilqg(x0, N, maxIter = 100, eps=1e-6, noise = T, logg = T, mode = "control", rand = F)
}

# velocity profiles
plot(c(res_1$u[2,],0), type='l', ylim = c(-0.35,0),lty =2, col=1, ylab = NA, yaxt = "n", xlab = "time", main = "ang. vel. profiles")
for(i in 1:length(tau_l)){
    res_t = res_l[[i]]
    par(new=T)
    plot(res_t$x[5,], type='l', ylim = c(-0.35,0), col=(i), yaxt="n", ylab = NA, xlab = NA)
}
axis(2,at= c(0,-0.1,-0.2,-0.3))
leg = c("vel. control")
for(i in 1:length(tau_l)){
    leg = c(leg, paste("tau =", tau_l[i]))
}
legend(x = 1, y= 1,legend = leg, lty=c(2,rep(1,7)), col = c(1,1:7), cex= 0.5)


vec_plot(res_1$x[1:2,], res_1$x[3,],xlim = c(0,2.2), ylim = c(0,2.2), col=1, asp=1, lty=2)
for(i in 1:length(tau_l)){
    res_t = res_l[[i]]
    par(new=T)
    vec_plot(res_t$x[1:2,], res_t$x[3,],xlim = c(0,2.2), ylim = c(0,2.2), col=i, asp=1)
}
legend(x = 1, y= 1,legend = leg, lty=c(2,rep(1,7)), col = c(1,1:7), cex= 0.5)
title("trajectory")


amax_l = c(0.1,0.5,1,2,5)
res_l3 = list()
for(i in 1:length(amax_l)){
    amax = amax_l[i]
    N = opt_time()[[1]]
    res_l3[[i]] = ilqg(x0, N, maxIter = 100, eps=1e-6, noise = T, logg = T, mode = "control", rand = F)
}

res_t = res_l3[[1]]
plot(res_t$u[1,], type='l', xlim = c(0,200),ylim = c(-1,1),lty =1, col=1, ylab = NA, yaxt = "n", xlab = "time")
for(i in 2:length(amax_l)){
    res_t = res_l3[[i]]
    par(new=T)
    plot(res_t$u[1,], type='l', xlim = c(0,200), ylim = c(-1,1), col=(i), yaxt="n", ylab = NA, xlab = NA)
}
axis(2)
