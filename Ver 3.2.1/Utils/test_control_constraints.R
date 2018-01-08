u = 0
Vseq = 1000 * c(0.5,1,2,5,8,10,12,15,17.5,20,25,30,40,50,60,90,100,120,150,200,250,400,500)
Nopt = 0
for(k in 1:10){    
    i= 1
    cost = 0
    Nseq = seq(10,100, by=2)
    V = Vseq[k]
    res = list()
    for(N in Nseq){
        Qf = V*diag(2)
        Kx = V/100
        res[[i]] = ilqg(x0, N, maxIter = 100, eps=1e-2, noise = T, logg = F, mode = "control", rand = F)
        cost[i] = (res[[i]]$cost - Kx)/N
        print(i)
        i = i+1
    }
    Nopt[k] = Nseq[which.min(cost)]
    res_opt = res[[which.min(cost)]]
    vel[k] = mean(res_opt$u[1,])
    print(k)
}
u=0
Nopt = 0
for(k in 1:length(Vseq)){
    V = Vseq[k]
    Qf = V*diag(2)
    Nseq = seq(10,100,10) ## 10 iterations
    ans = optimal_N(Nseq, kx = V/100)
    Nseq = seq(ans[[1]]-10, ans[[1]]+10, 5) #5 here
    ans = optimal_N(Nseq, kx = V/100)
    Nseq = seq(ans[[1]]-5, ans[[1]]+5, 1) #10 here
    ans = optimal_N(Nseq, kx = V/100)
    res = ans$res
    u[k] = mean(sqrt(res$u[2,]^2 + res$u[1,]^2))
    print(k)
}