# find optimal time
costps = function(N,kx=50){
    N = round(N)
    res = ilqg(x0, N, maxIter = 100, eps=1e-5,logg = F, rand = F)
    print(".")
    (res$cost - kx)/(N*dt)
}

###########  function
optimal_N = function(Nseq, kx = 50, dt = 1){
    dt <<- dt
    i=1
    cost=  0
    res = list()
    for(N in Nseq){
        res[[i]] = ilqg(x0, N, maxIter = 100, eps=1e-6, logg = F, rand = F)
        cost[i] = (res[[i]]$cost)# - kx)/N## rate
        #print(i*100/length(Nseq))
        i = i+1
    }
    imin = which.min(cost)
    #print("fin.")
    list(Nopt = Nseq[imin], cost, res = res[[imin]])
}


## heirachical approach
opt_time_old = function(){
    Nseq = seq(10,100,10) ## 10 iterations
    ans = optimal_N(Nseq, kx=50,dt=0.1)
    
    Nseq = seq(ans[[1]]-9, ans[[1]]+9, 5) #4 here
    ans = optimal_N(Nseq, kx = 50, dt=0.1)
    
    Nseq = seq(ans[[1]]-4, ans[[1]]+4, 1) #9 here
    ans = optimal_N(Nseq, kx = 50, dt=0.1)
    ans
}

opt_time = function(){
    # dt = 1
    Nseq = seq(2,20,5) # 20 30 40 50 60 70 80 90 100
    ans = optimal_N(Nseq,kx=50,dt=1)
    ans[[1]] = ans[[1]]-1
    # dt = 0.5 
    Nseq = seq((ans[[1]]*2 -4),(ans[[1]]*2 +4),2) # 50 55 60 65 70
    ans = optimal_N(Nseq,kx=50,dt=0.5)
    ans[[1]] = ans[[1]] -1
    # dt = 0.25
    Nseq = seq(ans[[1]]*2-2, ans[[1]]*2+2, 2) #5 here
    ans = optimal_N(Nseq, kx = 50, dt = 0.25)
    ans[[1]] = ans[[1]] -1
    
    temp = round(ans[[1]]*2.5)
    Nseq = seq(temp-3, temp+3, 1) #6 here
    ans = optimal_N(Nseq, kx = 50, dt = 0.1)
    ans
}