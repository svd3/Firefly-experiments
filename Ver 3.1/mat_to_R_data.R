ans =T
subs = list()
k = 1
for(i in seq(1,18,3)){
    n = length(dat$subjects[[i]])/2
    subs[[k]] = list()
    subs[[k]]$name = paste("subject",k)
    subs[[k]]$nTrials = n
    subs[[k]]$j.speed = c(v.max = dat$subjects[[i]][,,j]$prs[[2]][1,1], w.max = dat$subjects[[i]][,,j]$prs[[2]][1,2])
    subs[[k]]$trial = list()
    for(j in 1:n){
        subs[[k]]$trial[[j]] = list()
        #ans = ans & (dat$subjects[[i]][,,j]$prs[[2]][,1] == 200) & (dat$subjects[[i]][,,j]$prs[[2]][,2] == 90)
        xfly = dat$subjects[[i]][,,j]$coord[[4]]
        yfly = dat$subjects[[i]][,,j]$coord[[5]]
        fly  = c(x.fly = xfly,y.fly = yfly)
        x = dat$subjects[[i]][,,j]$coord[[1]]
        y = dat$subjects[[i]][,,j]$coord[[2]]
        phi = dat$subjects[[i]][,,j]$coord[[3]]
        tr  = matrix(cbind(x,y,phi), ncol=3); colnames(tr) = c("x","y","phi")
        subs[[k]]$trial[[j]]$fly = fly
        subs[[k]]$trial[[j]]$trj = tr
        subs[[k]]$trial[[j]]$floor.dens = as.numeric(dat$subjects[[i]][,,j]$prs[[1]])
    }
    j=1
    k = k+1
}