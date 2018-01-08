## converting to World coordinates
to_World = function(res, init = F){
    ## x = Frame origin wrt World (Our position wrt World)
    ## O_r = World origin wrt relative frame
    
    ## res$x is the translation of target in self-relative coordinate
    xt_r = res$x ## target wrt relative frame
    u = res$u
    if(init){
        xt_r = res$init_x
        u = res$init
    }
    
    N = dim(res$x)[2]
    dt = res$dt
    ## first calculate change of world origin in self-relative frame
    O_r = array(0, c(nx, N)) # initial world origin same as frame origin
    
    ## calculate shift in origin when control is applied
    for(k in 1:(N-1)){
        O_r[,k+1] = O_r[,k] + f(O_r[,k], u[,k])*dt
    }
    
    ## now to find theta or cos(theta) and sin(theta)
    A = ginv(matrix(c(x0,-x0[2],x0[1]), ncol=2))
    
    CS = A %*% (xt_r - O_r)  ## CS[,k] = [cos(theta), sin(theta)] @ time k
    ## normalize CS (not must)
    
    x = array(0, c(nx,N))
    th = 0
    for(k in 1:N){
        CS[,k] = CS[,k]/Norm(CS[,k])  ## normalizing (not needed)
        th[k] = atan2(CS[2,k],CS[1,k]) ## atan2(y,x)
        R = matrix(c(-CS[,k],CS[2,k], -CS[1,k]),ncol=2) # R = [-c, s; -s, -c]
        x[,k] = ginv(R)%*%O_r[,k]
    }
    th = pi/2 - th
    list(x=x, heading = th, CS = CS)
}

vec_plot = function(W, xlim, ylim){
    #rx = range(W$x[1,]); ry = range(W$x[2,])
    N = dim(W$x)[2]
    plot(t(W$x), col=3, type='l', xlim = xlim, ylim = ylim, asp=1, xlab = NA, ylab = NA); grid()
    for(k in c(seq(1,N, by=round(N/20)),N)){
        x = W$x[,k]
        th = W$heading[k] ## th = pi/2 - th ?
        quiver(x[1],x[2], cos(th), sin(th))
        #quiver(x[1],x[2], x[1]+0.1*W$CS[2,k], x[2]+0.1*W$CS[1,k])
    }
}