## converting to World coordinates
## updated 10-05-2017
## vec_plot updated to be more generic

to_World = function(u, init_pos = c(0,0,pi/2)){
    N = dim(u)[2] + 1
    pos = array(0,c(3,N))
    pos[,1] = init_pos
    for(k in 1:(N-1)){
        pos[1:2,k+1] = pos[1:2,k] + u[1,k]*dt*c(cos(pos[3,k]), sin(pos[3,k]))
        pos[3,k+1] = pos[3,k] - u[2,k]*dt
    }
    rownames(pos) = c("x","y","theta")
    list(x = pos[1:2,], heading = pos[3,])
}

vec_plot = function(pos, ang, xlim, ylim, col=3){
    #rx = range(W$x[1,]); ry = range(W$x[2,])
    #N = dim(W$x)[2]
    N = length(ang)
    #plot(t(W$x), col=col, type='l', xlim = xlim, ylim = ylim, asp=1, xlab = NA, ylab = NA); grid()
    plot(t(pos), col=col, type='l', xlim = xlim, ylim = ylim, asp=1, xlab = NA, ylab = NA); grid()
    if(N<20){
        time_seq = 1:N
    } else{
        time_seq = c(seq(1,N, by=round(N/20)),N)
    }
    for(k in time_seq){
        x = pos[,k] #W$x[,k]
        th = ang[k] # W$heading[k]
        quiver(x[1],x[2], cos(th), sin(th), col=col)
    }
}