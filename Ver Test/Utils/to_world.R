## converting to World coordinates
## updated 24-08-2017
## vec_plot updated to be more generic
vec_plot = function(pos, ang, col =3,...){
    N = length(ang)
    plot(t(pos), type='l', xlab = NA, ylab = NA, col = col, ...);
    abline(h = seq(-10,10,1), v = seq(-10,10,1), lty=3, col=8)
    if(N<20){
        time_seq = 1:N
    } else{
        time_seq = c(seq(1,N, by=round(N/20)),N)
    }
    for(k in time_seq){
        x = pos[,k] #W$x[,k]
        th = ang[k] # W$heading[k]
        quiver(x[1],x[2], cos(th), sin(th),col=col,...)
    }
    points(x=pos[1,N],y=pos[2,N], col=col, pch=19)
}

rotMat = function(ang){
    c = cos(ang); s = sin(ang)
    return (matrix(c(c, -s, s, c), c(2,2)))
}

### Obsolete
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