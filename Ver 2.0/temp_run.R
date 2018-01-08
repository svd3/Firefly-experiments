
ani.options(interval = 0.1)

animate  = function(){
    plot(res$x,type='n', col=2, xlim = c(0,3),ylim = c(0,3),asp=1, xlab =NA, ylab =NA, axes=F)
    box(); axis(1, at = c(0,1,2,3)); axis(2, at = c(0,1,2,3)); grid()
    for(k in 1:(N-1)){
        #Sys.sleep(0.1)
        plot(res$x,type='n', col=2, xlim = c(0,3),ylim = c(0,3),asp=1, xlab =NA, ylab =NA, axes=F)
        box(); axis(1, at = c(0,1,2,3)); axis(2, at = c(0,1,2,3)); grid()
        for(i in 1:k){
            arrows(res$x[1,i],res$x[2,i], res$x[1,i+1], res$x[2,i+1], length = 0.05, col=2)
        }
        #par(new=T);
    }
}
saveVideo(animate())


## last pt is target
pts = pts_gen(100,30,500)
pts = rbind(pts, x0)

for(k in 1:(N-1)){
    pts = imap(pts)
    x= pts[,1]; y = pts[,2]
    pts[,1] = x - w*y*dt
    pts[,2] = y - (v - w*x)*dt
    pts = map(pts)
    plot(pts, ylim = c(-0.4,0.2), xlim = c(-0.5, 0.5), pch = 19, xaxt ='n', yaxt ='n', xlab = NA, ylab = NA)
    points(x0[1],x0[2], col=2, pch=19)
    abline(h=0)
}
