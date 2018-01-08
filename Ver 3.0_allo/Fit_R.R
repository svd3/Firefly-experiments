r = runif(10,1.5,4)
th = runif(10,40,140)*pi/180
target = matrix(c(r*cos(th),r*sin(th)), ncol=2)

R0 = array(c(1,-0.75,-0.75,1), c(2,2))

## make R
makeR = function(a,b){ 
    R = diag(2)
    R[1,1] = a; R[2,1] = R[1,2] = b/2;
    R
}

create_trajectory = function(n,progress= T){
    if(progress){
        pb <- txtProgressBar(min = 0, max = n, style = 3)
    }
    x = array(0,c(n*N,2))
    for(i in 1:n){
        x0 = target[i,]
        res = ilqg(x0, N, maxIter = 100, eps=1e-4, noise = T, logg = F, mode = "control", rand = F)
        #x1[[i]] = t(res$x)
        x[(N*i-N+1):(N*i),] = t(res$x) #
        if(progress){
            setTxtProgressBar(pb, i)
        }
    }
    x
}

myfun = function(beta){
    R <<- makeR(beta[1],beta[2])
    #x = create_trajectory(3,T)
    res = ilqg(x0, N, maxIter = 100, eps=1e-3, noise = T, logg = T, mode = "control", rand = F)
    sum((res$x - xo)^2)
}

op = optim(c(1,0),myfun, control = list(maxit=1),method = "BFGS")

v=0
for(j in 1:100){
    print(paste("iter =", j))
    #i = sample(1:10,1)
    x = create_trajectory(1)
    #J = sum((x - x1[[i]])^2)
    J = sum((x - xo)^2)
    Jl[j] = J
    plot(j, J, xlim = c(1,100), ylim = c(0,3), pch = 19,cex=0.5)
    if(j!=100){ par(new=T) }
    print(paste("J =" , round(J,6)))
    R = makeR(a+0.01,b)
    x = create_trajectory(1)
    ga = (sum((x - xo)^2) - J) #(sum((x - x1[[i]])^2) - J)
    R = makeR(a,b+0.01)
    x = create_trajectory(1)
    gb =  (sum((x - xo)^2) - J) #(sum((x - x1[[i]])^2) - J)

    v = 0.1*v + 0.1*c(ga,gb)
    ##
    #e = lr(ga,gb); print(paste("lr =", round(e,6)))
    a = a - v[1]#lr[j]*ga
    b = b - v[2]#lr[j]*gb
    R = makeR(a,b)
    print(R)
}