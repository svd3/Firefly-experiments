sigm = function(a){
    2*exp(a)/(1+exp(a)) - 1
    
}
bellFunc = function(x,c,s){
    exp(-s*sum((x-c)^2))
}

rbfFitting = function(x,y, k=10, scale, model = NA, centers = NA){
    if(!is.na(model)[1]){
        N = nrow(x)
        genX = x[sample(N,round(0.8*N)),]
        genX = genX + rnorm(length(genX),0,0.1)
        genY = rbfEval(genX,model)
        x = rbind(x,genX)
        y = matrix(c(y,genY))
    }
    N  <- dim(x)[1] # number of observations
    if(is.na(centers)[1]){
        repeat {
        km <- kmeans(x, k)  # let's cluster K centers out of the dataset
        if (min(km$size)>0) # only accept if there are no empty clusters
            break
    }
        mus = km$centers
    }
    else{
        mus = centers
    }
    nF = nrow(mus)
    x2 = matrix(0, N,nF)

    for(i in 1:nF){
        x2[,i] = apply(x, 1, bellFunc, c= mus[i,], s = scale)
    }
    ls = lm(y~x2)
    w = as.numeric(ls$coefficients)
    err = mean(ls$residuals^2)
    
    list(weights=w, centers=mus, scale=scale, mse = err)
}

phi = function(s){
    
}
rbfEval = function(x, rbfObj){
    centers = rbfObj$centers
    w = as.numeric(rbfObj$weights)
    scale = rbfObj$scale
    N = dim(x)[1]
    if(is.null(N)){
        N = 1
        x = matrix(x, ncol = length(x))
    }
    nF = nrow(centers)
    pred <- rep(w[1],N)  # we need to init to a value, so let's start with the bias
    
    nF = nrow(centers)
    x2 = matrix(0, N,nF)
    
    for(i in 1:nF){
        x2[,i] = apply(x, 1, bellFunc, c= centers[i,], s = scale)
    }
    
    pred = pred + x2%*%w[-1]
    
    pred
}

# state space
data = mvrnorm(10000,rep(0,5), diag(c(5,5,3,1,1)))
y = rep(0,10000)



s0 = c(2,2,pi/2)
#episode
s=s0
reward = 0
X = matrix(0, nrow=10000, ncol=5)
for(i in 1:nrow(X)){
    if(i%%100==0)
        print(i)
    a = policy(s)
    X[i,] = c(s,a)
    reward[i] = R(s,a)
    s = nextState(s,a)
    if(isGoal(s)){
        print("Goal")
    }
}

y1 = apply(X,1,rbfEval,rbfObj = model)
y2 = cumsum(rev(reward))
y2 = 0.5*y1 + 0.5*y2
model = rbfFitting(X,y2,100,2, centers = model$centers)



policy = function(s){
    #op = optim(c(0,0),opf, s=s)
    a = mvrnorm(100,c(0,0),diag(2))
    temp = matrix(0, 100,5)
    temp[1:100,1:3] = rep(s,each=100)
    temp[,4:5] = a
    tempy = rbfEval(temp,model)
    a[which.max(tempy),]
}

nextState = function(s, a) {
    s1 = c(a[1]*cos(s[3]), a[1]*sin(s[3]), a[2]) * 0.1 + s + rnorm(3,0,0.01)
    s1
}

opf = function(a,s){
    -rbfEval(c(s,a),model)
}

isGoal = function(s){
    norm(s[1:2],"2") <= 0.5
}
R = function(s,a){
    if(isGoal(s))
        R = 10
    else
        R = -1
    R
}


    



for(i in 1:1000){
    x = matrix(c(runif(250,0,2), runif(250,-1,1)), 250,2)
    y = (x[,1]+x[,2])/(1+x[,1]^2 + x[,2]^2)
    model = rbfFitting(x,y,10,1,model)
}
ytest2 = rbfEval(xtest, model)

plot(ytact,type='l', ylim = c(-0.6,0.7))
par(new=T)
plot(ytest2,col=3,type='l', ylim = c(-0.6,0.7))
ytact = (xtest[,1]+xtest[,2])/(1+xtest[,1]^2 + xtest[,2]^2)
