# z(x)  = 1 + (g[k]/sqrt(d[k])) * exp(- (x-b)^2/ c^2)
# exp(- a[k] * (x-b)^2 / d[k])

z = function(x,k){
    sum = 1
    for(i in 1:M){
        sum = sum + g[k,i] * exp(- (x - b[k])^2 / c2[k,i])
    }
    sum
}

u = function(x,xp,k){
    dnorm(x, A[k-1] * xp, sd = sig) * z(x,k)/z(xp,k-1)
}

p = function(x, xp, k){
    dnorm(y,A[k-1] * xp, sd = sig)
}

# SYSTEM
x = array(0,N)
u_l = array(0,N-1)
A = array(0.98,N-1)
B = array(1,N-1) ## x[k+1] = A[k]*x[k] + B[k]*u[k] + e

sig = 0.1 ## noise
b0 = 2 ## target

## final cost  =  -M log(1 + w0 exp(- (x-b0)^2 / c0^2) )
c0 = 0.1; w0 = 1  ## cost function parameters

# z(x,k) =  1 + { (w/sqrt(d[k])) exp (- a[k]*(x-b[k])^2 / d[k]) }
# d[k] / a[k] = c2[k]
# g[k] = w / sqrt(d[k])
# z(x,k) =  1 + { g[k] exp (- (x-b[k])^2 / c2[k]) }
# w = nCr * c0 / sqrt(r)

## declaration of variables for desirability functions
a = b = array(0, N)
c2 = d = g = array(0, c(N,M))

## initializations for desirability functions
w = w0^(1:M); # w = array(1,M); 
for(m in 1:(M/2)){
    # nCm
    w[m] = w[M-m] = factorial(M)/ (factorial(m)*factorial(M-m))
}
w = w*c0/sqrt(1:M)

# z(x,N) = 1 + sum_m { nCm * (w0^m) * exp(- (x-b0)^2 / (c0^2 / m ) )   }
a[N] = 1
b[N] = b0
d[N,] = c0^2/(1:M)
c2[N,] = d[N,]/a[N]
g[N,] = w / sqrt(d[N,])

### BACKWARD PASS
for(k in (N-1):1){ 
    d[k,] = d[k+1,] + 2*(sig^2)*(a[k+1]^2) #(A^(2*k-2))
    a[k] = a[k+1]*(A[k]^2)
    b[k] = b[k+1]/A[k]  ## check A[k] = 0 condition!!
    c2[k,] = d[k,]/a[k]
    g[k,] = w/sqrt(d[k])
}

y = seq(-8,8,0.001) ## for sampler

### FORWARD PASS
for(k in 1:(N-1)){
    u_val = u(y,x[k],k+1)
    p_val = p(y,x[k],k+1)
    z_val = z(y,k+1)
    if(F){
        plot(y,z_val,type='l', col=2)
        #par(new=T)
        plot(y,p_val,type='l', ylim = c(0,6), col=1, xlim = c(-1,5))
        par(new=T)
        plot(y,u_val,type='l', ylim = c(0,6), col=2, xlim = c(-1,5))
    }
    x[k+1] = mean(sample(y,10000,replace = T, prob = u_val))
}
