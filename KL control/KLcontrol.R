x = 1
s=0.01; s2 = s^2
A=1.1
c =0.1
a = c/sqrt(c^2 + 2*s2)
u = function(y){dnorm(y,A*x[i], sd = s) * (1+100*exp(-y^2/c^2))/(1+100*a*exp(-(A^2*x[i]^2)/(c^2+2*s2)))}
p = function(y){dnorm(y,A*x[i], sd = s)}
z  = function(y){ 1 + 100*exp(-y^2/c^2) }

x = array(0,5); n = 4
x[4] = 1 # x0 = 10
for(k in 4:1){
    den =  c^2 + 2*s2*(A^(2*k) - 1)/(A^2 -1)
    #z[n-k+1] = 1 + ( c/sqrt(den) )*exp(-A^(2*k)*x[n-k+1]^2/den)
    k1 = k-1
    den1 = c^2 + 2*s2*(A^(2*k1) - 1)/(A^2 -1)
    u = function(y){dnorm(y,A*x[n-k+1], sd = s) * 
            (1 + 1000*(c/sqrt(den1))*exp(-A^(2*k1)*y^2/den1))/
            (1 + 1000*(c/sqrt(den))*exp(-A^(2*k)*x[n-k+1]^2/den) ) }
    p = function(y){dnorm(y,A*x[n-k+1], sd = s)}
    
    y = seq(-2,2,0.001)
    u_val = u(y)
    p_val = p(y)
    z_val = z(y)
    plot(y,p_val,type='l', ylim = c(0,6), col=2)
    par(new=T)
    plot(y,u_val,type='l', ylim = c(0,6), col=6)
    
    x[n-k1+1] = sample(y,1,replace = T,prob=u_val)
    
    # z(x)  = 1 + (a[k]/sqrt(d[k])) * exp((x-b)^2/ c^2)
    
    # N-k+1 
    a[k] = a0 * c0
    d[k] = d[k-1] + 2*(s^2)*(A^(2*k-2))
    
    
    kl = function(u_val,p_val){
        y = u_val
        id = which(p_val==0)
        y[id] = 0
        y[-id] = log(u_val[-id]/p_val[-id])
        sum(y*u_val)
    }
    
    