library(pracma)

root = function(f,y,range,...){
    # f(x) - y = 0
    val =0; nr = 0; x0 =0 
    x = seq(range[1], range[2], by=diff(range)/100)
    p_val = f(x[1],...)-y
    for(i in 1:length(x)){
    	val[i] = f(x[i],...)-y
    	if((val[i]>=0 & p_val <0) | (val[i]<=0 & p_val >0)){
    		nr = nr + 1
    		x0[nr] = x[i]
    	}
    	p_val = val[i]
    }
    roots = 0
    for(i in 1:nr){
    	xi = x0[i]
	    converged = F
	    while(!converged){
	        f1 = jacobian(f,x0=xi,...)
	        x1 = xi - (f(xi,...)-y)/f1
	        xi = x1
	        if(abs(f(xi,...)-y) < 10e-12){
	            converged = T
	        }
	    }
	    roots[i] = xi
	}
	roots
}

root2 = function(f,y, x0,eps = 1e-12,...){
  converged = F
  while(!converged){
    f1 = jacobian(f,x0=x0,...)
    #print(f1)
    x1 = x0 - ginv(f1)%*%(f(...,x0)-y)
    x0 = x1
    if(sum(abs(f(...,x0)-y)) < eps){
      converged = T
    }
  }
  x0
}
