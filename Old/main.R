source("gen_iLQG.R")


M = 1
nx = 2
nu = 2 
ny = 2
nw = 2
nv = 2
N = 20;
d = list(nx = nx, nu = nu, ny = ny, nw = nw, nv = nv, N =N)
dt = 0.1
x0 = array(c(5,5),2)
C0 = array(0, c(nx,nw))
D0 = array(0, c(ny,nv))

R = array(c(10,0,0,2), c(2,2))
Qf = Q = diag(2)

res = list()
x_n = x_p = array(0, c(nx,N,M))
u_n = u_p = array(0, c(nu,N-1,M))
Lx = array(0, c(nu,nx,N-1,M))
lx = array(0, c(nu,N-1,M))
for(m in 1:M){
	print("main loop :: ")
    soln = gen_ilqg(f,Fn,g,G,h,l,x0,d)
    x_p = t(soln[[1]])
    u_p = t(soln[[2]])
    
}

plot(x_p, asp=1)


