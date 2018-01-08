soln = run(100)
control_cost(t(init_u))
control_cost(t(soln$u))

state_cost(t(init_x))
state_cost(t(soln$x))


cost_ps = 0;
i=1;
for(N in 5:20){
    res = gen_trajectory(x0,N);
    #X = run(100,F,F)
    cost_ps[i] = res$cost/N; 
    i=i+1 
}
