Version 1.1

Updated 06-03-2017

For EXPERIMENTATION
Added file gen_trajectory.R
Implements “iterative ilqg”

Code working

LQR (own) implementation, no noise
Can be replaced by Kalman_lqg (but its slow) # checked and verified
#Both functions give same L matrix (verified)
So for noiseless case switch to lqr.R

Interesting finds

The iterative ilqg works better than ilqg
Better optimal cost
Reduces control costs greatly
Curved trajectory


