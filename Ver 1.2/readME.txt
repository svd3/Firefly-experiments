Version 1.2

Updated 22-03-2017

Changed state augmentation

x~ = (x, u[k-1], 1)
u~ = u

Code working. Perfected noiseless implementation

LQR (own) implementation, no noise
Can be replaced by Kalman_lqg (but its slow) # checked and verified
#Both functions give same L matrix (verified)
So for noiseless case switch to lqr.R




