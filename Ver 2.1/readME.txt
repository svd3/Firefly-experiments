Version 2.1

Updated 11-04-2017

Linearize 2 updated!! Dimension matching and partial derivative correction

Added framework for state and control dependent noise in both process and observation

And also observation can be a function of state and control now

i.e. in linearization : y(k) = H x(k) + E u(k) ## important

both incorporated in 

Added:
 
file: 				function:

linearize2.R	———————————>	linearize() ## new implementation
new_kalman.R	———————————>	lqg_all()

Optimal_finite_time.R —— optimizes finite horizon problem


Deletes:

gen_trajectory.R
temp_run.R
kalman_lqg.R #### kalman_lqg2.R is better so we are using that





