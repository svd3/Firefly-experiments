J=0
pb <- txtProgressBar(min = 0, max = 361, style = 3)
for(i in 1:dim(ab)[1]){
    R = makeR(ab[i,1],ab[i,2])
    x = create_trajectory(3)
    J[i] = sum((x - xo)^2)
    setTxtProgressBar(pb, i)
}