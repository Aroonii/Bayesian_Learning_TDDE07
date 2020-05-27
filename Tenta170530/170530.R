library(Bessel)

#compute log density 
log.post = function(x, theta, v){ 
  log_like = sum(log(x/v)*(-(x^2+theta^2)) + log(besselI(x*theta/v, 0)))
return(log_like)
}

#log.posterior = function(theta = 1, data ){
#  sum(log.dens(x = data ,theta,v = 1))
#}


like = c()
riceData <- c(1.556, 1.861, 3.135, 1.311, 1.877, 0.622, 3.219, 0.768, 2.358, 2.056)
sequence = seq(0.1, 5, 0.1)
theta = 1
v = 1
prior = 1

# Plot the posterior over a grid of theta values
gridWidth <- 0.01
thetaGrid <- seq(0.01, 3, by = gridWidth)
logRicePostGrid <- rep(0,length(thetaGrid))
count <- 0
for (theta in thetaGrid){
  count <- count + 1
  logRicePostGrid[count] <-log.post(riceData, theta, v)
}
exp(logRicePostGrid)
plot(gridWidth, exp(logRicePostGrid))
