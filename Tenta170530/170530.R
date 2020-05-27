library(Bessel)

#compute log posterior 
# Since we don't have any information about the prior we assuem it to be constant. Ex 1 => log(constant) = constant
# when normalazing this all vlues will be affected equally hence it does not make any difference
log.post = function(theta, x, v){ 
  log_like = sum(log(x/v)*(-(x^2+theta^2)) + log(besselI(x*theta/v, 0))) + 2
return(log_like)
}

riceData <- c(1.556, 1.861, 3.135, 1.311, 1.877, 0.622, 3.219, 0.768, 2.358, 2.056)
sequence = seq(0.1, 5, 0.1)
theta = 1
v = 1
prior = 1

# Plot the posterior over a grid of theta values. Gives us the probability depending on different theta values 
gridWidth <- 0.01
thetaGrid <- seq(0.01, 3, by = gridWidth)
logRicePostGrid <- rep(0,length(thetaGrid))
count <- 0
for (theta in thetaGrid){
  count <- count + 1
  logRicePostGrid[count] <-log.post(theta, riceData, v)
}

posterior = exp(logRicePostGrid)
posterior_normalized = (1/gridWidth)*posterior/(sum(posterior))

plot(posterior_normalized)


#b) Normal approximation of the posterior distr. of theta 

# Different starting values. Ideally, any random starting value gives you the same optimum (i.e. optimum is unique)
initVal <- c(1); 

# function which optmizes over expression log.posterior with respect to its first argument (betas). 
# returns optimal values for beta (mode), and hessian in the mode
OptParams<-optim(initVal,log.post,gr=NULL,riceData, v ,method=c("L-BFGS-B"),lower = 0, control=list(fnscale=-1),hessian=TRUE)

my_posterior = OptParams$par
# takes negative so that the posterior can be approx. as normal
# J = -second derivate evaluated in theta_hat
hessian_posterior = -OptParams$hessian

# take inverse for using it in the formula
hessian_posterior = (solve(hessian_posterior))
sigma = sqrt(hessian_posterior[1,1])

#Draw samples from Betas posterior distribution. 
set.seed(12345)

#We now have a posterior distributin which is ~N(my_posterior, sigma_square). We want to do the same thing as before
# plot the pdf over a grid of values => Since we have it in a good form we can do it immedeatly and directly get
# the normalized values 

approximate_density_distribution = dnorm(thetaGrid, mean = my_posterior, sd = sigma)
plot(thetaGrid, posterior_normalized)
lines(thetaGrid, approximate_density_distribution, col = "red", type ="b")


#C) 
#1. Explain on paper how the predictive distribution for a new observation is ocmputed by integration. Only
#The general dormula
#2. Compute the predictive distribution for a new observation by simulation. USe the approximate posterior from b
# simulator for r_rice is given by exam file

rRice <-function(n = 1, theta = 1, psi = 1){
  x <- rnorm(n = n, mean = 0, sd = sqrt(psi))
  y <- rnorm(n = n, mean = theta, sd = sqrt(psi))
  return(sqrt(x^2+y^2))
}
#Draw from the posterior distribution and use these draws to plug in and make draw from test data
outcome =c()
set.seed(12345)
for (i in 1:1000){
  #nake draw from posterior distr.
  theta_draw = rnorm(n = 1, mean = my_posterior, sd = sigma)
  #plug in draw and make fraw from predictions distr. 
  rice = rRice(1,theta_draw,1)
  outcome = append(outcome, rice)
}

#Histrogram of the predictive 
hist(outcome, main="Predictive distribution x_tilde", xlab="rice", ylab="Acumulation of rice out of 1000 draws")


#Assignemnt 2

# eBay bids data
load(file = 'bids.RData')    # Loading the vector 'bids' into workspace
bids = bids
bidsCounts <- table(bids)  # data2Counts is a frequency table of counts.
xGrid <- seq(min(bids),max(bids))  # A grid used as input to GibbsMixPois.R over which the mixture density is evaluated.

#a
#Copute the posterior distribution for theta and plot it => the posterior pdf
n = auctions = length(bids)
alpha = 1
beta = 1
amount_bids = sum(bids)

#Posterior parameters for gamma distriution
alpha_posterior = amount_bids + alpha
beta_posterior = n + beta

#Plot the posterior distribution for theta (gamma funktion) 
# => draws from gamma density distribution over a seequence of different values
gamma_grid = seq(3.4,3.9, length = 1000)
plot(gamma_grid, dgamma(gamma_grid, alpha_posterior, beta_posterior), xlab = "Theta", ylab = "Density")

#b
#Does the poisson model describe the distribution of the data well ?
#=> Model evaluation => Make draws from the posterior of theta. Use these to plug in
# to the distribution of the data (in this casepossion) and make new replica draws
data_density = density(bids)
plot(data_density)

xGrid = seq(min(bids),max(bids))

#Manually create the density of the data
#Normalize
data_dens = bidsCounts/sum(bidsCounts)
plot(xGrid, data_dens, type = 'l')

#Make repliica draws using draws from posterior

set.seed(12345)
theta_draw = rgamma(1000, alpha_posterior, beta_posterior)
#For each draw theta we want to compute the possion distribution (distribution for
# x values)

density_mean = rep(0,12)
for (i in 1:n){
#Compute poission distribution over Xgrid = över de x värden vi vill distribuera 
# proba
density_mean = density_mean + dpois(xGrid, theta_draw[i])
}
#Make a T function (condition) that we can compare with the original data. in this 
#case we compute the mean density for each grid value
density_mean = density_mean/n

plot(xGrid, data_dens, type = 'l' )

lines(density_mean)







