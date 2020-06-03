library(rstan)


t = seq(1,200, 1)
my = 10
sigma_square = 2
phi = -0.95

time_points = c()
set.seed(12345)

AR_process = function(t, phi){
  time_points = c()

for (i in 1:length(t)){
  if(i == 1){
    time_points[i] = my
    next
  }
  error_term = rnorm(1, mean = 0, sd = sqrt(sigma_square))
  
  time_points[i] = my + phi*(time_points[i-1] - my) + error_term
}
return(time_points)
}


plot(t, AR_process(t, 0.5), main = "phi = 0.5", ylab = "Xt", type = 'l')

#1b)

time_serie_0.95 = AR_process(t, phi = 0.95)
time_serie_0.3 = AR_process(t, phi = 0.3)


stanModel = ' data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real my;
  real phi;
  real<lower=0> sigma;
}

model {
  for (n in 2:N)
    y[n] ~ normal(my + phi *( y[n-1] - my), sigma);
}'

burnin = 1000
niter = 2000

# for phi = 0.95
data = list(N=length(t), y=time_serie_0.95)

fit_0.95 = stan(model_code=stanModel,data=data,
           warmup=burnin,iter=niter,chains=1)

# Print the fitted model and can se the numebr of efficient draws
print(fit_0.95,digits_summary=3)

# for phi = 0.3
data = list(N=length(t), y=time_serie_0.3)
fit_0.3 = stan(model_code=stanModel,data=data,
           warmup=burnin,iter=niter,chains=4)

# Print the fitted model and can se the nr of efficient draws
print(fit_0.3,digits_summary=3)



# Extract posterior samples
postDraws_0.95 <- extract(fit_0.95)
postDraws_0.3 <- extract(fit_0.3)

par(mfrow = c(1,1))
plot(postDraws_0.95$my,type="l",ylab="Phi", main="Sample convergence for phi when phi = 0.95")
plot(postDraws_0.3$my,type="l",ylab="Phi", main="Sample convergence for phi when phi = 0.3")

plot(postDraws_0.95$phi,type="l",ylab="My", main="Sample convergence for my when phi = 0.95")
plot(postDraws_0.3$phi,type="l",ylab="My", main="Sample convergence for my when phi = 0.3")


plot(postDraws_0.3$phi, postDraws_0.3$my, type="l",xlab = "phi", ylab="My", main="Sample convergence for my when phi = 0.3")
plot(postDraws_0.95$phi, postDraws_0.95$my, type="l",xlab = "phi", ylab="My", main="Sample convergence for my when phi = 0.95")

cat("Mean_my (0.3)", mean(postDraws_0.3$my))
cat("Mean_my (0.95)", mean(postDraws_0.95$my))
cat("Mean_my (0.3)", mean(postDraws_0.3$phi))
cat("Mean_my (0.95)", mean(postDraws_0.95$phi))


# Assignment 3

c_data = read.table('campy.dat', header = TRUE)



C_DataRStan<-
  list(N = nrow(c_data),
       c = c_data$c) 

# The poissons model
C_Model = '
data {
  int<lower=0> N;
  int c[N];
}



parameters {
  real alpha;
  real <lower=-1, upper= 1> beta;
  real<lower=0> sigma;
  real Xt[N];
}

model {

  // Prior
  for (n in 2:N){
    Xt[n] ~ normal(alpha + beta * Xt[n-1], sigma);
  }
  // Model/likelihood
  for(n in 1:N){
  
  c[n] ~ poisson(exp(Xt[n]));

  }                 
}'

fit_C_Model<-stan(model_code=C_Model,
               data=C_DataRStan,
               warmup=1000,
               iter=2000,
               chains=4)

print(fit_C_Model,digits=4)

# Plot some results
res<-extract(fit_C_Model)


poster_means_vector = fit_C_Model@.MISC[["summary"]][["msd"]][-c(1,2,3)]

cred_lower = fit_C_Model@.MISC[["summary"]][["c_quan"]][,1,1][-c(1,2,3,144)]
cred_upper = fit_C_Model@.MISC[["summary"]][["c_quan"]][,9,1][-c(1,2,3,144)]

N = nrow(c_data)
xGrid = seq(1,140,1)
plot(xGrid, c_data$c, col = "red", main = "Possion model & Data", ylab = "nr of Cases", xlab = "Draw")
lines(exp(poster_means_vector[1:N]), ylim = c(0,60), col = "black", type = 'l')
lines(exp(cred_lower), col = "blue")
lines(exp(cred_upper), col = "green")

legend('topleft',legend = c('Data', 'Mean','Lower','Upper'), col = c('red', 'black', 'blue', 'green'), lwd=2)





##2d)



# The poissons model
C_Model = '
data {
  int<lower=0> N;
  int x[N];
}



parameters {
  real alpha;
  real <lower=-1, upper= 1> beta;
  real<lower=0> sigma;
  real Xt[N];
}

model {

  // Prior for sigma
  sigma~inv_gamma(1,1)
  
  // Prior
  for (n in 2:N){
    Xt[n] ~ normal(alpha + beta * Xt[n-1], sigma);
  }
  // Model/likelihood
  for(n in 1:N){
  
  c[n] ~ poisson(exp(Xt[n]));

  }                 
}'

fit_C_Model<-stan(model_code=C_Model,
                  data=C_DataRStan,
                  warmup=1000,
                  iter=2000,
                  chains=4)

print(fit_C_Model,digits=4)

# Plot some results
res<-extract(fit_C_Model)


poster_means_vector = fit_C_Model@.MISC[["summary"]][["msd"]][-c(1,2,3)]

cred_lower = fit_C_Model@.MISC[["summary"]][["c_quan"]][,1,1][-c(1,2,3,144)]
cred_upper = fit_C_Model@.MISC[["summary"]][["c_quan"]][,9,1][-c(1,2,3,144)]

N = nrow(c_data)
xGrid = seq(1,140,1)
plot(xGrid, c_data$c, col = "red", main = "Possion model & Data", ylab = "nr of Cases", xlab = "Draw")
lines(exp(poster_means_vector[1:N]), ylim = c(0,60), col = "black", type = 'l')
lines(exp(cred_lower), col = "blue")
lines(exp(cred_upper), col = "green")

legend('topleft',legend = c('Data', 'Mean','Lower','Upper'), col = c('red', 'black', 'blue', 'green'), lwd=2)





