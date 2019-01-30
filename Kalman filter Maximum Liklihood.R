# Lesson demo 2
# EKF = Extended Kalman filter
# how to fit (to estmiate) unknown parameters of an EKF model by means of likelihood maximization.
# Credit to Stefan Gelissen
# http://blogs2.datall-analyse.nl/2016/02/15/rcode_extended_kalman_filter_maximum_likelihood/

library(dlm)
library(numDeriv)

#The R functions for implementing the extended Kalman filter (EKF) and its
#smoother can be downloaded from: 
source("http://www.datall-analyse.nl/EKF.R")
#Take a look at the functions dlmExtFilter (=the EKF) and dlmExtSmooth (=the
#RTS smoother), and notice that at the beginning of the scripts you will find
#a description of the functions' arguments.
dlmExtFilter
dlmExtSmooth


##EXAMPLE 1

#In this example we will apply the Extended Kalman Filter (EKF) to the
#logistic population model. The logistic population model describes the
#growth of a population over a time period.
#See the following Wikipedia article for more information on this model:
#en.wikipedia.org/wiki/Logistic_function#In_ecology:_modeling_population_growth


##TRUE STATES

#parameters
K <- 100 #carrying capacity
P0 <- 10 #population size at t=0
r <- .2 #growth rate

#sample time
dT <- .1
#you may change the value of dT and see how it influences
#the behavior of the extended Kalman filter

#observation times
t <- seq(0.1, 25, dT)

#simulate true state for population size (=P) at the observation times
P <- K*P0*exp(r*t) / (K + P0*(exp(r*t) - 1))

#plot simulated population size
plot(t, P, type="l", xlab="Time", ylab="Population size")



##MEASUREMENTS

#assume that the measurement noise is normally distributed
#with mean=0 and variance=25
varPop <- 25

#measurement of the population sizes at the specified observation times
Pm <- P + rnorm(length(t), 0, sqrt(varPop))

#plot the measured population sizes
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(t, Pm, type="l", col=gray(level=.7))
legend("topleft", lty=c(2, 1),
       col=c("red", col=gray(level=.7)), legend=c("true state", "measured"),
       bty="n", y.intersp=1.2, cex=.7)

#store measurement data
dataEx1 <- Pm



##ESTIMATING THE PARAMETERS OF AN EXTENDED KALMAN FILTER (EKF) MODEL

#Specify the state transition function for the EKF
#WARNING: always use arguments x and k when specifying the GGfunction
GGfunction <- function(x, k){
  r <- x[1]; P <- x[2]
  c(r,
    K * P * exp(r*dT) / (K + P * (exp(r*dT) - 1)))}

#Specify the observation/measurement function for the EKF
#WARNING: always use arguments x and k when specifying the FFfunction
FFfunction <- function (x, k){
  r <- x[1]; P <- x[2]
  c(P)}


#Function for computing the negative log-likelihood of an EKF model
llikss <- function (x, data) {
  #Specify dynamic model having 2 states, namely [growth rate, population size].
  mod <- list(
    #Initial state estimates growth rate and population size.
    #Notice that the specified initial state estimate (at t=0) 
    #for population size deviates from the value that was used for
    #generating the data (i.e., 10).
    #You may change these initial state estimates too and see how
    #they influence the behavior of the extended Kalman filter.
    m0=c(0.2, 5),
    #Error covariances of the initial state estimates:
    #this matrix reflects the uncertainty in our initial state estimates.
    #You may change the values in this matrix and see how they influence
    #the behavior of the Kalman filter.
    C0=diag(c(100, 100)),
    #Measurement noise (to be estimated)
    V=exp(x[1]),
    #Process noise
    W=diag(rep(0,2)))
  #Note that all covariances in the process noise matrix (W) are set to zero.
  #This makes sense since the change in the growth rate and population size
  #at each time step is fully explained by our logistic population model.
  
  #compute and return the negative log-likelihood
  dlmExtFilter(y=data, mod=mod,
               FFfunction=FFfunction, GGfunction=GGfunction,
               simplify=TRUE, logLik=TRUE)$logLik
}

#Minimize the negative log-likelihood:
#instead of maximizing the log-likelihood, we will minimize the negative
#log-likelihood.
#Use some arbitrary starting value for the to be estimated measurement noise
#(here we will use log(10), but it is possible to use some other
#arbitrary starting value).
mod <- optim(log(10), llikss,
             lower=log(1e-8),
             method="L-BFGS-B", hessian=TRUE,
             data=dataEx1)
#In case the BFGS method generates an error, then try using the
#Brent method for this one-dimensional optimization problem.

#Check for convergence
mod$message
mod$convergence #successful optimization yields value 0

#Estimated measurement noise
exp(mod$par)

#Construct 95% confidence interval for the estimated measurement variance
seParms <- sqrt(diag(solve(mod$hessian)))
exp(mod$par + qnorm(.05/2)*seParms%o%c(1,-1))

#True value for the measurement noise
varPop



##Compute the filtered (a posteriori) state estimates

#Using the estimated measurement variance
estVarPop <- exp(mod$par)

#Specify dynamic model having 2 states, namely [growth rate, population size].
ex1 <- list(m0=c(.2, 5), #initial state estimates
            #error covariances of the initial state estimates:
            C0=diag(c(100, 100)),
            #estimated measurement noise
            V=estVarPop,
            #process noise
            W=diag(rep(0,2)))

#Filtered (a posteriori) state estimates
ext1 <- dlmExtFilter(y=dataEx1, mod=ex1,
                     GGfunction=GGfunction, FFfunction=FFfunction)

#plot the filtered state estimates 
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(t, Pm, type="l", col=gray(level=.7))
lines(c(0,t), ext1$m[,2], lty=2, col="blue", cex=.7)
legend("topleft", lty=c(2, 1, 2), 
       col=c("red", gray(level=.7), "blue"),
       legend=c("true state", "measured", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)

#compute the confidence intervals for the filtered state estimates
varcovFilteredState <- dlmSvd2var(ext1$U.C, ext1$D.C)
#95% ci for growth rate
seFR <- sqrt(unlist(lapply(varcovFilteredState, function(x) x[1,1])))
ciFR <- ext1$m[,1] + qnorm(.05/2)*seFR%o%c(1, -1)
#95% ci for population size
seFP <- sqrt(unlist(lapply(varcovFilteredState, function(x) x[2,2])))
ciFP <- ext1$m[,2] + qnorm(.05/2)*seFP%o%c(1, -1)

#plot growth rate
plot(t, rep(r, length(t)), type="l", xlab="Time", ylab="Growth rate",
     ylim=c(min(ciFR[,1]), max(ciFR[,2])), col=gray(level=.7), lwd=1.5)
lines(c(0,t), ciFR[,1], lty=2, col="blue")
lines(c(0,t), ciFR[,2], lty=2, col="blue")
legend("topright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "smoothed state"),
       bty="n", y.intersp=1.2, cex=.7)

#plot population size
plot(t, P, type="l", xlab="Time", ylab="Population size",
     ylim=c(min(ciFP[,1]), max(ciFP[,2])), col=gray(level=.7), lwd=1.5)
lines(c(0,t), ciFP[,1], lty=2, col="blue")
lines(c(0,t), ciFP[,2], lty=2, col="blue")
legend("topleft", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "smoothed state"),
       bty="n", y.intersp=1.2, cex=.7)


##Compute smoothed state estimates

#apply fixed-interval smoothing
smoothState <- dlmExtSmooth(ext1)

#plot the smoothed state estimates for population size 
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(t, Pm, type="l", col=gray(level=.7))
lines(c(0,t), smoothState$s[,2], lty=2, col="blue", cex=.7)
legend("topleft", lty=c(2, 1, 2), 
       col=c("red", gray(level=.7), "blue"),
       legend=c("true state", "measured", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)

#confidence intervals for smoothed state estimates
varcovSmoothState <- dlmSvd2var(smoothState$U.S, smoothState$D.S)
#95% ci for growth rate
seSR <- sqrt(unlist(lapply(varcovSmoothState, function(x) x[1,1])))
ciSR <- smoothState$s[,1] + qnorm(.05/2)*seSR%o%c(1, -1)
#95% ci for population sizt
seSP <- sqrt(unlist(lapply(varcovSmoothState, function(x) x[2,2])))
ciSP <- smoothState$s[,2] + qnorm(.05/2)*seSP%o%c(1, -1)

#plot growth rate
plot(t, rep(r, length(t)), type="l", xlab="Time", ylab="Growth rate",
     ylim=c(.05, .5), col=gray(level=.7), lwd=1.5)
lines(c(0,t), ciSR[,1], lty=2, col="blue")
lines(c(0,t), ciSR[,2], lty=2, col="blue")
legend("topleft", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "smoothed state"),
       bty="n", y.intersp=1.2, cex=.7)

#plot population size
plot(t, P, type="l", xlab="Time", ylab="Population size",
     ylim=c(0,100), col=gray(level=.7), lwd=1.5)
lines(c(0,t), ciSP[,1], lty=2, col="blue")
lines(c(0,t), ciSP[,2], lty=2, col="blue")
legend("topleft", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "smoothed state"),
       bty="n", y.intersp=1.2, cex=.7)





##EXAMPLE 2

#In this example we will apply the Extended Kalman Filter (EKF) to the
#Lotka-Volterra system of nonlinear differential equations. The Lotka-Volterra
#equations describe how prey and predators interact.
#See the following Wikipedia article for more information on these equations:
#en.wikipedia.org/wiki/Lotka-Volterra_equations



##TRUE STATES

#Generate the amount of prey and predators based on the Lotka-Volterra equations.
#Data are generated using the forward Euler method.
dt <- 1/5000 #time step for Euler integration
tt <- 7 #upper bound of time window
st <- seq(0, tt, by=dt) #lower time bounds of the integration intervals
ns <- length(st) #number of Euler integrations
x <- matrix(0, ncol=2, nrow=ns) #specify matrix for states
x[1,] <- c(400, 200) #respective amounts of prey and predators at t=0
colnames(x) <- c("Prey", "Predators")

#parameters in the Lotka-Volterra model
alpha <- 1
beta <- 1/300
delta <- 1/200
gamma <- 1

#simulate true states
for (i in 2:ns) {
  #prey population
  x[i,1] <- x[i-1,1] + (alpha*x[i-1,1] - beta*x[i-1,1]*x[i-1,2])*dt
  #predator population
  x[i,2] <- x[i-1,2] + (delta*x[i-1,1]*x[i-1,2] - gamma*x[i-1,2])*dt
}

#phase-space plot of simulated predator-prey system
plot(x[,1], x[,2], type="l", xlab="Prey", ylab="Predators")



##MEASUREMENTS

#Take measurements with a sample time of .01
dT <- .01 #sample time for measurements
#you may change the value of dT and see how it influences
#the behavior of the extended Kalman filter
nm <- tt/dT #total number of measurements
mt <- seq(dT, tt, dT) #measurement times

#standard deviations for the measurement noise
sigmaPrey <- 7 #prey
sigmaPred <- 10 #predators

#measurements at specified measurement times
yprey <- sapply(1:nm, function(i) x[i*((ns-1)/nm) + 1, 1] + rnorm(1, 0, sigmaPrey))
ypred <- sapply(1:nm, function(i) x[i*((ns-1)/nm) + 1, 2] + rnorm(1, 0, sigmaPred))

#store measurement data
dataEx2 <- cbind(yprey, ypred)

#plot the generated measurements
plot(x[,1], x[,2], type="l", xlab="Prey", ylab="Predators",
     xlim=range(yprey), ylim=range(ypred))
points(yprey, ypred, cex=.5, col="red")

#plot of time against measurements
par(mfrow=c(2, 1))
plot(st, x[,1], type="l", xlab="Time", ylab="Prey",
     ylim=range(yprey))
lines(mt, yprey, col="darkgreen")

plot(st, x[,2], type="l", xlab="Time", ylab="Predators",
     ylim=range(ypred))
lines(mt, ypred, col="brown")
par(mfrow=c(1, 1))



##ESTIMATING THE PARAMETERS OF AN EXTENDED KALMAN FILTER (EKF) MODEL

#Specify the state transition function for the EKF:
#note that we will use as state functions the difference equations
#given by Euler's forward method. These difference equations will yield valid 
#estimates for the amounts of prey and predators at each time step as long
#as the specified value for dT above is relatively small.
#WARNING: always use arguments x and k when specifying the GGfunction
GGfunction <- function (x, k){
  prey <- x[1]; pred <- x[2]
  c(prey + (alpha*prey - beta*prey*pred)*dT,
    pred + (delta*prey*pred - gamma*pred)*dT)}

#Specify the observation/measurement function for the EKF
#WARNING: always use arguments x and k when specifying the FFfunction
FFfunction <- function (x, k){
  prey <- x[1]; pred <- x[2]
  c(prey, pred)}


#Function for computing the negative log-likelihood
llikss <- function (x, data) {
  #Specify dynamic model having 2 states, namely
  #[amount of prey, amount of predators].
  mod <- list(
    #Initial state estimates for amount of prey and predators.
    #Notice that the specified initial state estimates (at t=0) below for prey
    #and predators deviate from the values that were used for generating the
    #data (i.e., 400 for prey and 200 for predators).
    #You may change these initial state estimates too and see how
    #they influence the behavior of the extended Kalman filter.
    m0=c(390, 220),
    #Error covariances of the initial state estimates:
    #this matrix reflects the uncertainty in our initial state estimates.
    #You may change the values in this matrix and see how they influence
    #the behavior of the Kalman filter.
    C0=diag(100, 2),
    #Measurement noise (to be estimated)
    V=matrix(c(exp(x[1]), 0,
               0, exp(x[2])), nrow=2, byrow=TRUE),
    #Process noise
    W=diag(rep(0,2)))
  #Note that all covariances in the process noise matrix (W) are set to zero.
  #This makes sense since the change in the amount of prey and predators at
  #each time step is fully explained by our Lotka-Volterra model describing the
  #prey-predator interaction.
  
  #Before computing the negative log-likelihood for the current parameter
  #estimates, first check if the resulting covariance matrix
  #for the measurement noise (=V) is positive semidefinite.
  #If the resulting covariance matrix appears not to be positive
  #semidefinite, then assign NaN to the negative log-likelihood value.
  if (all(eigen(mod$V)$values >= 0))
    dlmExtFilter(data, mod,
                 FFfunction=FFfunction, GGfunction=GGfunction,
                 simplify=TRUE, logLik=TRUE)$logLik
  else NaN
}

#Minimize the negative log-likelihood.
#Use some arbitrary starting values for the to be estimated measurement
#covariances (here we will use log(50), but it is possible
#to use some other arbitrary starting value).
mod <- optim(rep(log(50), 2), llikss,
             lower=rep(log(1), 2),
             method="L-BFGS-B", hessian=TRUE,
             data=dataEx2)

#Check for convergence
mod$message
mod$convergence #successful optimization yields value 0

#Estimated measurement covariances
exp(mod$par)

#Construct 95% confidence interval for the estimated measurement covariances
seParms <- sqrt(diag(solve(mod$hessian)))
exp(mod$par + qnorm(.05/2)*seParms%o%c(1,-1))

#True values for the measurement covariances
c(sigmaPrey^2, sigmaPred^2)

#Check by hand if the resulting V-matrix is indeed positive semidefinite
v <- matrix(c(exp(mod$par)[1], 0,
              0, exp(mod$par)[2]), nrow=2, byrow=TRUE)
all(eigen(v)$values >= 0)



#In case the BFGS method generates an error, then try using the
#more robust Nelder-Mead method
mod <- optim(rep(log(50), 2), llikss,
             method="Nelder-Mead", hessian=TRUE,
             data=dataEx2)
#In case this Nelder-Mead generates an error, then try using
#the same Nelder-Mead but with hessian=FALSE

#Check for convergence
mod$convergence #successful optimization yields value 0

#Estimated measurement covariances
exp(mod$par)

#Construct 95% confidence interval for the estimated measurement covariances
seParms <- sqrt(diag(solve(mod$hessian)))
exp(mod$par + qnorm(.05/2)*seParms%o%c(1,-1))