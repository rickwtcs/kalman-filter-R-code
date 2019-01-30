# Lesson demo 3
# how to predict future system states and observations with the EKF.
# Credit to Stefan Gelissen
# http://blogs2.datall-analyse.nl/2016/02/16/rcode_extended_kalman_filter_forecasting/

library(dlm)
library(numDeriv)

#The R functions for the extended Kalman filter (EKF) and
#forecasting with the EKF can be downloaded from: 
source("http://www.datall-analyse.nl/EKF.R")
#Take a look at the functions dlmExtFilter (=the EKF) and dlmExtForecast (=the
#forecast function), and notice that at the beginning of the scripts you
#will find a description of the functions' arguments.
dlmExtFilter
dlmExtForecast


##EXAMPLE

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

#measurement of the population sizes at the observation times
Pm <- P + rnorm(length(t), 0, sqrt(varPop))

#measured population sizes
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(t, Pm, type="l", col=gray(level=.7))
legend("topleft", lty=c(2, 1),
       col=c("red", col=gray(level=.7)), legend=c("true state", "measured"),
       bty="n", y.intersp=1.2, cex=.7)

#store measurement data
dataEx1 <- Pm



##FORECASTING WITH THE EXTENDED KALMAN FILTER (EKF)

#Specify dynamic model having 2 states, namely [growth rate, population size].
ex1 <- list(
  #Initial state estimates growth rate and population size.
  #Notice that the specified initial state estimate (at t=0) 
  #for population size deviates from the value that was used for
  #generating the data (i.e., 10).
  #You may change change these initial state estimates too and see how
  #they influence the behavior of the extended Kalman filter.
  m0=c(.2, 5),
  #Error covariances of the initial state estimates:
  #this matrix reflects the uncertainty in our initial state estimates.
  #You may change the values in this matrix and see how they influence
  #the behavior of the Kalman filter.
  C0=diag(c(100, 100)),
  #Measurement noise
  V=varPop,
  #Process noise
  W=diag(rep(0,2)))
#Note that all covariances in the process noise matrix (W) are set to zero.
#This makes sense since the change in the growth rate and population size
#at each time step is fully explained by our logistic population model.

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

##Compute the filtered (a posteriori) state estimates
ext1 <- dlmExtFilter(y=dataEx1, mod=ex1,
                     GGfunction=GGfunction, FFfunction=FFfunction)

#plot the filtered state estimates 
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(c(0,t), ext1$m[,2], lty=2, col="blue", cex=.7)
legend("topleft", lty=c(2, 2), 
       col=c("red", "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)
#Notice that the EKF's behavior at earlier observation times
#(at approx. Time<10) is less stable than the behavior at later
#observation times (at approx. Time>10).
#Below we will see if and how this stability might influence
#the EKF's forecasts.


##Compute the filtered (a posteriori) state estimates
#EKF for the first 50 observations
#Note that these 50 observations comprise "earlier observation times"
#as stated above.
ext2 <- dlmExtFilter(y=dataEx1[1:50], mod=ex1,
                     GGfunction=GGfunction, FFfunction=FFfunction)

#Compute forecasts for 50 time steps ahead
fc <- dlmExtForecast(ext2, nAhead=50)


#95% confidence intervals for forecasts

#Future system states
seFoA <- sqrt(unlist(lapply(fc$R, function(x) x[2,2])))
ciFoA <- fc$a[,2] + qnorm(.05/2)*seFoA%o%c(1, -1)

#Plot future system states including confidence intervals
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(t[1:50], ext2$m[-1,2], col=gray(level=.5), lwd=1)
lines(t[51:100], fc$a[,2], col="blue", lwd=2)
lines(t[51:100], ciFoA[,1], col="blue", lty=2, lwd=1)
lines(t[51:100], ciFoA[,2], col="blue", lty=2, lwd=1)
legend("topleft", lty=c(2, 1), 
       col=c("red", "blue"),
       legend=c("true state", "future state"),
       bty="n", y.intersp=1.2, cex=.7)
#Notice that the confidence intervals get wider for
#forecasts further ahead in time.


#Future observations
seFoX <- sqrt(unlist(lapply(fc$Q, function(x) x[1,1])))
ciFoX <- fc$f[,1] + qnorm(.05/2)*seFoX%o%c(1, -1)

#Plot future observations including prediction intervals
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(t[1:50], ext2$f[,1], col=gray(level=.1), lwd=1)
lines(t, Pm, type="l", col=gray(level=.7))
lines(t[51:100], fc$f[,1], col="blue", lwd=2)
lines(t[51:100], ciFoX[,1], col="blue", lty=2, lwd=1)
lines(t[51:100], ciFoX[,2], col="blue", lty=2, lwd=1)
legend("topleft", lty=c(2, 1, 1), 
       col=c("red", col=gray(level=.7), "blue"),
       legend=c("true state", "measurements", "future observations"),
       bty="n", y.intersp=1.2, cex=.7)
#Notice again that the prediction intervals get wider for
#forecasts further ahead in time.


#EKF for the first 150 observations
#Note that these 150 observations comprise "later observation times"
#as stated above.
ext3 <- dlmExtFilter(y=dataEx1[1:150], mod=ex1,
                     GGfunction=GGfunction, FFfunction=FFfunction)

#Compute forecasts for 50 time steps ahead
fc <- dlmExtForecast(ext3, nAhead=50)


#95% confidence intervals for forecasts

#Future system states
seFoA <- sqrt(unlist(lapply(fc$R, function(x) x[2,2])))
ciFoA <- fc$a[,2] + qnorm(.05/2)*seFoA%o%c(1, -1)

#Plot future system states including confidence intervals
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(t[1:150], ext3$m[-1,2], col=gray(level=.5), lwd=1)
lines(t[151:200], fc$a[,2], col="blue", lwd=2)
lines(t[151:200], ciFoA[,1], col="blue", lty=2, lwd=1)
lines(t[151:200], ciFoA[,2], col="blue", lty=2, lwd=1)
legend("topleft", lty=c(2, 1), 
       col=c("red", "blue"),
       legend=c("true state", "future state"),
       bty="n", y.intersp=1.2, cex=.7)
#Notice that all the confidence intervals have approximately the
#same width (also for forecasts further ahead in time).


#Future observations
seFoX <- sqrt(unlist(lapply(fc$Q, function(x) x[1,1])))
ciFoX <- fc$f[,1] + qnorm(.05/2)*seFoX%o%c(1, -1)

#Plot future observations including prediction intervals
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(t[1:150], ext3$f[,1], col=gray(level=.1), lwd=1)
lines(t, Pm, type="l", col=gray(level=.7))
lines(t[151:200], fc$f[,1], col="blue", lwd=2)
lines(t[151:200], ciFoX[,1], col="blue", lty=2, lwd=1)
lines(t[151:200], ciFoX[,2], col="blue", lty=2, lwd=1)
legend("topleft", lty=c(2, 1, 1), 
       col=c("red", col=gray(level=.7), "blue"),
       legend=c("true state", "measurements", "future observations"),
       bty="n", y.intersp=1.2, cex=.7)
#Notice again that all the prediction intervals have approximately
#the same width (also for forecasts further ahead in time).

#In the above example we saw that our observed stability of the
#EKF influenced the width of the forecasts' confidence/prediction
#intervals.
#That is, for the "stable regime" (for Time>10), the widths of the 
#confidence/prediction intervals are similar over time.
#In contrast, for the "unstable regime" (for Time<10), the
#forecasts' confidence/prediction intervals get wider over time.



#In the following example we will make the EKF more stable at
#earlier observation times. As a consequence, we now expect that
#the width of the confidence/prediction intervals for forecasts at
#"earlier observation times" will also not vary substantially over time.

#Specify dynamic model having 2 states, namely [growth rate, population size].
ex2 <- list(
  #Initial state estimates growth rate and population size.
  m0=c(.2, 10),
  #Error covariances of the initial state estimates:
  #setting these error covariances to 0 will make the EKF
  #more stable at "earlier observation times".
  C0=diag(c(0, 0)),
  #Measurement noise
  V=varPop,
  #Process noise
  W=diag(rep(0,2)))

##Compute the filtered (a posteriori) state estimates
ext4 <- dlmExtFilter(y=dataEx1, mod=ex2,
                     GGfunction=GGfunction, FFfunction=FFfunction)

#plot the filtered state estimates 
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(c(0,t), ext4$m[,2], lty=2, col="blue", cex=.7)
legend("topleft", lty=c(2, 2), 
       col=c("red", "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)
#Notice that the EKF's behavior at "earlier observation times"
#is now much more stable.


##Compute the filtered (a posteriori) state estimates
#EKF for the first 50 observations
#Note that these 50 observations comprise "earlier observation times"
ext5 <- dlmExtFilter(y=dataEx1[1:50], mod=ex2,
                     GGfunction=GGfunction, FFfunction=FFfunction)

#Compute forecasts for 50 time steps ahead
fc <- dlmExtForecast(ext5, nAhead=50)


#95% confidence intervals for forecasts

#Future system states
seFoA <- sqrt(unlist(lapply(fc$R, function(x) x[2,2])))
ciFoA <- fc$a[,2] + qnorm(.05/2)*seFoA%o%c(1, -1)

#Plot future system states including confidence intervals
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(t[1:50], ext5$m[-1,2], col=gray(level=.5), lwd=1)
lines(t[51:100], fc$a[,2], col="blue", lwd=2)
lines(t[51:100], ciFoA[,1], col="blue", lty=2, lwd=1)
lines(t[51:100], ciFoA[,2], col="blue", lty=2, lwd=1)
legend("topleft", lty=c(2, 1), 
       col=c("red", "blue"),
       legend=c("true state", "future state"),
       bty="n", y.intersp=1.2, cex=.7)
#Notice that the width of the confidence intervals does
#not change that much over time.


#Future observations
seFoX <- sqrt(unlist(lapply(fc$Q, function(x) x[1,1])))
ciFoX <- fc$f[,1] + qnorm(.05/2)*seFoX%o%c(1, -1)

#Plot future observations including prediction intervals
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(t[1:50], ext5$f[,1], col=gray(level=.1), lwd=1)
lines(t, Pm, type="l", col=gray(level=.7))
lines(t[51:100], fc$f[,1], col="blue", lwd=2)
lines(t[51:100], ciFoX[,1], col="blue", lty=2, lwd=1)
lines(t[51:100], ciFoX[,2], col="blue", lty=2, lwd=1)
legend("topleft", lty=c(2, 1, 1), 
       col=c("red", col=gray(level=.7), "blue"),
       legend=c("true state", "measurements", "future observations"),
       bty="n", y.intersp=1.2, cex=.7)
#Notice again that the prediction intervals do not get
#substantially wider for forecasts further ahead in time.


#EKF for the first 150 observations
#Note that these 150 observations comprise "later observation times"
ext6 <- dlmExtFilter(y=dataEx1[1:150], mod=ex2,
                     GGfunction=GGfunction, FFfunction=FFfunction)

#Compute forecasts for 50 time steps ahead
fc <- dlmExtForecast(ext6, nAhead=50)


#95% confidence intervals for forecasts

#Future system states
seFoA <- sqrt(unlist(lapply(fc$R, function(x) x[2,2])))
ciFoA <- fc$a[,2] + qnorm(.05/2)*seFoA%o%c(1, -1)

#Plot future system states including confidence intervals
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(t[1:150], ext6$m[-1,2], col=gray(level=.5), lwd=1)
lines(t[151:200], fc$a[,2], col="blue", lwd=2)
lines(t[151:200], ciFoA[,1], col="blue", lty=2, lwd=1)
lines(t[151:200], ciFoA[,2], col="blue", lty=2, lwd=1)
legend("topleft", lty=c(2, 1), 
       col=c("red", "blue"),
       legend=c("true state", "future state"),
       bty="n", y.intersp=1.2, cex=.7)
#Notice that all the confidence intervals have approximately the
#same width (also for forecasts further ahead in time).


#Future observations
seFoX <- sqrt(unlist(lapply(fc$Q, function(x) x[1,1])))
ciFoX <- fc$f[,1] + qnorm(.05/2)*seFoX%o%c(1, -1)

#Plot future observations including prediction intervals
plot(t, P, type="l", lty=2, col="red", ylim=range(Pm),
     xlab="Time", ylab="Population size")
lines(t[1:150], ext6$f[,1], col=gray(level=.1), lwd=1)
lines(t, Pm, type="l", col=gray(level=.7))
lines(t[151:200], fc$f[,1], col="blue", lwd=2)
lines(t[151:200], ciFoX[,1], col="blue", lty=2, lwd=1)
lines(t[151:200], ciFoX[,2], col="blue", lty=2, lwd=1)
legend("topleft", lty=c(2, 1, 1), 
       col=c("red", col=gray(level=.7), "blue"),
       legend=c("true state", "measurements", "future observations"),
       bty="n", y.intersp=1.2, cex=.7)
#Notice again that all the prediction intervals have approximately
#the same width (also for forecasts further ahead in time).
