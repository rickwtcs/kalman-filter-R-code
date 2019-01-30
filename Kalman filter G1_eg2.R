# Source http://blogs2.datall-analyse.nl/2016/02/11/rcode_extended_kalman_filter/
# R Code, simulations and modeling
# by Stefan Gelissen

library(dlm)
library(numDeriv)

#The R functions for the extended Kalman filter (EKF) and its smoother
#can be downloaded from: 
source("http://www.datall-analyse.nl/EKF.R")
#Take a look at the functions dlmExtFilter (=the EKF) and dlmExtSmooth (=the
#RTS smoother), and notice that at the beginning of the scripts you will find
#a description of the functions' arguments.
dlmExtFilter
dlmExtSmooth
#You may save these functions to your computer using the save function in R.


##EXAMPLE 2

#Similar to Example 1, but in this example we assume that we only
#measured the number of prey at each time step.

#Data for Example 2
dataEx2 <- yprey

#Dynamic model:
#specifying 2 states, namely [number of prey, number of predators].
#Note that the specified initial state estimates (at t=0) for prey and predators
#deviate from the values that were used for generating the data (i.e., 400
#for prey and 200 for predators).
ex2 <- list(m0=c(390, 220), #initial state estimates
            #error covariances of the initial state estimates:
            #this matrix reflects the uncertainty in our initial state estimates
            #you may change the values in this matrix and see how they influence
            #the behavior of the Kalman filter
            C0=diag(rep(100,2)),
            #measurement noise
            V=sigmaPrey^2,
            #process noise
            W=diag(rep(0,2)))

#Specify the state transition function
#REMEMBER: always use arguments x and k when specifying the GGfunction
GGfunction <- function (x, k){
  prey <- x[1]; pred <- x[2]
  c(prey + (alpha*prey - beta*prey*pred)*dT,
    pred + (delta*prey*pred - gamma*pred)*dT)}

#Specify the observation function
#REMEMBER: always use arguments x and k when specifying the FFfunction
FFfunction <- function (x, k){
  prey <- x[1]
  c(prey)}



##Compute the filtered (a posteriori) state estimates
ext2 <- dlmExtFilter(y=dataEx2, mod=ex2,
                     GGfunction=GGfunction, FFfunction=FFfunction)

#plot the filtered state estimates
plot(x[,1], x[,2], type="l", lwd=2, xlab="Prey", ylab="Predators",
     xlim=range(yprey), ylim=range(ypred), col=gray(level=.5))
lines(ext2$m[,1], ext2$m[,2], lty=2, col="blue", lwd=2)
legend("topright", lty=c(1, 2), lwd=c(2, 2),
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)


#compute the confidence intervals for the filtered state estimates
varcovFilteredState <- dlmSvd2var(ext2$U.C, ext2$D.C)
#95% ci for prey
seFY <- sqrt(unlist(lapply(varcovFilteredState, function(x) x[1,1])))
ciFY <- ext2$m[,1] + qnorm(.05/2)*seFY%o%c(1, -1)
#95% ci for predators
seFD <- sqrt(unlist(lapply(varcovFilteredState, function(x) x[2,2])))
ciFD <- ext2$m[,2] + qnorm(.05/2)*seFD%o%c(1, -1)

#prey population
plot(st, x[,1], type="l", xlab="time", ylab="Prey",
     col=gray(level=.7), lwd=1.5)
lines(c(0,mt), ciFY[,1], lty=2, col="blue")
lines(c(0,mt), ciFY[,2], lty=2, col="blue")
legend("bottomright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)

#predator population
plot(st, x[,2], type="l", xlab="time", ylab="Prey",
     col=gray(level=.7), lwd=1.5)
lines(c(0,mt), ciFD[,1], lty=2, col="blue")
lines(c(0,mt), ciFD[,2], lty=2, col="blue")
legend("topright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)



##Compute smoothed state estimates

#apply fixed-interval smoothing
smoothState <- dlmExtSmooth(ext2)
plot(x[,1], x[,2], type="l", lwd=2, xlab="Prey", ylab="Predators",
     xlim=range(yprey), ylim=range(ypred), col=gray(level=.5))
lines(ext2$m[,1], ext2$m[,2], type="l", lty=2, col="blue")
lines(smoothState$s[,1], smoothState$s[,2], type="l", lty=2, col="darkgreen")
legend("topright", lty=c(1, 2), lwd=c(1, 1),
       col=c("black", "blue", "darkgreen"),
       legend=c("true state", "filtered state", "smoothed state"),
       bty="n", y.intersp=1.2, cex=.7)

#confidence intervals for smoothed state estimates
varcovSmoothState <- dlmSvd2var(smoothState$U.S, smoothState$D.S)
#95% ci for prey population
seSY <- sqrt(unlist(lapply(varcovSmoothState, function(x) x[1,1])))
ciSY <- smoothState$s[,1] + qnorm(.05/2)*seSY%o%c(1, -1)
#95% ci for predator population
seSD <- sqrt(unlist(lapply(varcovSmoothState, function(x) x[2,2])))
ciSD <- smoothState$s[,2] + qnorm(.05/2)*seSY%o%c(1, -1)

#plot prey
plot(st, x[,1], type="l", xlab="time", ylab="Prey",
     col=gray(level=.7), lwd=1.5)
lines(c(0,mt), ciSY[,1], lty=2, col="blue")
lines(c(0,mt), ciSY[,2], lty=2, col="blue")
legend("bottomright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "smoothed state"),
       bty="n", y.intersp=1.2, cex=.7)

#plot predators
plot(st, x[,2], type="l", xlab="time", ylab="Predator",
     col=gray(level=.7), lwd=1.5)
lines(c(0,mt), ciSD[,1], lty=2, col="blue")
lines(c(0,mt), ciSD[,2], lty=2, col="blue")
legend("topright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "smoothed state"),
       bty="n", y.intersp=1.2, cex=.7)

# End of example 2