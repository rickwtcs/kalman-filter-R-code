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





#EXAMPLE 3

#In the above examples we assumed that the values for the parameters
#beta and gamma were somehow known to us. Let us now include these
#two parameters as states in our model and have the EKF estimate their values.

#Data
dataEx3 <- cbind(yprey, ypred)

#Dynamic model:
#specifying 4 states, namely [number of prey, number of predators, beta, delta].
#Note that the specified initial state estimates (at t=0) deviate from the
#values that were used for generating the data.
#Remember that the values for beta and delta used in generating the data
#were 1/300 and 1/200, respectively.
ex3 <- list(m0=c(390, 220, 1/400, 1/100), #initial state estimates
            #error covariances of the initial state estimates:
            #this matrix reflects the uncertainty in our initial state estimates
            #you may change the values in this matrix and see how they influence
            #the behavior of the Kalman filter
            C0=diag(rep(1000,4)),
            #measurement noise
            V=diag(c(sigmaPrey^2, sigmaPred^2)),
            #process noise
            W=diag(rep(0,4)))

#Specify the state transition function
#REMEMBER: always use arguments x and k when specifying the GGfunction
GGfunction <- function (x, k){
  prey <- x[1]; pred <- x[2]; beta <- x[3]; delta <- x[4]
  c(prey + (alpha*prey - beta*prey*pred)*dT,
    pred + (delta*prey*pred - gamma*pred)*dT,
    beta,
    delta)}

#Specify the observation function
#REMEMBER: always use arguments x and k when specifying the FFfunction
FFfunction <- function (x, k){
  prey <- x[1]; pred <- x[2]
  c(prey, pred)}



##Compute the filtered (a posteriori) state estimates
ext3 <- dlmExtFilter(y=dataEx3, mod=ex3,
                     GGfunction=GGfunction, FFfunction=FFfunction)

#plot the filtered state estimates
plot(x[,1], x[,2], type="l", lwd=2, xlab="Prey", ylab="Predators",
     xlim=range(yprey), ylim=range(ypred), col=gray(level=.5))
lines(ext3$m[,1], ext3$m[,2], lty=2, col="blue", lwd=2)
legend("topright", lty=c(1, 2), lwd=c(2, 2),
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)

#compute the confidence intervals for the filtered state estimates
varcovFilteredState <- dlmSvd2var(ext3$U.C, ext3$D.C)
#95% ci for prey
seFY <- sqrt(unlist(lapply(varcovFilteredState, function(x) x[1,1])))
ciFY <- ext3$m[,1] + qnorm(.05/2)*seFY%o%c(1, -1)
#95% ci for predators
seFD <- sqrt(unlist(lapply(varcovFilteredState, function(x) x[2,2])))
ciFD <- ext3$m[,2] + qnorm(.05/2)*seFD%o%c(1, -1)
#95% ci for beta
seFB <- sqrt(unlist(lapply(varcovFilteredState, function(x) x[3,3])))
ciFB <- ext3$m[,3] + qnorm(.05/2)*seFB%o%c(1, -1)
#95% ci for delta
seFE <- sqrt(unlist(lapply(varcovFilteredState, function(x) x[4,4])))
ciFE <- ext3$m[,4] + qnorm(.05/2)*seFE%o%c(1, -1)

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

#beta
plot(x=NA, y=NA, type="n", xlab="time", ylab="Beta",
     xlim=c(0,7), ylim=c(0, .005))
abline(h=beta, lty=2, col=gray(level=.5))
lines(c(0,mt), ciFB[,1], lty=2, col="blue")
lines(c(0,mt), ciFB[,2], lty=2, col="blue")
legend("bottomright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)

#delta
plot(x=NA, y=NA, type="n", xlab="time", ylab="Gamma",
     xlim=c(0,7), ylim=c(.004, .006))
abline(h=delta, lty=2, col=gray(level=.5))
lines(c(0,mt), ciFE[,1], lty=2, col="blue")
lines(c(0,mt), ciFE[,2], lty=2, col="blue")
legend("bottomright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)



##Compute smoothed state estimates

#apply fixed-interval smoothing
smoothState <- dlmExtSmooth(ext3)

#confidence intervals for smoothed state estimates
varcovSmoothState <- dlmSvd2var(smoothState$U.S, smoothState$D.S)
#95% ci for beta
seSD <- sqrt(unlist(lapply(varcovSmoothState, function(x) x[3,3])))
ciSD <- smoothState$s[,3] + qnorm(.05/2)*seSD%o%c(1, -1)
#95% ci for delta
seSE <- sqrt(unlist(lapply(varcovSmoothState, function(x) x[4,4])))
ciSE <- smoothState$s[,4] + qnorm(.05/2)*seSE%o%c(1, -1)

#beta
plot(x=NA, y=NA, type="n", xlab="time", ylab="Beta",
     xlim=c(0,7), ylim=c(0, .005))
abline(h=beta, lty=2, col=gray(level=.5))
lines(c(0,mt), ciSD[,1], lty=2, col="blue")
lines(c(0,mt), ciSD[,2], lty=2, col="blue")
legend("bottomright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "smoothed state"),
       bty="n", y.intersp=1.2, cex=.7)

#delta
plot(x=NA, y=NA, type="n", xlab="time", ylab="Gamma",
     xlim=c(0,7), ylim=c(.004, .006))
abline(h=delta, lty=2, col=gray(level=.5))
lines(c(0,mt), ciSE[,1], lty=2, col="blue")
lines(c(0,mt), ciSE[,2], lty=2, col="blue")
legend("bottomright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "smoothed state"),
       bty="n", y.intersp=1.2, cex=.7)