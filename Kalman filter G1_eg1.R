# Source http://blogs2.datall-analyse.nl/2016/02/11/rcode_extended_kalman_filter/
# R Code, simulations and modeling
# by Stefan Gelissen
install.packages("dlm")
library(dlm) # YOU MUST INSTALL THESE 2 PACKAGES
library(numDeriv)

#The R functions for the extended Kalman filter (EKF) and its smoother
#can be downloaded from: 
source("http://www.datall-analyse.nl/EKF.R")
#Take a look at the functions dlmExtFilter (=the EKF) and dlmExtSmooth (=the
#RTS smoother), and notice that at the beginning of the scripts you will find
#a description of the functions' arguments.

dlmExtFilter
dlmExtSmooth

#You may save these functions FROM THE SOURCE to your computer using the save function in R.

##EXAMPLE 1

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
dataEx1 <- cbind(yprey, ypred)

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

##EXTENDED KALMAN FILTER (EKF)

#Dynamic model:
#specifying 2 states, namely [amount of prey, amount of predators].
#We will first use the dlm function from the dlm package for specifying the
#dynamic model. The dlm function checks for us if the dimensions
#of our specified matrices are correct, and if the specified covariance matrices
#are positive semidefinite.
library(dlm)
ex1 <- dlm(m0=c(400, 200), #initial state estimates for prey and predators
           #error covariances of the initial state estimates
           #this matrix reflects the uncertainty in our initial state estimates
           C0=diag(rep(100,2)),
           #observation matrix:
           #we will not use this FF matrix in the extended Kalman filter,
           #so we set the values in this matrix to zero
           FF=matrix(0, nrow=2, ncol=2),
           #measurement noise
           V=diag(c(sigmaPrey^2, sigmaPred^2)),
           #state transition matrix:
           #we will not use this GG matrix in the extended Kalman filter,
           #so we also set the values in this matrix to zero
           GG=matrix(0, nrow=2, ncol=2),
           #process noise
           W=diag(rep(0,2)))

#For the EKF we will use a list-object (instead of the dlm-object above),
#and remove the FF and GG matrix (which we do not need anyway for the EKF).
#Note that the specified initial state estimates (at t=0) below for prey and
#predators deviate from the values that were used for generating the data
#(i.e., 400 for prey and 200 for predators).
#You may change these initial state estimates too and see how they
#influence the behavior of the extended Kalman filter.
ex1 <- list(m0=c(390, 220), #initial state estimates
            #error covariances of the initial state estimates:
            #this matrix reflects the uncertainty in our initial state estimates
            #you may change the values in this matrix and see how they influence
            #the behavior of the Kalman filter
            C0=diag(rep(100,2)),
            #measurement noise
            V=diag(c(sigmaPrey^2, sigmaPred^2)),
            #process noise
            W=diag(rep(0,2)))
#Note that all covariances in the process noise matrix (W) are set to zero.
#This makes sense since the change in the amount of prey and predators at
#each time step is fully explained by our Lotka-Volterra model describing the
#prey-predator interaction.

#Specify the state transition function:
#note that we will use as state functions the difference equations
#given by Euler's forward method. These difference equations will yield valid 
#estimates for the amounts of prey and predators at each time step as long
#as the specified value for dT above is relatively small.
#WARNING: always use arguments x and k when specifying the GGfunction
GGfunction <- function (x, k){
  prey <- x[1]; pred <- x[2]
  c(prey + (alpha*prey - beta*prey*pred)*dT,
    pred + (delta*prey*pred - gamma*pred)*dT)}

#Specify the observation/measurement function
#WARNING: always use arguments x and k when specifying the FFfunction
FFfunction <- function (x, k){
  prey <- x[1]; pred <- x[2]
  c(prey, pred)}


##Compute the filtered (a posteriori) state estimates
## dlmExtFilter below HAS ERROR ; CANNOT FIND FUNCTION

ext1 <- dlmExtFilter(y=dataEx1, mod=ex1,
                     GGfunction=GGfunction, FFfunction=FFfunction)

#Instead of relying on a numerical method for approximating the Jacobians
#in the EKF, it is also possible to calculate the Jacobians by hand and
#subsequently use these in the EKF.
#WARNING: always use arguments x and k when specifying the GGjacobian
GGjacobian <- function (x, k){
  prey <- x[1]; pred <- x[2]
  c(1 + alpha*dT - beta*pred*dT, -beta*prey*dT,
    delta*pred*dT, 1 + delta*prey*dT - gamma*dT)}

#WARNING: always use arguments x and k when specifying the FFjacobian
FFjacobian <- function (x, k){
  prey <- x[1]; pred <- x[2]
  c(1, 0,
    0, 1)}

#Use these latter Jacobians in the EKF
ext1 <- dlmExtFilter(y=dataEx1, mod=ex1,
                     GGfunction=GGfunction, FFfunction=FFfunction,
                     GGjacobian=GGjacobian, FFjacobian=FFjacobian)

#plot the filtered state estimates
plot(x[,1], x[,2], type="l", lwd=2, xlab="Prey", ylab="Predators",
     xlim=range(yprey), ylim=range(ypred), col=gray(level=.5))
points(yprey, ypred, col="red", cex=.5)
lines(ext1$m[,1], ext1$m[,2], lty=2, col="blue", lwd=2)
legend("topright", pch=c(NA, 1, NA), lty=c(1, NA, 2), lwd=c(2, NA, 2),
       col=c(gray(level=.5), "red", "blue"),
       legend=c("true state", "measured", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)

#compute the confidence intervals for the filtered state estimates
varcovFilteredState <- dlmSvd2var(ext1$U.C, ext1$D.C)
#95% ci for prey
seFY <- sqrt(unlist(lapply(varcovFilteredState, function(x) x[1,1])))
ciFY <- ext1$m[,1] + qnorm(.05/2)*seFY%o%c(1, -1)
#95% ci for predators
seFD <- sqrt(unlist(lapply(varcovFilteredState, function(x) x[2,2])))
ciFD <- ext1$m[,2] + qnorm(.05/2)*seFD%o%c(1, -1)

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

##Diagnostics

#Diagnostics below are computed for the prey estimates only
#(diagnostics for the predator estimates are calculated in a similar way)

#raw residuals (=innovations)
plot(mt, yprey-ext1$f[,1], type="l", col="gray", xlab="time", ylab="Raw residuals")
abline(h=0, lty=2, col="blue")

#standardized residuals
varcovPredState <- dlmSvd2var(ext1$U.R, ext1$D.R)
Qt <- sapply(1:nm, function (i)
  ex1$V + ext1$dFF.dx[[i]]%*%varcovPredState[[i]]%*%t(ext1$dFF.dx[[i]]),
  simplify=FALSE)
sqrtQTx <- sqrt(unlist(lapply(Qt, function(x) x[1,1])))
serrx <- (yprey-ext1$f[,1])/sqrtQTx
plot(mt, serrx, type="l", col="gray", xlab="time", ylab="Standardized residuals")
abline(h=0, lty=2, col="blue")

#qq-plot for standardized residuals
qqnorm(serrx, cex=.5)
abline(a=0, b=1, col="darkgray", lty=2, lwd=2)

#acf plot for standardized residuals
acf(serrx, lag.max=10, main="Standardized residuals")


##Compute smoothed state estimates

#apply fixed-interval smoothing
smoothState <- dlmExtSmooth(ext1)
plot(x[,1], x[,2], type="l", lwd=2, xlab="Prey", ylab="Predators",
     xlim=range(yprey), ylim=range(ypred), col=gray(level=.5))
lines(ext1$m[,1], ext1$m[,2], type="l", lty=2, col="blue")
lines(smoothState$s[,1], smoothState$s[,2], type="l", lty=2, col="darkgreen")
legend("topright", lty=c(1, 2, 2), lwd=c(1, 1, 1),
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

# End of example 1
