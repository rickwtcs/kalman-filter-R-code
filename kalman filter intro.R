#Example of a Kalman filter for estimating a fixed value with 
#measurement error from the Welch and Bishop's "An Introduction 
#to the Kalman Filter"  University of North Carolina at Chapel Hill, 
#Department of Computer Science TR 95-041

#Code ported from: Andrew D. Straw's implementation
# at: wiki.scipy.org/Cookbook/KalmanFiltering

# https://gist.github.com/mathew-hall/2ca753c68a594e2c37b1


count = 50
true_value = -0.37727 #actual value
z = rnorm(count,mean=true_value, sd=0.1) #observations


Q=1e-5 #process variance

#Allocate space:
xhat = rep(0,count) #a posteri estimate at each step
P = rep(0,count)  #a posteri error estimate
xhatminus=rep(0,count) #a priori estimate
Pminus=rep(0,count) #a priori error estimate
K=rep(0,count) #gain

#estimate of measurement variance
R = 1**2

#initialise guesses: assume true_value=0, error is 1.0
xhat[1] <- 0
P[1] <- 1

for (k in 2:count){
	#time update
	xhatminus[k] <- xhat[k-1]
	Pminus[k] <- P[k-1] + Q
	
	#measurement update
	K[k] = Pminus[k] / (Pminus[k] + R)
	xhat[k] = xhatminus[k] + K[k] * (z[k] - xhatminus[k])
	P[k] = (1-K[k]) * Pminus[k]
}


par(mfrow=c(2,1))

plot(xhat, 
		type="l", 
		col="blue", 
		xlab="Iteration", 
		ylab="Volts",
		ylim=1.05 * c(min(z,true_value,xhat),max(z,true_value,xhat))
		)
points(z,  pch=2)

abline(true_value,0)

plot(Pminus[2:count],type="l",xlab="Iteration", ylab=substitute(Variance~(Volts^2 )))