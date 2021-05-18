###adding arms to a simple clinical trial

##we require the multiivariate normal functions
require(mvtnorm)

#########################################################################################################

##Construct a Dunnett intersection p-value, this is specific to the 2 hypothesis case
##inputs are
#	maxz = maximum z-value
#	rho  = correlation between the two endpoints
##output is the dunnet p-value for a given maximum and correltion matrix
dunp <- function(maxZ,rho)
{
	1-pmvnorm(upper=c(maxZ,maxZ),corr = matrix(c(1,rho,rho,1),nrow=2))
}
#vectorize this for future use
dunpv <- Vectorize(dunp, vectorize.args = "maxZ")

#########################################################################################################

##Directly compute the probabilities of rejecting the Dunnett proecdure for intersection test
##inputs are
#	alphaplot = alpha to test the intersection at
#	rho	    = correlation between the two endpoints
#	xi1	    = true parameter value for Z1
#	xi2	    = true parameter value for Z2
#	tau	    = proportion of first sample before interim analysis
##output is the probability of rejection
dunpow <- function(alphaplot,rho,xi1,xi2,tau)
{
	#compute the rejection region for the Dunnett test
	dunrej <- qmvnorm(1-alphaplot[1],corr = matrix(c(1,rho,rho,1),nrow=2))$quantile
	#compute the probability of rejection
	condp  <- 1-pmvnorm(upper=c(dunrej,dunrej),mean=c(xi1*sqrt(1-tau),xi2),
		corr = matrix(c(1,rho,rho,1),nrow=2))
	#output the probability
	condp[[1]]
}
#vectorize for future use
dunpowv <- Vectorize(dunpow, vectorize.args = "alphaplot")

#########################################################################################################

##compute probabilities of rejection under different testing schemes
##inputs are
#	xi1	= ratio of treatment effect and sample size for treatment 1
#	xi2	= ratio of treatment effect and sample size for treatment 2
#	tau	= proportion of first sample before interim analysis
#	rho	= correlation between first and second treatment effect in stage 2
#	alpha = nominal type 1 and familywise error rate
#	n	= sample size for simulation
##output is null, probabilities output as messages
rejprobs <- function(xi1,xi2,tau,rho,alpha,n)
{
	message(paste("xi1 = ",xi1,"xi2 = ",xi2))	

	#the first treatment
	#simulate the stage 1 data
	Z1_1 = rnorm(n,xi1*sqrt(tau),1)
	#simulate the stage 2 data
	Z1_2 = rnorm(n,xi1*sqrt(1-tau),1)
	#the combined p-value
	Z1_c = sqrt(tau)*Z1_1 +sqrt(1-tau)*Z1_2
	
	#conditional error rate for first null
	alphap = 1-pnorm(sqrt(2)*qnorm(1-alpha)-Z1_1,0,1)
	
	#the introduced treatment
	Z2 = rnorm(n,xi2 + (Z1_2-xi1*sqrt(1-tau))/2,sqrt(1-rho^2))
	
	#max z-value to apply Dunnett test to
	maxZ <- Z1_2*(Z1_2 >= Z2) + Z2*(Z1_2 < Z2)
	#find each of the p-values
	P1_2 <- 1 - pnorm(Z1_2)
	P2   <- 1 - pnorm(Z2)
	P12  <- dunpv(maxZ,rho)
	
	#probabilities of rejecting component hypotheses of closed testing procedure
	message("Probabilities of rejecting the components of the closed testing procedure")
	message(c("    The probability of rejecting H_01 is: ",mean(P1_2 < alphap)))
	message(c("    The probability of rejecting H_02 is: ",mean(P2 < alpha)))
	message(c("    The probability of rejecting H_012 at alphap is: ",mean(P12 < alphap)))
	message(c("    The FWER is: ",
		mean((xi1<=0)*(xi2<=0)*(P1_2 < alphap)*(P2 < alpha)*(P12 < alphap)+
		(xi1<=0)*(P1_2 < alphap)*(P12 < alphap)+
		(xi2<=0)*(P2 < alpha)*(P12 < alphap))
	))
	message(" ")

	#probabilities of rejection, new testing procedure
	R1 <- mean((P1_2 < alphap)*(P2 > alpha)*(P12 < alphap))
	R2 <- mean((P1_2 > alphap)*(P2 < alpha)*(P12 < alphap))
	RB <- mean((P1_2 < alphap)*(P2 < alpha)*(P12 < alphap))
	
	message("Probabilities of rejection under proposed testing procedure (alphap for intersection)")
	message(c("Probability of rejecting H_01 only: ",R1))
	message(c("Probability of rejecting H_02 only: ",R2))
	message(c("Probability of rejecting both: ",RB))
	message(c("Probability of rejecting any: ",R1+R2+RB))
	message(" ")
	
	#probabilities of rejection when testing intersection based only on stage 2 at alpha 
	R1 <- mean((P1_2 < alphap)*(P2 > alpha)*(P12 < alpha))
	R2 <- mean((P1_2 > alphap)*(P2 < alpha)*(P12 < alpha))
	RB <- mean((P1_2 < alphap)*(P2 < alpha)*(P12 < alpha))
	
	message("Probabilities of rejection using alpha for intersection)")
	message(c("Probability of rejecting H_01 only: ",R1))
	message(c("Probability of rejecting H_02 only: ",R2))
	message(c("Probability of rejecting both: ",RB))
	message(c("Probability of rejecting any: ",R1+R2+RB))
	message(" ")
	
	#probabilities of rejection when testing intersection based only on H_{01}
	R1 <- mean((P1_2 < alphap)*(P2 > alpha))
	R2 <- NA
	RB <- mean((P1_2 < alphap)*(P2 < alpha))
	
	message("Probabilities of rejection under H_01 gatekeeping procedure")
	message(c("Probability of rejecting H_01 only: ",R1))
	message(c("Probability of rejecting H_02 only: ",R2))
	message(c("Probability of rejecting both: ",RB))
	message(c("Probability of rejecting any: ",R1+RB))
}
#Vecotrize this function by xi1 and xi2 to make long runs easier
rejprobsv <- Vectorize(rejprobs,vectorize.args = c("xi1","xi2"))

#########################################################################################################

##produce plots of conditional error vs probability of rejection
##inputs are
#	xi1	= ratio of treatment effect and sample size for treatment 1
#	xi2	= ratio of treatment effect and sample size for treatment 2
#	main  = titles for treatment effects as expressions
#	m 	= number of points to be used to produce the plots
#	rho	= correlation between first and second treatment effect in stage 2
#	alpha = nominal type 1 and familywise error rate
#	tau	= proportion of first sample before interim analysis
##output is a plot demonstrating conditional probabilities of rejection for corresponding z-values
plotcond <- function(xi1,xi2,main,m,rho,alpha,tau)
{
	#vectors of observations to plot over
	#90% coverage for the Z-value
	Z1plot <- seq(from=qnorm(0.05,xi1*sqrt(tau),1),to=qnorm(0.95,xi1*sqrt(tau),1),length.out=m)
	#corresponding conditional errors
	alphaplot = 1-pnorm(sqrt(1/tau)*qnorm(1-alpha)-Z1plot,0,1)
	
	#compute corresponding rejection probabilities
	intrej <- dunpowv(alphaplot,rho,xi1,xi2,tau)
	
	#produce the plot
	#leave appropriate space for the labels
	par(mar = c(6, 6, 3, 6))
	#plot prob of rejection
	plot(alphaplot,intrej,type="l",lwd=2,col="red",xlab="",ylab="",ylim=c(0,1))
	#maintain same plot window
	par(new=TRUE)
	#plot density function on same plot
	plot(alphaplot,dnorm(Z1plot,xi1/sqrt(2),1),type="l",
		lty=2,col="gray",xaxt = "n", yaxt = "n",ylab = "", xlab = "")
	#add in the axis for density on the opposite side
	axis(side = 4)

	#add all the labels, note the first is dictated by main and should include the treatment effects
	mtext(main, side = 3,line=1)
	mtext(expression(f( z[1]^(1) )), side = 4, line = 3)
	mtext(expression(A( z[1]^(1) )), side = 1, line = 3)
	mtext(expression(P(Reject ~ H[O1] ~ "|" ~ z[1]^(1))),side = 2, line = 3)
}
#vectorize xi1, xi2 and the treatment effect labels
plotcondv <- Vectorize(plotcond, vectorize.args = c("xi1","xi2","main"))

