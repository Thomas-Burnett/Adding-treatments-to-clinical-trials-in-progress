###Examples adding arms to trials in progress

##required packages
require(MAMS)
require(parallel)
require(mvtnorm)
##source functions
funloc <- ()
source(funloc)

#####################################################################################################

#Defining parameters
#The existing MAMS trial
#choose the parameters
#the number of treatments
K <- 2
#vecotr of the number of treatments in each of the 
X <- c(1:K)
#number of interim analyses
J=3
#vectors of allocation ratios for each experimental treatment
r=1:3
#ratios on control
r0=1:3
#target FWER
alpha=0.05
#power
power = 0.9
#effect size parameters for the trial
p = 0.75
p0 = 0.5

#####################################################################################################

#set up the original MAMS trial

#maximum intersection to set the lower boundaries 
trial = mamsc(K=K,
              J=J,
              alpha=alpha,
              r=r,
              r0=r0,
              power=power,
              p=p,
              p0=p0,
              ushape="triangular",
              lshape="triangular") 
lfix = trial$l[-J]

#remove the maximum an compute the boundaries for the remaining tests
X = X[-K]

#construct all other tests
#set up parallelisation cores
#Calculate the number of cores
#note the parallelisation isn't needed for low numbers of treatments
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

#what are all the required components of a closed MAMS trial

htests <- parLapply(cl=cl,
                     X = X,
                     fun = mamsc,
                     J=J,
                     alpha=alpha,
                     r=r,
                     r0=r0,
                     power=power,
                     p=p,
                     p0=p0,
                     ushape="triangular",
                     lshape="closed",
                     lfix=lfix,
                     sample.size=FALSE)

# Stop the cluster
stopCluster(cl)

#add in the full intersection test
htests[[(length(htests)+1)]] = trial

#####################################################################################################
#conditional adding

#the number of arms to be added
Kadd <- 2

#some example observations
z <- data.frame(tr1=c(2),tr2=c(1.5))

#defining all required hypothesis tests

#matrix of sample sizes
rtestsm <- t(combn(K+Kadd,1,FUN=tabulate,nbins=K+Kadd))
for(i in 2:(K+Kadd)){
  rtestsm <- rbind(rtestsm,t(combn(K+Kadd,i,FUN=tabulate,nbins=K+Kadd)))
}
rtests <- lapply(seq_len(nrow(rtestsm)), function(i) rtestsm[i,])
#we will loop through rtests in parallel

#what test is the full intersection
test <- rtests[[length(rtests)]]
#full intersection off the null hypotheses conditional test (from which to take the lower testing boundaries)
ftestc <- condtests(test,htests,z,K,Kadd,J,alpha,r,r0,power,p,p0,ushape="triangular",lshape="triangular",lfix=NULL,funloc)
#what are the lower boundaries
lfix = ftestc$l
lfix = lfix[-length(lfix)]

#remove the full intersection from the set and compute all other necessary boundaries
comptests = rtests
comptests[[length(rtests)]] = NULL

# Initiate cluster
cl = makeCluster(no_cores)

#what are all the required components of a closed MAMS trial
ctests = parLapply(cl=cl,
                    X = comptests,
                    fun = condtests,
                    htests,z,K,Kadd,J,alpha,r,r0,power,p,p0,ushape="triangular",lshape="closed",lfix=lfix,funloc)

# Stop the cluster
stopCluster(cl)

#add the full intersection
ctests[[length(rtests)]] = ftestc

#####################################################################################################

##simulating operating characteristics
#configuration of treatments, c(0.5,0.5,0.5,0.5) is the global null
pv = c(0.5,0.5,0.5,0.5)
#number of simulations to use (note this is low due to computational intensity)
nsim = 1000

#simulating the expected behavior of the conditionally updated trial under a given configuration
simcond(nsim=nsim,rtestsm,ctests,htests[[length(htests)]],z,pv)

#simulating the full trial with a conditional update, note this uses parallel computation
cl = makeCluster(no_cores)

simfullpar(pv,htests,Kadd,alpha,funloc,nsim=nsim,p=p,p0=p0,cl)

stopCluster(cl)