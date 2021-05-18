#MAMS arm adding functions

#where possible we will make use of the pre-existing MAMS package and functions 
require(mvtnorm)

###conderror
##function for computing the conditional FWER 
#inputs are
#	emams = the details of the existing MAMS trial
#	obs	= the observations to date in the trial that has been conducted so far
# N defines the mesh for the numerical integration
#outputs are
#	the conditional FWER under the global null
#	note the FWER is maximised at this point meaning all sub tests will be conservative
conderror = function(emams,obs,N = 20)
{
	#any necessary unpacking
	K  = emams$K
	J  = emams$J
	r  = emams$r[2,]
	r0 = emams$r[1,] 
	u  = emams$u
	l  = emams$l
	z  = obs

	#we require some standard functions from the MAMS package
	#with minor adjustment to allow for the conditioning on previous obs
	######
	 mesh <- function(x, d, w = 1/length(x) + x * 0) {
  	      n <- length(x)
   	     W <- X <- matrix(0, n^d, d)
  	      for (i in 1:d) {
   	         X[, i] <- x
  	          W[, i] <- w
  	          x <- rep(x, rep(n, length(x)))
 	           w <- rep(w, rep(n, length(w)))
  	      }
   	     w <- exp(rowSums(log(W)))
   	     list(X = X, w = w)
  	  }
	#products are no longer over the same thing for each endpoint, each one has a unique boundary 
	#since we never made use of the symetry (except computationally)
	#this is why this is maximised under the global null 
   	 prodsum <- function(x, lc, uc, r, r0, r0diff, J, K, Sigma) {
		#product over density functions for integration
      	  int <- prod(sapply(x, dnorm))

		#each conditional set of error boundaries will correspond to  different L and U
		for(i in 1:K)
		{
		  l <- lc[,i]
		  u <- uc[,i]
	        L <- sqrt(r[1])/r0[1] * sqrt(r0diff[1]) * x[1] + l[1] * 
	            sqrt(1 + r[1]/r0[1])
	        U <- sqrt(r[1])/r0[1] * sqrt(r0diff[1]) * x[1] + u[1] * 
	            sqrt(1 + r[1]/r0[1])
	        insum <- pnorm(L)
	        if (J > 1) {
      	      for (j in 2:J) {
	                L[j] <- sqrt(r[j])/r0[j] * (sqrt(r0diff[1:j]) %*% 
	                  x[1:j]) + l[j] * sqrt(1 + r[j]/r0[j])
	                U[j] <- sqrt(r[j])/r0[j] * (sqrt(r0diff[1:j]) %*% 
	                  x[1:j]) + u[j] * sqrt(1 + r[j]/r0[j])
	                insum <- insum + pmvnorm(lower = c(L[1:j - 1], 
	                  -Inf), upper = c(U[1:j - 1], L[j]), sigma = Sigma[1:j, 
	                  1:j])[1]
	            }
	        }
	        int <- int * insum
		}
	        return(int)
	    }
	#####
	##now from these we compute the probability of a type 1 error given exisiting obs
	#first find the conditional error boundaries 
	
	#what stage of the trial has we observed up to 
	os <- nrow(z)

	#compute weight of current data 
	w1 <- rep(0,(J-os)) 
	for(i in 1:J-os){
		w1[i] <- sqrt((r0[os] + r[os])/(r0[os+i]+r[os+i]))
	}

	#hence the weight for the planned data in reference trial
	w2 <- sqrt(1-w1*w1)

	#what are the remaining testing boundaries we are interest in
	u <- u[-(1:os)]
	l <- l[-(1:os)]	

	#conditional upper testing boundaries, for each experimental treatment
	uc <- data.frame(matrix(u,ncol=ncol(z),nrow=(J-os)))
	for(i in 1:ncol(z))
	{
		uc[,i] <- (u-w1*z[os,i])*(1/w2)
	}
	
	#conditional lower testing boundaries, for each experimental treatment
	lc <- data.frame(matrix(l,ncol=ncol(z),nrow=(J-os)))
	for(i in 1:ncol(z))
	{
		lc[,i] <- (l-w1*z[os,i])*(1/w2)
	}

	##given conditional testing boundaries what is conditional error probability
	#compute some crucial facts for the integral from the original trial setup
	#taken from mams funciton in the mams package, these are the essential details
	#in setting up the covariance matrix for the computation
	#####
	bottom <- matrix(r, J, J)
	top <- matrix(rep(r, rep(J, J)), J, J)
	top[upper.tri(top)] <- t(top)[upper.tri(top)]
	bottom[upper.tri(bottom)] <- t(bottom)[upper.tri(bottom)]
	Sigma <- sqrt(top/bottom)
	r0lag1 <- c(0, r0[1:J - 1])
	r0diff <- r0 - r0lag1
	#####

	#mesh size for integration, note we only need mesh for stages left to come
	mmp <- mesh((1:N - 0.5)/N * 12 - 6, J-os, rep(12/N, N))
	#evaluate value of the function at each integration point
	evs <- apply(mmp$X, 1, prodsum, lc = lc, uc = uc, r = r, 
		r0 = r0, r0diff = r0diff, J = J-os, K = K, Sigma = Sigma)
	condalpha <- 1 - mmp$w %*% evs

	#output conditional error
	as.numeric(condalpha)
}


########################################################################################################
###condmams
##function for computing MAMS testing boundaries adjusting for exisiting results
##note this is useful for the unrestricted case where the conditional error is lower 
##than the original alpha 
##inputs are 
#	emams = the details of the existing MAMS trial
#	obs	= the observations to date in the trial that has been conducted so far
# Kadd = the number of treatments to be added to the trial
# alpha = the target familywise error rate
# the other inputs match those given for mams() function in the mams package in R
##outputs as class MAMS are
# l = lower bounds to be used
# u = upper bounds to be used
# n = sample size
# rMat = matrix of randomisation proportions for the remainder of the trial
# N = the maximum possible sample size
# K = the total number of arms in the rest of the trial (including the added arms)
# J = the remaining number of stages of the trial
# alpha = the familywise error rate
condmams <- function(emams,obs,Kadd,alpha, p = 0.75, p0 = 0.5, delta = NULL, delta0 = NULL, sd = NULL, 
                     ushape = "obf", lshape = "fixed", ufix = NULL, lfix = 0, 
                     nstart = 1, nstop = NULL, sample.size = TRUE, N = 20, type = "normal")
{
  #some unpacking
  K  = emams$K
  J  = emams$J
  r  = emams$r[2,]
  r0 = emams$r[1,] 
  u  = emams$u
  l  = emams$l
  n  = emams$n
  z  = obs
  
  #we require some standard functions from the MAMS package
  #with minor adjustment to allow for the conditioning on previous obs
  ######
  mesh <- function(x, d, w = 1/length(x) + x * 0) {
    n <- length(x)
    W <- X <- matrix(0, n^d, d)
    for (i in 1:d) {
      X[, i] <- x
      W[, i] <- w
      x <- rep(x, rep(n, length(x)))
      w <- rep(w, rep(n, length(w)))
    }
    w <- exp(rowSums(log(W)))
    list(X = X, w = w)
  }
  
  prodsum <- function(x, lc, uc, r, r0, r0diff, J, K, Sigma) {
    #product over density functions for integration
    int <- prod(sapply(x, dnorm))
    
    #each conditional set of error boundaries will correspond to  different L and U
    for(i in 1:K)
    {
      l <- lc[,i]
      u <- uc[,i]
      L <- sqrt(r[1])/r0[1] * sqrt(r0diff[1]) * x[1] + l[1] * 
        sqrt(1 + r[1]/r0[1])
      U <- sqrt(r[1])/r0[1] * sqrt(r0diff[1]) * x[1] + u[1] * 
        sqrt(1 + r[1]/r0[1])
      insum <- pnorm(L)
      if (J > 1) {
        for (j in 2:J) {
          L[j] <- sqrt(r[j])/r0[j] * (sqrt(r0diff[1:j]) %*% 
                                        x[1:j]) + l[j] * sqrt(1 + r[j]/r0[j])
          U[j] <- sqrt(r[j])/r0[j] * (sqrt(r0diff[1:j]) %*% 
                                        x[1:j]) + u[j] * sqrt(1 + r[j]/r0[j])
          insum <- insum + pmvnorm(lower = c(L[1:j - 1], 
                                             -Inf), upper = c(U[1:j - 1], L[j]), sigma = Sigma[1:j, 
                                                                                               1:j])[1]
        }
      }
      int <- int * insum
    }
    return(int)
  }
  #this function requires the most 
  typeI <- function(C, alpha, N, r, r0, r0diff, J, K, z, Sigma, 
                    ushape, lshape,Jc,rc,r0c, lfix = NULL, ufix = NULL) {
    
    #set boundaries based on the supplied root value
    #note this is based on the stages left to be run
    if (!is.function(ushape)) {
      if (ushape == "obf") {
        u <- C * sqrt(rc[Jc]/rc)
      }
      else if (ushape == "pocock") {
        u <- rep(C, Jc)
      }
      else if (ushape == "fixed") {
        u <- c(rep(ufix, Jc - 1), C)
      }
      else if (ushape == "triangular") {
        u <- C * (1 + rc/rc[Jc])/sqrt(rc)
      }
    }
    else {
      u <- C * ushape(Jc)
    }
    if (!is.function(lshape)) {
      if (lshape == "obf") {
        l <- c(-C * sqrt(rc[Jc]/rc[1:(Jc - 1)]), u[Jc])
      }
      else if (lshape == "pocock") {
        l <- c(rep(-C, Jc - 1), u[Jc])
      }
      else if (lshape == "fixed") {
        l <- c(rep(lfix, Jc - 1), u[Jc])
      }
      else if (lshape == "closed") {
        l <- c(lfix,u[Jc])
        #for optimisation need to ensure u>l
        u <- (u>=l)*u + (u<l)*l
      }
      else if (lshape == "triangular") {
        if (ushape == "triangular") {
          l <- -C * (1 - 3 * rc/rc[Jc])/sqrt(rc)
        }
        else {
          l <- -C * (1 - 3 * rc/rc[Jc])/sqrt(rc)/(-1 * (1 - 
                                                          3)/sqrt(Jc))
        }
      }
    }
    else {
      l <- c(C * lshape(Jc)[1:(Jc - 1)], u[Jc])
    }
    
    #up to this point we have used the same structure as the unconditional function
    #now we must add the conditional error boundaries for the existing treatments
    
    ###############################################
    
    #what stage of the trial has we observed up to 
    os <- nrow(z)
    
    #compute weight of current data 
    w1 <- rep(0,(J-os)) 
    for(i in 1:J-os){
      w1[i] <- sqrt((r0[os] + r[os])/(r0[os+i]+r[os+i]))
    }
    
    #hence the weight for the planned data in reference trial
    w2 <- sqrt(1-w1*w1)	
    
    #conditional upper testing boundaries, for each experimental treatment
    uc <- data.frame(matrix(u,ncol=ncol(z),nrow=(J-os)))
    for(i in 1:ncol(z))
    {
      uc[,i] <- (u-w1*z[os,i])*(1/w2)
    }
    
    #conditional lower testing boundaries, for each experimental treatment
    lc <- data.frame(matrix(l,ncol=ncol(z),nrow=(J-os)))
    for(i in 1:ncol(z))
    {
      lc[,i] <- (l-w1*z[os,i])*(1/w2)
    }
    
    #combine all the required error boundaries, conditional and added
    #the nrow adds a new column of unadjusted boundaries for each endpoint to be added
    up  <- data.frame(uc,matrix(u,ncol=(Kadd),nrow=Jc))
    lp  <- data.frame(lc,matrix(l,ncol=(Kadd),nrow=Jc))
    
    ##given conditional testing boundaries what is conditional error probability
    #compute some crucial facts for the integral from the original trial setup
    #taken from mams funciton in the mams package, these are the essential details
    #in setting up the covariance matrix for the computation
    #####
    #note this only concerns the remainder of the trial now
    bottom <- matrix(rc, Jc, Jc)
    top <- matrix(rep(rc, rep(Jc, Jc)), Jc, Jc)
    top[upper.tri(top)] <- t(top)[upper.tri(top)]
    bottom[upper.tri(bottom)] <- t(bottom)[upper.tri(bottom)]
    Sigma <- sqrt(top/bottom)
    r0lag1 <- c(0, r0c[1:Jc - 1])
    r0diff <- r0c - r0lag1
    #####
    
    #mesh size for integration, note we only need mesh for stages left to come
    mmp <- mesh((1:N - 0.5)/N * 12 - 6, J-os, rep(12/N, N))
    #evaluate value of the function at each integration point
    evs <- apply(mmp$X, 1, prodsum, lc = lp, uc = up, r = rc, 
                 r0 = r0c, r0diff = r0diff, J = J-os, K = K+Kadd, Sigma = Sigma)
    truealpha <- 1 - mmp$w %*% evs
    return(truealpha - alpha)
    
    ###############################################
  }
  #####
  #this should allow the computation of updated testing boundaries for the conditional procedure
  
  ###Error checking
  ##kept from the basic mams funtion
  
  if (K%%1 > 0 | J%%1 > 0) {
    stop("K and J need to be integers.")
  }
  if (K < 1 | J < 1) {
    stop("The number of stages and treatments must be at least 1.")
  }
  if (N <= 3) {
    stop("Number of points for integration by quadrature to small or negative.")
  }
  if (N > 3 & N <= 10) {
    warning("Number of points for integration by quadrature is small which may result in inaccurate solutions.")
  }
  if (alpha < 0 | alpha > 1 ) {
    stop("Error rate not between 0 and 1.")
  }
  if (length(r) != length(r0)) {
    stop("Different length of allocation ratios on control and experimental treatments.")
  }
  if (length(r) != J) {
    stop("Length of allocation ratios does not match number of stages.")
  }
  if (is.numeric(p) & is.numeric(p0) & is.numeric(delta) & 
      is.numeric(delta0) & is.numeric(sd)) {
    stop("Specify the effect sizes either via (p, p0) or via (delta, delta0, sd) and set the other parameters to NULL.")
  }
  if (is.numeric(p) & is.numeric(p0)) {
    if (p < 0 | p > 1 | p0 < 0 | p0 > 1) {
      stop("Treatment effect parameter not within 0 and 1.")
    }
    if (p <= p0) {
      stop("Interesting treatment effect must be larger than uninteresting effect.")
    }
    if (p0 < 0.5) {
      warning("Uninteresting treatment effect less than 0.5 which implies that reductions in effect over placebo are interesting.")
    }
  }
  else {
    if (is.numeric(delta) & is.numeric(delta0) & is.numeric(sd)) {
      if (sd <= 0) {
        stop("Standard deviation must be positive.")
      }
    }
    else {
      stop("Specify the effect sizes either via (p, p0) or via (delta, delta0, sd).")
    }
  }
  if (is.function(ushape) & is.function(lshape)) {
    warning("You have specified your own functions for both the lower and upper boundary. Please check carefully whether the resulting boundaries are sensible.")
  }
  #####End of error checking, we are now ready to apply the new design	
  
  #the boundaries should be set based only on the number of stages left to conduct
  os <- nrow(z)
  rc  <- r[(os+1):J]-r[os]
  r0c <- r0[(os+1):J]-r0[os]
  Jc <- J-os
  
  if (!is.function(ushape)) {
    if (!ushape %in% c("pocock", "obf", "triangular", "fixed")) {
      stop("Upper boundary does not match the available options.")
    }
    if (ushape == "fixed" & is.null(ufix)) {
      stop("ufix required when using a fixed upper boundary shape.")
    }
  }
  else {
    b <- ushape(Jc)
    if (!all(sort(b, decreasing = TRUE) == b)) {
      stop("Upper boundary shape is increasing.")
    }
  }
  if (!is.function(lshape)) {
    if (!lshape %in% c("pocock", "obf", "triangular", "closed", "fixed")) {
      stop("Lower boundary does not match the available options.")
    }
    if (lshape == "fixed" & is.null(lfix)) {
      stop("lfix required when using a fixed lower boundary shape.")
    }
  }
  else {
    b <- lshape(Jc)
    if (!all(sort(b, decreasing = FALSE) == b)) {
      stop("Lower boundary shape is decreasing.")
    }
  }
  if (is.numeric(p) & is.numeric(p0)) {
    delta <- sqrt(2) * qnorm(p)
    delta0 <- sqrt(2) * qnorm(p0)
    sig <- 1
  }
  else {
    delta <- delta
    delta0 <- delta0
    p0 <- pnorm(delta0/sqrt(2 * sd^2))
    sig <- sd
  }
  h <- min(c(rc, r0c))
  rc <- rc/h
  r0c <- r0c/h
  
  bottom <- matrix(rc, Jc, Jc)
  top <- matrix(rep(rc, rep(Jc, Jc)), Jc, Jc)
  top[upper.tri(top)] <- t(top)[upper.tri(top)]
  bottom[upper.tri(bottom)] <- t(bottom)[upper.tri(bottom)]
  Sigma <- sqrt(top/bottom)
  r0lag1 <- c(0, r0c[1:Jc - 1])
  r0diff <- r0c - r0lag1
  if (!is.function(lshape)) {
    if (Jc == 1 & lshape == "obf") {
      lshape <- "pocock"
    }
  }
  
  #place some limits on the possible uJ
  root <- max(qnorm(1 - alpha)/2,0.4)
  uJ <- NULL
  try(uJ <- uniroot(typeI, c(root, 5), alpha = alpha, 
                    N = N, r = r, r0 = r0, r0diff = r0diff, J = J, K = K, z = z, 
                    Sigma = Sigma, ushape = ushape, lshape = lshape,Jc=Jc,rc=rc,r0c=r0c, lfix = lfix, 
                    ufix = ufix, tol = 0.001)$root, silent = TRUE)
  
  if (is.null(uJ)) {
    stop("No boundaries can be found.")
  }
  if (!is.function(ushape)) {
    if (ushape == "obf") {
      u <- uJ * sqrt(rc[Jc]/rc)
    }
    else if (ushape == "pocock") {
      u <- rep(uJ, Jc)
    }
    else if (ushape == "fixed") {
      u <- c(rep(ufix, Jc - 1), uJ)
    }
    else if (ushape == "triangular") {
      u <- uJ * (1 + rc/rc[Jc])/sqrt(rc)
    }
  }
  else {
    u <- uJ * ushape(Jc)
  }
  if (!is.function(lshape)) {
    if (lshape == "obf") {
      l <- c(-uJ * sqrt(rc[Jc]/rc[1:(Jc - 1)]), u[Jc])
    }
    else if (lshape == "pocock") {
      l <- c(rep(-uJ, Jc - 1), u[Jc])
    }
    else if (lshape == "fixed") {
      l <- c(rep(lfix, Jc - 1), u[Jc])
    }
    else if (lshape == "closed") {
      l <- c(lfix,u[Jc])
      #for optimisation need to ensure u>l
      u <- (u>=l)*u + (u<l)*l
    }
    else if (lshape == "triangular") {
      if (ushape == "triangular") {
        l <- -uJ * (1 - 3 * rc/rc[Jc])/sqrt(rc)
      }
      else {
        l <- -uJ * (1 - 3 * rc/rc[Jc])/sqrt(rc)/(-1 * (1 - 
                                                         3)/sqrt(Jc))
      }
    }
  }
  else {
    l <- c(uJ * lshape(Jc)[1:(Jc - 1)], u[Jc])
  }
  
  #find the alpha spent by each stage of the trial (conditional) 
  alpha.star <- numeric(Jc)
  alpha.star[1] <- typeI(u[1], alpha = 0, N = N, r = r, 
                         r0 = r0, r0diff = r0diff[1], J = os+1, K = K, z = z, Sigma = Sigma, 
                         ushape = "fixed", lshape = "fixed",Jc=1,rc=rc[1],r0c=r0c[1], lfix = NULL, ufix = NULL)
  
  if (Jc > 1) {
    for (j in 2:Jc) {
      alpha.star[j] <- typeI(u[j], alpha = 0, N = N, r = r, 
                             r0 = r0, r0diff = r0diff[1:j], J = os+j, K = K, z = z,
                             Sigma = Sigma, ushape = "fixed", lshape = "fixed", Jc=j,rc=rc[1:j],r0c=r0c[1:j],
                             lfix = l[1:(j - 1)], ufix = u[1:(j - 1)])
    }
  }
  #####this should compute the boundaries now output the results
  
  
  #output the results
  
  res <- NULL
  res$l <- l
  res$u <- u
  res$n <- n
  res$rMat <- rbind(r0c, matrix(rc, ncol = Jc, nrow = K+Kadd, byrow = TRUE))
  res$N <- sum(ceiling(res$rMat[, Jc] * res$n))
  res$K <- K+Kadd
  res$J <- Jc
  res$alpha <- alpha
  res$alpha.star <- alpha.star
  res$power <- NA
  res$type <- type
  class(res) <- "MAMS"
  return(res)
}	

########################################################################################################

###simerror
##simulate only the first stage of the trial to show conditional errors will intgrate to alpha as required
##inputs are
#	emams = the details of the existing MAMS trial
# delta = the treatment effect (should be zero to show the error rate as required)
##this outputs a vector of conditional error rates to be used for the remainder of the trial
simerror <- function(emams, delta, nsim = 1000, sig=1) 
{	
#the basic principle is that we simulate the whole trial
#and then see what has happened (easier than conditional simulation
	#number of trial stages
	R <- emams$r[-1,1]
	r0<-emams$r[1,1]
	u <- emams$u[1]
	l <- emams$l[1]
	n <- emams$n
	J <- 1
	K <- length(R-1)

	error <- numeric(nsim)

	for(i in 1:nsim){
		eff <- 0; fut <- 0; calpha <- 0
		#treatment means in the experimental arms
		#NOTE this is where the simulation is actually happening
		mukhats <- (rnorm(K)*sig*sqrt(n) + n*delta)/(R*n)

		#treatment effect in the control treatment
	        mu0hats <- rnorm(1, 0, sig * sqrt(n))/(r0*n)
		#compute the corresponding zvalues, 
	        zks <- (mukhats - mu0hats)/(sig * sqrt((R + r0)/(R * 
	            r0 * n)))

		#what do we do with this
		fut <- 1-(any(zks>l))
		eff <- (any(zks>u))
	
		if((1-fut) & (1-eff)){
			cond <- which(zks > l)
			temp <- emams
			temp$K <- length(cond)
			calpha <- conderror(temp,matrix(zks[cond],nrow=1))
		}
		error[i] <- eff + calpha
	}
	return(error)
}

########################################################################################################

###simcond
##function for finding out the characteristics of the trial via simulation
##we are interested in probability of rejection and 
##inputs are
# nsim = the number of simulations to be used for the remainder of the trial
# rtestsm = the randomisations as a matrix
# ctests = conditional testing boundaries
# emams = the details of the existing MAMS trial
#	obs	= the observations to date in the trial that has been conducted so far
# pv = vector of treatment effect probabilities
##outputs are
# hmat = prbabiilities of rejection of each null hypothesis
# proprej = probabilities of rejecting multiple null hypotheses
# exss = expected sample size
simcond <- function(nsim=1000,rtestsm,ctests,emams,obs,pv=rep(0.5,4),deltav=NULL,sd=NULL,ptest=1) 
{
  pv <- as.numeric(pv)
  
  #first we need to unpack what we have been supplied
  amams <- ctests[[length(ctests)]]
  
  ###from supplied boundaries find corresponding conditional boundaries
  #what stage of the trial has we observed up to 
  #general adding of variables	
  R  <- t(amams$r[-1,])
  r0 <- amams$r[1,]
  l  <- amams$l
  K  <- amams$K
  z  <- obs
  os <- nrow(z)
  Jc <- amams$J
  J  <- Jc + os
  rc <- emams$r[2,os]
  r0c<- emams$r[1,os]
  n  <- emams$n
  r  <- amams$r[2,]
  
  #upper testing boundaries for each test at each stage
  u = list()
  
  for(j in 1:(J-os)){
    u[[j]] = rtestsm
    for(i in 1:nrow(rtestsm)){
      u[[j]][i,] = rtestsm[i,] * rep(ctests[[i]]$u[j],K)
    }}
  
  nMat <- t(amams$r * n)
  
  #####Error checking
  if (is.numeric(pv) & is.numeric(deltav) & is.numeric(sd)) {
    stop("Specify the effect sizes either via pv or via deltav and sd, and set the other parameters to NULL.")
  }
  if (is.numeric(pv)) {
    if (any(pv < 0) | any(pv > 1)) {
      stop("Treatment effect parameter not within 0 and 1.")
    }
    if (length(pv) != (ncol(nMat) - 1)) 
      stop("Length of pv is not K.")
  }
  else {
    if (is.numeric(deltav) & is.numeric(sd)) {
      if (length(deltav) != (ncol(nMat) - 1)) 
        stop("Length of deltav is not K.")
      if (sd <= 0) {
        stop("Standard deviation must be positive.")
      }
    }
    else {
      stop("Specify the effect sizes either via pv or via deltav and sd.")
    }
  }
  #####end of error checking
  
  
  #####simulate a single realisation of the trial
  sim <- function(n, l, u, R, r0, delta, sig) 
  {	
    #the basic principle is that we simulate the whole trial
    #and then see what has happened (easier than conditional simulation
    #number of trial stages
    J <- dim(R)[1]
    K <- dim(R)[2]
    #change in sample prooportion between stages
    Rdiff <- R - rbind(0, R[-J, ])
    r0diff <- r0 - c(0, r0[-J])
    #treatment means in the experimental arms
    #NOTE this is where the simulation is actually happening
    #also note simple no centrality parameter 
    mukhats <- apply(matrix(rnorm(J * K), J, K) * sig * sqrt(Rdiff * 
                                                               n) + Rdiff * n * matrix(delta, nrow = J, ncol = K, 
                                                                                       byrow = TRUE), 2, cumsum)/(R * n)
    #treatment effect in the control treatment
    mu0hats <- cumsum(rnorm(J, 0, sig * sqrt(r0diff * n)))/(r0 * 
                                                              n)
    #compute the corresponding zvalues, 
    zks <- (mukhats - mu0hats)/(sig * sqrt((R + r0)/(R * 
                                                       r0 * n)))
    #for endpoints already in the trial incorporate existing information
    ###from supplied boundaries find corresponding conditional boundaries
    #what stage of the trial has we observed up to 
    #compute weight of current data 
    w1 <- rep(0,(J)) 
    for(i in 1:J){
      w1[i] <- sqrt((r0c + rc)/(r0[i]+r[i]+r0c + rc))
    }
    #hence the weight for the planned data in reference trial
    w2 <- sqrt(1-w1*w1)	
    
    #original endpoints
    zkso <- data.frame(zks[,1:dim(z)[2]])
    #weight together with exisiting data
    for(j in 1:J)
    {
      zkso[j,] <- w1[j]*z + w2[j]*zkso[j,]
    }
    #added endpoints are straight forward
    zksa <- zks[,(1+dim(z)[2]):K]
    
    zks  <- data.frame(zkso,zksa)
    
    #indicators for efficacy, futility and sample size
    eff <- 0
    fut <- 0
    ss <- 0
    #start with the first stage of the trial
    j <- 1
    #matrix of results for each arm in each stage ready to be populated
    hmat <- matrix(0, nrow = J, ncol = K)
    remaining <- rep(T, K)
    #which of the treatments are still in the trial
    #while we have not stopped either for efficacy or futility
    #or the end of the trial do this
    
    while ((eff == 0) && (fut == 0) && (j <= J)) {
      #current sample size
      ss <- sum((n * Rdiff[j, ])[remaining]) + n * r0diff[j] + ss
      
      #turn observations for the current stage into matrix (used for checking rejections)
      Zmat = t(matrix(zks[j,],nrow=length(zks[j,]),ncol=nrow(rtestsm)))
      
      #which hypotheses can be rejected globally
      #on the left indicators of which tests have been rejected for which hypotheses
      #this checks which tests have been rejected, then iff all tests for a corresponding hypotheses are rejected
      rej = (colSums(rtestsm*(rowSums(rtestsm*(Zmat > u[[j]]))>0)) == colSums(rtestsm))*remaining
      eff = sum(rej)
      
      fut <- (max(zks[j, remaining]) < l[j])
      
      hmat[j,]  = rej
      
      remaining = (zks[j, ] > l[j]) & remaining
      
      #increment stage of the trial
      j <- j + 1
      
    }
    
    #indicator if something was rejected
    rej <- eff
    #how many endpoints were rejected?
    pow <- (remaining[1]) && eff
    if (pow) {
      pow <- (which(zks[j - 1, remaining] == max(zks[j - 
                                                       1, remaining])) == 1)
    }
    if (all(remaining == 0)) {
      rem <- 0
    }
    else {
      rem <- as.numeric(which(remaining))
    }
    #return what happened in the trial
    return(list(rej = rej, pow = pow, ess = ss, hmat = hmat))
  }
  #####end of simulation of a single realisation of the trial
  
  ##a little work to set up the appropriate matricies
  r0 <- nMat[, 1]/nMat[1, 1]
  if (ncol(nMat) == 2) {
    R <- t(t(nMat[, -1]/nMat[1, 1]))
  }
  else {
    R <- nMat[, -1]/nMat[1, 1]
  }
  if (!is.matrix(R) && is.vector(R)) 
    R <- t(as.matrix(nMat[, -1]/nMat[1, 1]))
  n <- nMat[1, 1]
  if (is.numeric(pv)) {
    deltas <- sqrt(2) * qnorm(pv)
    sig <- 1
  }
  else {
    deltas <- deltav
    sig <- sd
  }
  
  ##now we can simulate multiple realisations of the trial
  reps <- sapply(rep(n, nsim), sim, l, u, R, r0, deltas, sig)
  rejdat <- matrix(numeric(nsim*K),nrow=nsim,ncol=K)
  for (i in 1:nsim) {
    rejdat[i,] <- colSums(reps["hmat", i][[1]])
  }   
  
  res <- NULL
  res$hmat <- colMeans(rejdat)
  res$proprej  <- summary(as.factor(rowSums(rejdat)))/nsim
  res$exss <- mean(unlist(reps["ess", ]))
  return(res)
}

#####################################################################################################

###simfull
##function for simulating the operating characteristics of the full trial should
##treatment arms be added during the trial
##inputs are
# pv = vector of treatment effect probabilities
# htests = all elements of the closed testing procedure for the orginal MAMS trial
# Kadd = the number of treatments to be added
# alpha = the familywise error rate
# funloc = the location of this functions script (required to parallelisation later)
# all other inputs are as described in the MAMS package
simfull <- function(pv,htests,Kadd,alpha,funloc,nsim = 1000, sd=1,
                    p = 0.75, p0 = 0.5, delta = NULL, delta0 = NULL, 
                    ushape = "obf", lshape = "fixed", ufix = NULL, lfix = 0, 
                    nstart = 1, nstop = NULL, sample.size = TRUE, N = 20, type = "normal") 
{	
  #the basic principle is that we simulate the whole trial
  #and then see what has happened (easier than conditional simulation
  emams <-  htests[[length(htests)]]
  #number of trial stages
  R <- emams$r[-1,1]
  r0 <-emams$r[1,1]
  u <- emams$u[1]
  l <- emams$l[1]
  n <- emams$n
  alpha <- emams$alpha
  J <- 1
  K <- emams$K
  addindex <- c((K+1):(K+Kadd))
  sig <- sd
  
  #upper testing boundaries for stage 1
  temp <- t(combn(K,1,FUN=tabulate,nbins=K))*htests[[1]]$u[1]
  #stage 1 required tests
  s1tests <- t(combn(K,1,FUN=tabulate,nbins=K))
  for(i in 2:K){
    temp <- rbind(temp,t(combn(K,i,FUN=tabulate,nbins=K))*htests[[1]]$u[1])
    s1tests <- rbind(s1tests,t(combn(K,i,FUN=tabulate,nbins=K)))
  }
  u <- temp
  rm(temp)
  
  if (is.numeric(pv)) {
    deltav <- sqrt(2) * qnorm(pv)
    sig <- 1
  }
  
  trials <- NULL
  for(i in 1:nsim){
    eff <- 0; fut <- 0; calpha <- 0
    #treatment means in the experimental arms
    #NOTE this is where the simulation is actually happening
    
    #need to sumulate full set of muhats
    mukhats <- numeric(K)
    for(k in 1:K){
      mukhats[k] <- (rnorm(1)*sig*sqrt(n) + n*deltav[k])/(R[1]*n)
    }
    #treatment effect in the control treatment
    mu0hats <- rnorm(1, 0, sig * sqrt(n))/(r0*n)
    #compute the corresponding zvalues, 
    zks <- (mukhats - mu0hats)/(sig * sqrt((R + r0)/(R * 
                                                       r0 * n)))
    #has the trial stopped for futility
    fut <- 1-(any(zks>l))
    
    #has the trial stopped for effficacy
    #turn observations for the current stage into matrix (used for checking rejections)
    Zmat = t(matrix(zks,nrow=length(zks),ncol=nrow(s1tests)))
    
    #which hypotheses can be rejected globally
    #on the left indicators of which tests have been rejected for which hypotheses
    #this checks which tests have been rejected, then iff all tests for a corresponding hypotheses are rejected
    rej = (colSums(s1tests*(rowSums(s1tests*(Zmat > u))>0)) == colSums(s1tests))
    eff = sum(rej)
    
    if((1-fut) & (1-eff)){
      cond <- which(zks > l)
      temp <- emams
      temp$K <- length(cond)
      obs <- matrix(zks[cond],nrow=1)
      calpha <- conderror(temp,obs)
      
      #simulate the remainder of the trial
      
      #all required hypothesis tests
      rtestsm <- t(combn(temp$K+Kadd,1,FUN=tabulate,nbins=temp$K+Kadd))
      for(j in 2:(temp$K+Kadd)){
        rtestsm <- rbind(rtestsm,t(combn(temp$K+Kadd,j,FUN=tabulate,nbins=temp$K+Kadd)))
      }
      rtests <- lapply(seq_len(nrow(rtestsm)), function(j) rtestsm[j,])
 
      test <- rtests[[length(rtests)]]
      #full intersection off the null hypotheses conditional test (from which to take the lower testing boundaries)
      ftestc <- condtests(test,htests,data.frame(obs),temp$K,Kadd,emams$J,alpha,emams$rMat[2,],emams$rMat[1,],power,p,p0,ushape="triangular",lshape="triangular",lfix=NULL,funloc)

      #what are the lower boundaries
      lfix = ftestc$l
      lfix = lfix[-length(lfix)]
      
      ctests <- list()

      #for all other tests     
      for(j in 1:(length(rtests)-1)){
        ctests[[j]] <- condtests(rtests[[j]],htests,data.frame(obs),temp$K,Kadd,emams$J,alpha,emams$rMat[2,],emams$rMat[1,],power,p,p0,ushape="triangular",lshape="closed",lfix=lfix,funloc)
      }
      
      ctests[[length(rtests)]] = ftestc

    #  print(list(nsim=1,rtestsm,ctests[[length(ctests)]]$K,htests[[length(htests)]],data.frame(obs),c(pv[which(as.vector(zks)>=l)],pv[(length(pv)-Kadd+1):length(pv)])))
      trials[[i]] <- simcond(nsim=1,rtestsm,ctests,htests[[length(htests)]],data.frame(obs),c(pv[which(as.vector(zks)>=l)],pv[(length(pv)-Kadd+1):length(pv)]))
      
      trials[[i]]$exss <- trials[[i]]$exss + n*K
      tempmat <- c(0,0,0,0)
      tempmat[-which(zks < l)] <- trials[[i]]$hmat
      trials[[i]]$hmat <- tempmat
    }
    else{ 
      res <- list()
      res$hmat    <- rej
      res$proprej <- eff
      res$exss 	<- n*K
      trials[[i]] <- res
    }	
    print(i)
  }
  return(trials)	
}

#a version of the above that allows parallel computing 
simfullpar <- function(pv,htests,Kadd,alpha,funloc,nsim = 1000,p = 0.75, p0 = 0.5,cl) 
{	
  #the basic principle is that we simulate the whole trial
  #and then see what has happened (easier than conditional simulation
  emams <-  htests[[length(htests)]]
  #number of trial stages
  R <- emams$r[-1,1]
  r0 <-emams$r[1,1]
  u <- emams$u[1]
  l <- emams$l[1]
  n <- emams$n
  J <- 1
  K <- emams$K
  addindex <- c((K+1):(K+Kadd))
  sig <- sd
  
  #upper testing boundaries for stage 1
  temp <- t(combn(K,1,FUN=tabulate,nbins=K))*htests[[1]]$u[1]
  #stage 1 required tests
  s1tests <- t(combn(K,1,FUN=tabulate,nbins=K))
  for(i in 2:K){
    temp <- rbind(temp,t(combn(K,i,FUN=tabulate,nbins=K))*htests[[i]]$u[1])
    s1tests <- rbind(s1tests,t(combn(K,i,FUN=tabulate,nbins=K)))
  }
  u <- temp
  rm(temp)
  
  if (is.numeric(pv)) {
    deltav <- sqrt(2) * qnorm(pv)
    sig <- 1
  }
  
  simfullloop <- function(nsim,R,r0,l,n,alpha,J,K,addindex,s1tests,u,deltav,sig,pv,htests,Kadd,funloc,p,p0,emams){
    source(funloc,local=TRUE)
    eff <- 0; fut <- 0; calpha <- 0
    #treatment means in the experimental arms
    #NOTE this is where the simulation is actually happening
    
    #need to sumulate full set of muhats
    mukhats <- numeric(K)
    for(k in 1:K){
      mukhats[k] <- (rnorm(1)*sig*sqrt(n) + n*deltav[k])/(R[1]*n)
    }
    #treatment effect in the control treatment
    mu0hats <- rnorm(1, 0, sig * sqrt(n))/(r0*n)
    #compute the corresponding zvalues, 
    zks <- (mukhats - mu0hats)/(sig * sqrt((R + r0)/(R * 
                                                       r0 * n)))
    #has the trial stopped for futility
    fut <- 1-(any(zks>l))
    
    #has the trial stopped for effficacy
    #turn observations for the current stage into matrix (used for checking rejections)
    Zmat = t(matrix(zks,nrow=length(zks),ncol=nrow(s1tests)))
    
    #which hypotheses can be rejected globally
    #on the left indicators of which tests have been rejected for which hypotheses
    #this checks which tests have been rejected, then iff all tests for a corresponding hypotheses are rejected
    rej = (colSums(s1tests*(rowSums(s1tests*(Zmat > u))>0)) == colSums(s1tests))
    eff = sum(rej)
    
    if((1-fut) & (1-eff)){
      cond <- which(zks > l)
      temp <- emams
      temp$K <- length(cond)
      obs <- matrix(zks[cond],nrow=1)
      calpha <- conderror(temp,obs)
      
      #simulate the remainder of the trial
      
      #all required hypothesis tests
      rtestsm <- t(combn(temp$K+Kadd,1,FUN=tabulate,nbins=temp$K+Kadd))
      for(j in 2:(temp$K+Kadd)){
        rtestsm <- rbind(rtestsm,t(combn(temp$K+Kadd,j,FUN=tabulate,nbins=temp$K+Kadd)))
      }
      rtests <- lapply(seq_len(nrow(rtestsm)), function(j) rtestsm[j,])
      
      test <- rtests[[length(rtests)]]
      #full intersection off the null hypotheses conditional test (from which to take the lower testing boundaries)
      ftestc <- condtests(test,htests,data.frame(obs),temp$K,Kadd,emams$J,alpha,emams$rMat[2,],emams$rMat[1,],power,p,p0,ushape="triangular",lshape="triangular",lfix=NULL,funloc)
      
      #what are the lower boundaries
      lfix = ftestc$l
      lfix = lfix[-length(lfix)]
      
      ctests <- list()
      
      test <- rtests[[1]]
      condtests(test,htests,data.frame(obs),temp$K,Kadd,emams$J,alpha,emams$rMat[2,],emams$rMat[1,],power,p,p0,ushape="triangular",lshape="closed",lfix=lfix,funloc)
      
      #for all other tests     
      for(j in 1:(length(rtests)-1)){
        ctests[[j]] <- condtests(test,htests,data.frame(obs),temp$K,Kadd,emams$J,alpha,emams$rMat[2,],emams$rMat[1,],power,p,p0,ushape="triangular",lshape="closed",lfix=lfix,funloc)
      }
      
      ctests[[length(rtests)]] = ftestc
      
      trials <- simcond(nsim=1,rtestsm,ctests,htests[[length(htests)]],data.frame(obs),c(pv[which(as.vector(zks)>=l)],pv[(length(pv)-Kadd+1):length(pv)]))
      
      trials$exss <- trials$exss + n*K
      tempmat <- c(0,0,0,0)
      tempmat[-which(zks < l)] <- trials$hmat
      trials$hmat <- tempmat
    }
    else{ 
      res <- list()
      res$hmat    <- c(rej,rep(0,Kadd))
      res$proprej <- eff
      res$exss 	<- n*K
      trials <- res
    }
    return(trials)
  }
  
  nsimlist = lapply(seq_len(nsim), function(i) i)
  
  #what are all the required components of a closed MAMS trial
  results = parLapply(cl=cl,
                      X = nsimlist,
                      fun = simfullloop,
                      R,r0,l,n,alpha,J,K,addindex,s1tests,u,deltav,sig,pv,htests,Kadd,funloc,p,p0,emams)
  
  return(results)
}

#####################################################################################################

###mamsc
##this is a slightly modified version of the mams() function from the mams package that allows
##lower bounds to be directly defined
mamsc <- function (K = 4, J = 2, alpha = 0.05, power = 0.9, r = 1:2, r0 = 1:2, 
                   p = 0.75, p0 = 0.5, delta = NULL, delta0 = NULL, sd = NULL, 
                   ushape = "obf", lshape = "fixed", ufix = NULL, 
                   lfix = 0, nstart = 1, nstop = NULL, sample.size = TRUE, N = 20, 
                   type = "normal") 
{
  require(mvtnorm)
  mesh <- function(x, d, w = 1/length(x) + x * 0) {
    n <- length(x)
    W <- X <- matrix(0, n^d, d)
    for (i in 1:d) {
      X[, i] <- x
      W[, i] <- w
      x <- rep(x, rep(n, length(x)))
      w <- rep(w, rep(n, length(w)))
    }
    w <- exp(rowSums(log(W)))
    list(X = X, w = w)
  }
  prodsum <- function(x, l, u, r, r0, r0diff, J, K, Sigma) {
    int <- prod(sapply(x, dnorm))
    L <- sqrt(r[1])/r0[1] * sqrt(r0diff[1]) * x[1] + l[1] * 
      sqrt(1 + r[1]/r0[1])
    U <- sqrt(r[1])/r0[1] * sqrt(r0diff[1]) * x[1] + u[1] * 
      sqrt(1 + r[1]/r0[1])
    insum <- pnorm(L)
    if (J > 1) {
      for (j in 2:J) {
        L[j] <- sqrt(r[j])/r0[j] * (sqrt(r0diff[1:j]) %*% 
                                      x[1:j]) + l[j] * sqrt(1 + r[j]/r0[j])
        U[j] <- sqrt(r[j])/r0[j] * (sqrt(r0diff[1:j]) %*% 
                                      x[1:j]) + u[j] * sqrt(1 + r[j]/r0[j])
        insum <- insum + pmvnorm(lower = c(L[1:j - 1], 
                                           -Inf), upper = c(U[1:j - 1], L[j]), sigma = Sigma[1:j, 
                                                                                             1:j])[1]
      }
    }
    int <- int * insum^K
    return(int)
  }
  typeI <- function(C, alpha, N, r, r0, r0diff, J, K, Sigma, 
                    ushape, lshape, lfix = NULL, ufix = NULL) {
    if (!is.function(ushape)) {
      if (ushape == "obf") {
        u <- C * sqrt(r[J]/r)
      }
      else if (ushape == "pocock") {
        u <- rep(C, J)
      }
      else if (ushape == "fixed") {
        u <- c(rep(ufix, J - 1), C)
      }
      else if (ushape == "triangular") {
        u <- C * (1 + r/r[J])/sqrt(r)
      }
    }
    else {
      u <- C * ushape(J)
    }
    if (!is.function(lshape)) {
      if (lshape == "obf") {
        l <- c(-C * sqrt(r[J]/r[1:(J - 1)]), u[J])
      }
      else if (lshape == "pocock") {
        l <- c(rep(-C, J - 1), u[J])
      }
      else if (lshape == "fixed") {
        l <- c(rep(lfix, J - 1), u[J])
      }
      else if (lshape == "closed") {
        l <- c(lfix,u[J])
        #for optimisation need to ensure u>l
        u <- (u>=l)*u + (u<l)*l
      }
      else if (lshape == "triangular") {
        if (ushape == "triangular") {
          l <- -C * (1 - 3 * r/r[J])/sqrt(r)
        }
        else {
          l <- -C * (1 - 3 * r/r[J])/sqrt(r)/(-1 * (1 - 
                                                      3)/sqrt(J))
        }
      }
    }
    else {
      l <- c(C * lshape(J)[1:(J - 1)], u[J])
    }
    mmp <- mesh((1:N - 0.5)/N * 12 - 6, J, rep(12/N, N))
    evs <- apply(mmp$X, 1, prodsum, l = l, u = u, r = r, 
                 r0 = r0, r0diff = r0diff, J = J, K = K, Sigma = Sigma)
    truealpha <- 1 - mmp$w %*% evs
    return(truealpha - alpha)
  }
  prodsum2 <- function(x, r, r0, l, u, K, delta, delta0, n, 
                       sig) {
    int <- dnorm(x)
    int <- int * pnorm(x + (delta - delta0) * sqrt(r[1] * 
                                                     n)/sig)^(K - 1) * pnorm(sqrt(r0[1]/r[1]) * (x + delta * 
                                                                                                   sqrt(r[1] * n)/sig - u[1] * sqrt(1 + r[1]/r0[1])))
    return(int)
  }
  prodsum3 <- function(x, l, u, r, r0, r0diff, J, K, delta, 
                       delta0, n, sig, Sigma, SigmaJ) {
    int <- prod(sapply(x, dnorm))
    L <- sqrt(r[1])/r0[1] * sqrt(r0diff[1]) * x[1] + l[1] * 
      sqrt(1 + r[1]/r0[1]) - delta0 * sqrt(r[1] * n)/sig
    U <- sqrt(r[1])/r0[1] * sqrt(r0diff[1]) * x[1] + u[1] * 
      sqrt(1 + r[1]/r0[1]) - delta0 * sqrt(r[1] * n)/sig
    insum <- pnorm(L)
    if (J > 2) {
      for (j in 2:(J - 1)) {
        L[j] <- sqrt(r[j])/r0[j] * (sqrt(r0diff[1:j]) %*% 
                                      x[1:j]) + l[j] * sqrt(1 + r[j]/r0[j]) - delta0 * 
          sqrt(r[j] * n)/sig
        U[j] <- sqrt(r[j])/r0[j] * (sqrt(r0diff[1:j]) %*% 
                                      x[1:j]) + u[j] * sqrt(1 + r[j]/r0[j]) - delta0 * 
          sqrt(r[j] * n)/sig
        insum <- insum + pmvnorm(lower = c(L[1:j - 1], 
                                           -Inf), upper = c(U[1:j - 1], L[j]), sigma = Sigma[1:j, 
                                                                                             1:j])[1]
      }
    }
    U[J] <- x[J] + (delta - delta0) * sqrt(r[J] * n)/sig
    insum <- insum + pmvnorm(lower = c(L, -Inf), upper = U, 
                             sigma = Sigma)[1]
    int <- int * insum^(K - 1)
    LJ <- sqrt(r[J]/(r[J] - r[1])) * (sqrt(r[1])/r0[1] * 
                                        sqrt(r0diff[1]) * x[1] + l[1] * sqrt(1 + r[1]/r0[1]) - 
                                        delta * sqrt(r[1] * n)/sig - sqrt(r[1]/r[J]) * x[J])
    UJ <- sqrt(r[J]/(r[J] - r[1])) * (sqrt(r[1])/r0[1] * 
                                        sqrt(r0diff[1]) * x[1] + u[1] * sqrt(1 + r[1]/r0[1]) - 
                                        delta * sqrt(r[1] * n)/sig - sqrt(r[1]/r[J]) * x[J])
    if (J > 2) {
      for (j in 2:(J - 1)) {
        LJ[j] <- sqrt(r[J]/(r[J] - r[j])) * (sqrt(r[j])/r0[j] * 
                                               sqrt(r0diff[1:j]) %*% x[1:j] + l[j] * sqrt(1 + 
                                                                                            r[j]/r0[j]) - delta * sqrt(r[j] * n)/sig - 
                                               sqrt(r[j]/r[J]) * x[J])
        UJ[j] <- sqrt(r[J]/(r[J] - r[j])) * (sqrt(r[j])/r0[j] * 
                                               sqrt(r0diff[1:j]) %*% x[1:j] + u[j] * sqrt(1 + 
                                                                                            r[j]/r0[j]) - delta * sqrt(r[j] * n)/sig - 
                                               sqrt(r[j]/r[J]) * x[J])
      }
    }
    int <- int * pmvnorm(lower = LJ, upper = UJ, sigma = SigmaJ)[1]
    int <- int * pnorm((r0[J]/sqrt(r[J]) * (x[J] + delta * 
                                              sqrt(r[J] * n)/sig - u[J] * sqrt(1 + r[J]/r0[J])) - 
                          (sqrt(r0diff[1:(J - 1)]) %*% x[1:(J - 1)]))/sqrt(r0diff[J]))
    return(int)
  }
  typeII <- function(n, beta, l, u, N, r, r0, r0diff, J, K, 
                     delta, delta0, sig, Sigma) {
    mmp <- mesh((1:N - 0.5)/N * 12 - 6, 1, rep(12/N, N))
    evs <- apply(mmp$X, 1, prodsum2, r = r, r0 = r0, l = l, 
                 u = u, K = K, delta = delta, delta0 = delta0, n = n, 
                 sig = sig)
    pi <- mmp$w %*% evs
    if (J > 1) {
      for (j in 2:J) {
        A <- diag(sqrt(r[j]/(r[j] - r[1:(j - 1)])), ncol = j - 
                    1)
        SigmaJ <- A %*% (Sigma[1:(j - 1), 1:(j - 1)] - 
                           Sigma[1:(j - 1), j] %*% t(Sigma[1:(j - 1), 
                                                           j])) %*% A
        mmp <- mesh((1:N - 0.5)/N * 12 - 6, j, rep(12/N, 
                                                   N))
        evs <- apply(mmp$X, 1, prodsum3, l = l, u = u, 
                     r = r, r0 = r0, r0diff = r0diff, J = j, K = K, 
                     delta = delta, delta0 = delta0, n = n, sig = sig, 
                     Sigma = Sigma[1:j, 1:j], SigmaJ = SigmaJ)
        pi <- pi + mmp$w %*% evs
      }
    }
    return(1 - beta - pi)
  }
  if (K%%1 > 0 | J%%1 > 0) {
    stop("K and J need to be integers.")
  }
  if (K < 1 | J < 1) {
    stop("The number of stages and treatments must be at least 1.")
  }
  if (N <= 3) {
    stop("Number of points for integration by quadrature to small or negative.")
  }
  if (N > 3 & N <= 10) {
    warning("Number of points for integration by quadrature is small which may result in inaccurate solutions.")
  }
  if (alpha < 0 | alpha > 1 | power < 0 | power > 1) {
    stop("Error rate or power not between 0 and 1.")
  }
  if (length(r) != length(r0)) {
    stop("Different length of allocation ratios on control and experimental treatments.")
  }
  if (length(r) != J) {
    stop("Length of allocation ratios does not match number of stages.")
  }
  if (is.numeric(p) & is.numeric(p0) & is.numeric(delta) & 
      is.numeric(delta0) & is.numeric(sd)) {
    stop("Specify the effect sizes either via (p, p0) or via (delta, delta0, sd) and set the other parameters to NULL.")
  }
  if (is.numeric(p) & is.numeric(p0)) {
    if (p < 0 | p > 1 | p0 < 0 | p0 > 1) {
      stop("Treatment effect parameter not within 0 and 1.")
    }
    if (p <= p0) {
      stop("Interesting treatment effect must be larger than uninteresting effect.")
    }
    if (p0 < 0.5) {
      warning("Uninteresting treatment effect less than 0.5 which implies that reductions in effect over placebo are interesting.")
    }
  }
  else {
    if (is.numeric(delta) & is.numeric(delta0) & is.numeric(sd)) {
      if (sd <= 0) {
        stop("Standard deviation must be positive.")
      }
    }
    else {
      stop("Specify the effect sizes either via (p, p0) or via (delta, delta0, sd).")
    }
  }
  if (is.function(ushape) & is.function(lshape)) {
    warning("You have specified your own functions for both the lower and upper boundary. Please check carefully whether the resulting boundaries are sensible.")
  }
  if (!is.function(ushape)) {
    if (!ushape %in% c("pocock", "obf", "triangular", 
                       "fixed")) {
      stop("Upper boundary does not match the available options.")
    }
    if (ushape == "fixed" & is.null(ufix)) {
      stop("ufix required when using a fixed upper boundary shape.")
    }
  }
  else {
    b <- ushape(J)
    if (!all(sort(b, decreasing = TRUE) == b)) {
      stop("Upper boundary shape is increasing.")
    }
  }
  if (!is.function(lshape)) {
    if (!lshape %in% c("pocock", "obf", "triangular", 
                       "fixed","closed")) {
      stop("Lower boundary does not match the available options.")
    }
    if (lshape == "fixed" & is.null(lfix)) {
      stop("lfix required when using a fixed lower boundary shape.")
    }
  }
  else {
    b <- lshape(J)
    if (!all(sort(b, decreasing = FALSE) == b)) {
      stop("Lower boundary shape is decreasing.")
    }
  }
  if (is.numeric(p) & is.numeric(p0)) {
    delta <- sqrt(2) * qnorm(p)
    delta0 <- sqrt(2) * qnorm(p0)
    sig <- 1
  }
  else {
    delta <- delta
    delta0 <- delta0
    p0 <- pnorm(delta0/sqrt(2 * sd^2))
    sig <- sd
  }
  h <- min(c(r, r0))
  r <- r/h
  r0 <- r0/h
  bottom <- matrix(r, J, J)
  top <- matrix(rep(r, rep(J, J)), J, J)
  top[upper.tri(top)] <- t(top)[upper.tri(top)]
  bottom[upper.tri(bottom)] <- t(bottom)[upper.tri(bottom)]
  Sigma <- sqrt(top/bottom)
  r0lag1 <- c(0, r0[1:J - 1])
  r0diff <- r0 - r0lag1
  if (!is.function(lshape)) {
    if (J == 1 & lshape == "obf") {
      lshape <- "pocock"
    }
  }
  uJ <- NULL
  try(uJ <- uniroot(typeI, c(qnorm(1 - alpha)/2, 5), alpha = alpha, 
                    N = N, r = r, r0 = r0, r0diff = r0diff, J = J, K = K, 
                    Sigma = Sigma, ushape = ushape, lshape = lshape, lfix = lfix, 
                    ufix = ufix, tol = 0.001)$root, silent = TRUE)
  
  if (is.null(uJ)) {
    stop("No boundaries can be found.")
  }
  if (!is.function(ushape)) {
    if (ushape == "obf") {
      u <- uJ * sqrt(r[J]/r)
    }
    else if (ushape == "pocock") {
      u <- rep(uJ, J)
    }
    else if (ushape == "fixed") {
      u <- c(rep(ufix, J - 1), uJ)
    }
    else if (ushape == "triangular") {
      u <- uJ * (1 + r/r[J])/sqrt(r)
    }
  }
  else {
    u <- uJ * ushape(J)
  }
  if (!is.function(lshape)) {
    if (lshape == "obf") {
      l <- c(-uJ * sqrt(r[J]/r[1:(J - 1)]), u[J])
    }
    else if (lshape == "pocock") {
      l <- c(rep(-uJ, J - 1), u[J])
    }
    else if (lshape == "fixed") {
      l <- c(rep(lfix, J - 1), u[J])
    }
    else if (lshape == "closed") {
      l <- c(lfix,u[J])
      #for optimisation need to ensure u>l
      u <- (u>=l)*u + (u<l)*l
    }
    else if (lshape == "triangular") {
      if (ushape == "triangular") {
        l <- -uJ * (1 - 3 * r/r[J])/sqrt(r)
      }
      else {
        l <- -uJ * (1 - 3 * r/r[J])/sqrt(r)/(-1 * (1 - 
                                                     3)/sqrt(J))
      }
    }
  }
  else {
    l <- c(uJ * lshape(J)[1:(J - 1)], u[J])
  }
  alpha.star <- numeric(J)
  alpha.star[1] <- typeI(u[1], alpha = 0, N = N, r = r[1], 
                         r0 = r0[1], r0diff = r0diff[1], J = 1, K = K, Sigma = Sigma, 
                         ushape = "fixed", lshape = "fixed", lfix = NULL, 
                         ufix = NULL)
  if (J > 1) {
    for (j in 2:J) {
      alpha.star[j] <- typeI(u[j], alpha = 0, N = N, r = r[1:j], 
                             r0 = r0[1:j], r0diff = r0diff[1:j], J = j, K = K, 
                             Sigma = Sigma, ushape = "fixed", lshape = "fixed", 
                             lfix = l[1:(j - 1)], ufix = u[1:(j - 1)])
    }
  }
  if (J == 1 & p0 == 0.5) {
    if (r0 > r) {
      r <- r/r0
      r0 <- r0/r0
    }
    rho <- r/(r + r0)
    corr <- matrix(rho, K, K) + diag(1 - rho, K)
    quan <- qmvnorm(1 - alpha, mean = rep(0, K), corr = corr)$quantile
    n <- ((quan + qnorm(power))/(qnorm(p) * sqrt(2)))^2 * 
      (1 + 1/r)
  }
  else {
    n <- nstart
    pow <- 0
    if (sample.size) {
      if (is.null(nstop)) {
        nx <- nstart
        po <- 0
        while (po == 0) {
          nx <- nx + 1
          po <- (typeII(nx, beta = 1 - power, l = l, 
                        u = u, N = N, r = r, r0 = r0, r0diff = r0diff, 
                        J = 1, K = K, delta = delta, delta0 = delta0, 
                        sig = sig, Sigma = Sigma) < 0)
        }
        nstop <- 3 * nx
      }
      while (pow == 0 & n <= nstop) {
        n <- n + 1
        pow <- (typeII(n, beta = 1 - power, l = l, u = u, 
                       N = N, r = r, r0 = r0, r0diff = r0diff, J = J, 
                       K = K, delta = delta, delta0 = delta0, sig = sig, 
                       Sigma = Sigma) < 0)
      }
      if ((n - 1) == nstop) {
        warning("The sample size was limited by nstop.")
      }
    }
    else {
      n <- NULL
    }
  }
  res <- NULL
  res$l <- l
  res$u <- u
  res$n <- n
  res$rMat <- rbind(r0, matrix(r, ncol = J, nrow = K, byrow = TRUE))
  res$N <- sum(ceiling(res$rMat[, J] * res$n))
  res$K <- K
  res$J <- J
  res$alpha <- alpha
  res$alpha.star <- alpha.star
  if (sample.size) {
    res$power <- power
  }
  else {
    res$power <- NA
  }
  res$type <- type
  class(res) <- "MAMS"
  return(res)
}

###mamsc.sim.out
##a modified version of mams.sim that simulates mams designs created by mamsc
mamsc.sim.out <- function (nsim = 1000, nMat = matrix(c(44, 88), nrow = 2, ncol = 5), rtestsm,
                           u = c(3.068, 2.169), l = c(0, 2.169), pv = rep(0.5, 4), deltav = NULL, 
                           sd = NULL, ptest = 1) 
{
  if (is.numeric(pv) & is.numeric(deltav) & is.numeric(sd)) {
    stop("Specify the effect sizes either via pv or via deltav and sd, and set the other parameters to NULL.")
  }
  if (is.numeric(pv)) {
    if (any(pv < 0) | any(pv > 1)) {
      stop("Treatment effect parameter not within 0 and 1.")
    }
    if (length(pv) != (ncol(nMat) - 1)) 
      stop("Length of pv is not K.")
  }
  else {
    if (is.numeric(deltav) & is.numeric(sd)) {
      if (length(deltav) != (ncol(nMat) - 1)) 
        stop("Length of deltav is not K.")
      if (sd <= 0) {
        stop("Standard deviation must be positive.")
      }
    }
    else {
      stop("Specify the effect sizes either via pv or via deltav and sd.")
    }
  }
  
  sim <- function(n, l, u, R, r0, delta, sig) {
    J <- dim(R)[1]
    K <- dim(R)[2]
    Rdiff <- R - rbind(0, R[-J, ])
    r0diff <- r0 - c(0, r0[-J])
    mukhats <- apply(matrix(rnorm(J * K), J, K) * sig * sqrt(Rdiff * 
                                                               n) + Rdiff * n * matrix(delta, nrow = J, ncol = K, 
                                                                                       byrow = TRUE), 2, cumsum)/(R * n)
    mu0hats <- cumsum(rnorm(J, 0, sig * sqrt(r0diff * n)))/(r0 * 
                                                              n)
    zks <- (mukhats - mu0hats)/(sig * sqrt((R + r0)/(R * 
                                                       r0 * n)))
    eff <- 0
    fut <- 0
    ss <- 0
    j <- 1
    hmat <- matrix(0, nrow = J, ncol = K)
    remaining <- rep(T, K)
    while ((eff == 0) && (fut == 0) && (j <= J)) {
      ss <- sum((n * Rdiff[j, ])[remaining]) + n * r0diff[j] + ss
      
      Zmat = t(matrix(zks[j,],nrow=length(zks[j,]),ncol=nrow(rtestsm)))
      
      rej = (colSums(rtestsm*(rowSums(rtestsm*(Zmat > u[[j]]))>0)) == colSums(rtestsm))*remaining
      eff = sum(rej)
      
      #    eff <- (max(zks[j, remaining]) > u[j])
      fut <- (max(zks[j, remaining]) < l[j])
      
      hmat[j,]  = rej
      
      remaining <- (zks[j, ] > l[j]) & remaining
      j <- j + 1
    }
    rej <- eff
    pow <- (remaining[1]) && eff
    if (pow) {
      pow <- (which(zks[j - 1, remaining] == max(zks[j - 
                                                       1, remaining])) == 1)
    }
    if (all(remaining == 0)) {
      rem <- 0
    }
    else {
      rem <- as.numeric(which(remaining))
    }
    return(list(rej = rej, pow = pow, ess = ss, hmat = hmat))
  }
  r0 <- nMat[, 1]/nMat[1, 1]
  if (ncol(nMat) == 2) {
    R <- t(t(nMat[, -1]/nMat[1, 1]))
  }
  else {
    R <- nMat[, -1]/nMat[1, 1]
  }
  if (!is.matrix(R) && is.vector(R)) 
    R <- t(as.matrix(nMat[, -1]/nMat[1, 1]))
  n <- nMat[1, 1]
  if (is.numeric(pv)) {
    deltas <- sqrt(2) * qnorm(pv)
    sig <- 1
  }
  else {
    deltas <- deltav
    sig <- sd
  }
  reps <- sapply(rep(n, nsim), sim, l, u, R, r0, deltas, sig)
  rej <- 0
  K <- dim(R)[2]
  hmat <- matrix(0,ncol=K,nrow=nsim)
  stage1 <- matrix(0,ncol=K,nrow=nsim)
  for (i in 1:nsim) {
    if (any(reps["hmat", i][[1]][, ptest] > 0)) {
      rej <- rej + 1
    }
    hmat[i,] <- colSums(reps["hmat", i][[1]])
    #should be genralised to the stages before addition
    stage1[i,] <- reps["hmat", i][[1]][1,]
  }
  
  res <- NULL
  res$hmat <- hmat
  res$stage1 <- stage1
  res$cont <- length(which(unlist(reps["ess", ])> n*(K+1)))/nsim
  res$l <- l
  res$u <- u
  res$n <- n
  res$K <- dim(R)[2]
  res$J <- dim(R)[1]
  res$rMat <- rbind(r0, t(R))
  res$N <- sum(res$rMat[, res$J] * res$n)
  res$nsim <- nsim
  res$typeI <- mean(unlist(reps["rej", ]))
  res$power <- mean(unlist(reps["pow", ]))
  res$ptest <- ptest
  res$prop.rej <- rej/nsim
  res$exss <- mean(unlist(reps["ess", ]))
  class(res) <- "MAMS.sim"
  return(res)
}

#####################################################################################################

###condtest
##function for computing all elements of the closed test for the conditionally updated trial
##this splits things between mamsc (for the tests that do not require conditioning) and 
##condmams (for tests that do)
##inputs and outputs are as previously described with each test of the required closed testing
##procedure defined
condtests <- function(test,htests,z,K,Kadd,J,alpha,r,r0,power,p,p0,ushape="triangular",lshape="triangular",lfix=NULL,funloc)
{
  #source in the required functions
  source(funloc,local=TRUE)
  #which hypotheses are in this test
  hyp <- which(test==1)
  os  <- nrow(z)
  
  #does this hypothesis require conditioning
  if(min(hyp)>K){
    #if not then construct as mams
    out <- mamsc(K=sum(test),J=(J-os),alpha=alpha,r=r[-c(1:os)],r0=r0[-c(1:os)],power=power,p=p,p0=p0,ushape=ushape,lshape=lshape,lfix=lfix,sample.size=FALSE) 
  }else{
    trial <- htests[sum(test[1:K])][[1]]
    
    zt <- as.data.frame(z[,which(test[1:K]==1)])
    
    calpha <- conderror(trial,zt,N = 20)
    
    out <- condmams(trial,
                    zt,
                    Kadd,
                    calpha,
                    p=p,
                    p0=p0,
                    ushape=ushape,
                    lshape=lshape,
                    lfix=lfix)
  }
  out$test <- test
  out
}

#####################################################################################################
#parallelisable versions of the coresponding functions above should the user wish to use parallel package in R
simcondp <- function(pv,nsim,rtestsm,ctests,htests,z,funloc){
  source(funloc,local=TRUE)
   simcond(nsim,rtestsm,ctests,htests[[length(htests)]],z,pv)
}
mamsc.sim.outp <- function (pv = rep(0.5, 4),nsim = 1000, nMat = matrix(c(44, 88), nrow = 2, ncol = 5), rtestsm,
                           u = c(3.068, 2.169), l = c(0, 2.169),  deltav = NULL, 
                           sd = NULL, ptest = 1,funloc) 
{
  source(funloc,local=TRUE)
  mamsc.sim.out(nsim,nMat,rtestsm,u,l,as.numeric(pv),deltav,sd,ptest) 
}