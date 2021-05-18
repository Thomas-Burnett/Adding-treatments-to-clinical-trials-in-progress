###running single stage arm adding examples
##add location of the functions file here
source()

#########################################################################################################
##choose the parameters to be used

#chooose delta = qnorm(0.95)+qnorm(0.9) gives power of 0.9 for individual hypotheses
delta = qnorm(0.95)+qnorm(0.9)
#timing of the arm adding
tau = 0.5
#define the xi as per the paper
xi1 = delta
xi2 = 0
#correlation structure
rho = 0.5
#nominal alpha
alpha = 0.05
#number of simulations (note set low here for demonstration of methods)
n = 1000

#########################################################################################################
##examine the rejection probabilities

rejprobsv(xi1,xi2,tau,rho,alpha,n)
#note the output is null while full results out output as messages

#########################################################################################################
##plotting rejection probabilities

#main title for the plot
main = expression(xi[1] ~ "=" ~ delta ~ ", " ~ xi[2] ~ "=" ~ 0)

#produce the plot
plotcondv(xi1,xi2,main,n,rho,alpha,tau)