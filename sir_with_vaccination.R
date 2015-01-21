##################################################################################
# An R script to solve ODE's of a Susceptible Infected Recovered (SIR) model 
# with vaccination that begins at some time time_vaccination_begins
# and ends at time_vaccination_ends
# During the vaccination period, the vaccination rate is rho
# (the fraction of the population that can be vaccinated per day)
#
# http://www.sherrytowers.com/sir_with_vaccination.R
#
# Author: Sherry Towers
#         admin@sherrytowers.com
# Created: Dec 1st, 2012
#
# Copyright Sherry Towers, 2012
#
# This script is not guaranteed to be free of bugs and/or errors.
#
# This script can be freely used and shared as long as the author and
# copyright information in this header remain intact.
#
##################################################################################
setwd("/Users/cesarmunayco/Documents/R codes/flu-scripts") 
source("sir_func.R")  # this file contains a function SIRfunc_vaccination that calculates
                      # the derivatives of S I R and Rvac wrt to time

##################################################################################
##################################################################################
# Let's set up some initial conditions at time t=0
##################################################################################
npop = 10000000
I_0 = 1       # put one infected person in the population
S_0 = npop-I_0
R_0 = 0       # initially no one has recovered yet
Rvac_0 = 0    # initially no one has been vaccinated

##################################################################################
# now the parameters of the model.  Note that in some posts on sherrytowers.com
# I refer to the recovery rate as k (here it is gamma), and the transmission
# rate as b (here is is beta).
##################################################################################
vt = seq(0,150,1)  # let's determine the values of S,I and R at times in vt
gamma = 1/3        # approximate average recovery period of influenza in days^{-1} 
R0    = 1.5        # R0 of a hypothetical strain of pandemic influenza 
beta  = R0*gamma   # transmission rate of influenza is estimated by solving R0=beta/gamma for beta

##################################################################################
# if iter = 1, do not include vaccination (set the vaccination rate to 0)
# if iter = 2, vaccinate with rho = 0.005
# if iter = 3, vaccinate with rho = 0.010
##################################################################################

par(mfrow=c(1,1)) # ensures one plot per page
quartz(width=10, height=6, pointsize=10)
for (iter in 1:3){
  time_vaccination_begins = 70   # vaccination begins on this day
  time_vaccination_ends   = 100  # vaccination ends this day   
  if (iter==1) rho = 0.000 # the fraction of the population that can be vaccinated per day
  if (iter==2) rho = 0.005 # the fraction of the population that can be vaccinated per day
  if (iter==3) rho = 0.010 # the fraction of the population that can be vaccinated per day
  
  vparameters = c(gamma=gamma,beta=beta,time_vaccination_begins=time_vaccination_begins,time_vaccination_ends=time_vaccination_ends,rho=rho)
  inits = c(S=S_0,I=I_0,R=R_0,Rvac=Rvac_0)
  
  sirmodel = as.data.frame(lsoda(inits, vt, SIRfunc_with_vaccination, vparameters))
  cat("The item names in the sirmodel object are:",names(sirmodel),"\n")

  ##################################################################################
  # now let's plot the results
  # the list of R color indices can be found at 
  # http://research.stowers-institute.org/efg/R/Color/Chart/
  ##################################################################################
  cat("The vaccination fraction is ",rho,"\n")
  cat("The total fraction infected at the end of the epidemic is ",max(sirmodel$R)/npop,"\n")

  if (iter==1){
     plot(sirmodel$time                         # abscissa
         ,sirmodel$I/npop                       # ordinate
         ,type="l"                              # plot a line (not points)
         ,xlab="Time, in days"                  # label for the x axis
         ,ylab="Fraction infected (prevalence)" # label for the y axis
         ,ylim=c(0,1.1*max(sirmodel$I/npop))    # set the max on the y axis to 1.1 times the maximum value being plotted
         ,lwd=4                                 # set the line width to be wide
         ,col=iter                              # col 1 is black, 2 is red, 4 is blue, 3 is green
         ,main=paste("Vaccination between days ",time_vaccination_begins," and ",time_vaccination_ends,sep=""))
  }else{
     lines(sirmodel$time,sirmodel$I/npop,lwd=4,col=iter)
  }
} # end loop over the vaccination rate assumptions

legend("topleft",legend=c("no vaccination","vaccination rate=0.005","vaccination rate=0.010"),bty="n",lwd=4,col=c(1,2,3))





