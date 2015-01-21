##################################################################################
##################################################################################
# an R script to fit the parameters of a simple SIR model to influenza epidemic
# data
#
# Author:  Sherry Towers
#          admin@sherrytowers.com
# Created: Dec 6, 2012
#
# Copyright Sherry Towers, 2012
#
# This script is not guaranteed to be free of bugs and/or errors
#
# This script can be freely used and shared as long as the author and
# copyright information in this header remain intact.
##################################################################################
require("chron")
require("deSolve")
par(pch=20)  # the dot style for the plots
setwd("/Users/cesarmunayco/Documents/R codes/flu-scripts") 
source("sir_func.R")

##################################################################################
# read in 2007 influenza surveillance data (confirmed cases) for
# the Midwest (CDC region 5) from
# http://www.cdc.gov/flu/weekly/regions2007-2008/datafinalHHS/whoregX.htm
# where X=5 is midwest
# X=1 is northeast
# X=2 is NY and NJ
# X=3 are eastern seabord states like PA, DE, etc
#
# the weeks are number of weeks relative to Jan 1, 2007
# week 1 in 2007 ended Jan 6, 2007
##################################################################################
adat = read.table("midwest_influenza_2007_to_2008.dat",header=T,sep=",")
cat("The data file contains: ",names(adat),"\n")
adat=subset(adat,week>=47)  # make sure we are far enough into the season that there is at least one case per week

##################################################################################
# The CDC data is weekly, with the date of each point corresponding to the
# date of the end of the week over which the data were collected.
# Let's convert these dates to time in days, relative to Jan 1, 2007
# We will be using this vector of dates, vday, to obtain the model estimates
# of the incidence at that time.
##################################################################################
adat$julian_day = julian(1,6,2007)+(adat$week-1)*7-julian(1,1,2007)
vday = append(min(adat$julian_day)-7,adat$julian_day)

##################################################################################
# set up the model parameters
##################################################################################
npop = 52000000 # this is approximately the population of IL IN MI MN OH WI (CDC region 5)
I_0 = 1      # put one infected person in the population
S_0 = npop-I_0
R_0 = 0      # initially no one has recovered yet

gamma = 1/3  # recovery period of influenza in days^{-1}
#R0 = 1.24    # the R0 hypothesis
#t0 = 233     # the time-of-introduction hypothesis, measured in days from Jan 1, 2007
R0 = 1.35    # the R0 hypothesis
t0 = 310     # the time-of-introduction hypothesis, measured in days from Jan 1, 2007

vt = seq(t0,550,1)  # let's determine the values of S,I and R at times in vt
beta  = R0*gamma
   
##################################################################################
# set up the parameter vector and initial values, then solve the SIR
# system of equations numerically with lsoda in the deSolve package.  Put
# the results in sirmodel
##################################################################################
vparameters = c(gamma=gamma,beta=beta)
inits = c(S=S_0,I=I_0,R=R_0)
   
sirmodel = as.data.frame(lsoda(inits, vt, SIRfunc, vparameters)) # this numerically solves the SIR model

##################################################################################
# the B influenza data is incidence, not prevalence (actually, the B influenza
# data is incidence times the fraction that the CDC actually confirms.)
#
# The incidence over some time step is the difference in S over that time
# step (in a model with no births or deaths, immigration or emigration, or 
# recurring susceptibility).
#
# The time step are derived from the dates at the end of the week of each
# data point (vday)
#
# sirmodel$time%in%vday returns the indices of the elements of simodel$time that
# are also found in vday
##################################################################################
susceptible = sirmodel$S[sirmodel$time%in%vday]  
incidence = -diff(susceptible)

##################################################################################
# from the model estimate of the incidence and the data, we can 
# estimate the fraction of cases that are confirmed
##################################################################################
frac_confirmed = sum(adat$B)/sum(incidence)
cat("The estimated fraction of confirmed cases among all infections is = ",frac_confirmed,"\n")

incidence = incidence*frac_confirmed # normalize the model prediction so area under curve
                                     # equals the sum of the data incidence
   
##################################################################################
# now let's overlay the model on the data, and calculate the least-squares
# statistic that compares the data to this model calculated
# under a particular hypothesis of R0 and t0
##################################################################################
if (length(incidence)>1&!is.na(sum(incidence))){
   if (min(incidence)>0&length(incidence)==length(adat$B)){
      ######################################################################## 
      # calculate the least squares statistic
      ######################################################################## 
      least_squares = sum((incidence-adat$B)^2)/mean(incidence)

      ######################################################################## 
      # plot the results
      ######################################################################## 
      par(mfrow=c(1,1)) # set the plotting area to one plot per page
      quartz(width=10, height=6, pointsize=10)
      plot(adat$week,adat$B,ylim=c(0,1.5*max(adat$B)),xlab="Time, in weeks relative to Jan 1, 2007",ylab="Incidence",cex=2,main="Confirmed B influenza cases, Midwest region, 2007-2008 season")
      lines(adat$week,incidence,col=2,lwd=5) # overlay the model

      ######################################################################## 
      # show the distance between the model and data for each point
      ######################################################################## 
      for (iind in 1:length(adat$week)){
         arrows(adat$week[iind],incidence[iind],adat$week[iind],adat$B[iind],code=3,lwd=1,length=0.10,col=4,angle=20)
      }
      legend("topleft",legend=c("Data","SIR model prediction","Point-by-point distance between data and model"),col=c(1,2,4),lwd=3,bty="n")
      ######################################################################## 
      # print the least squares info on the plot
      ######################################################################## 
      #text(47,180,paste("When R0 = ",R0,sep=""),adj=0)
      #text(47,160,paste(" and t0  = ",t0,",",sep=""),adj=0)
      text(47,260,paste("Model initial conditions and parameters:"),adj=0)
      text(47,240,paste(" population = ",npop,sep=""),adj=0)
      text(47,220,paste(" I0 = ",I_0," and rest of population susceptible",sep=""),adj=0)
      text(47,200,paste(" 1/gamma  = ",1/gamma," days",sep=""),adj=0)
      text(47,180,paste(" t0  = week ",round(t0/7),sep=""),adj=0)
      text(47,160,paste(" R0 = ",R0,sep=""),adj=0)
      text(47,120,paste(" The least-squares \n statistic is = ",round(least_squares,1),sep=""),adj=0)

   } # end check that the length of the incidence vector matches the length of the data vector
} # end check that the incidence vector actually has some data in it




