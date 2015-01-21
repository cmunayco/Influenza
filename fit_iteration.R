##################################################################################
##################################################################################
# an R script to fit the parameters of a simple SIR model to influenza epidemic
# data by minimizing the Least Squares statistic
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

##################################################################################
# iterate a whole bunch of times and sample a different value of R0 and t0
# at each iteration.  Compare the resulting SIR model to the data, and obtain
# the least squares statistic.  Put the information
# into vectors, to enable a nice summary plot at the end
##################################################################################
vR0 = numeric(0)
vt0 = numeric(0)
vbest_leastsq_fit = numeric(0) # this will contain the model estimate of the best least squares fit

vleastsq = numeric(0) # this will containt the least-squared statistic calculated for the R0 and t0 hypotheses

amin_least = 1e10 # this will be used to check if each model is a better least squares fit than obtained from other R0 and t0 hypothesis

niter = 10000  # number of different R0,t0 hypotheses we should test
for (iter in 1:niter){

   if (iter%%100==0) cat("Doing iteration ",iter," out of ",niter,"\n")  # inform user the script is doing something, and not hung

   R0 = runif(1,1.16,1.26)             # randomly sample R0 uniformly 
   t0 = as.integer(runif(1,160,260))   # randomly sample t0 uniformly 
   #R0 = rnorm(1,1.20,0.10) # preferential sampling close to a certain value
   #t0 = rnorm(1,200,20)    # preferential sampling close to a certain value
   
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

         vR0 = append(vR0,R0)
         vt0 = append(vt0,t0)
         vleastsq = append(vleastsq,least_squares)

         if (least_squares<amin_least){
            amin_least = least_squares
            vbest_leastsq_fit = incidence
            if (length(vR0)>100){
              ######################################################################## 
              # plot the results for the least squares fit
              ######################################################################## 
              quartz(width=10, height=6, pointsize=12)
              par(mfrow=c(2,2)) # set the plotting area to one plot per page
              amax = 100 # cut off for the least squares for the plot
              plot(vR0[vleastsq<amax],vleastsq[vleastsq<amax],xlab="R0 hypothesis",ylab="Least squares goodness of fit statistic")
              points(vR0[vleastsq==min(vleastsq)],vleastsq[vleastsq==min(vleastsq)],col=4,cex=1.5) # indicate the best-fit value on the plot
              legend("bottomleft",legend=c("best fit value"),col=c(4),lwd=3,bty="n",cex=0.7)
              
              plot(vt0[vleastsq<amax],vleastsq[vleastsq<amax],xlab="t0 hypothesis",ylab="Least squares goodness of fit statistic")
              points(vt0[vleastsq==min(vleastsq)],vleastsq[vleastsq==min(vleastsq)],col=4,cex=1.5) # indicate the best-fit value on the plot
              
              plot(adat$week,adat$B,ylim=c(0,1.5*max(adat$B)),xlab="Time, in weeks relative to Jan 1, 2007",ylab="Incidence",cex=2,main="Confirmed B influenza cases\n Midwest region, 2007-2008 season")
              lines(adat$week,vbest_leastsq_fit,col=2,lwd=5) # overlay the model
              points(adat$week,adat$B,cex=2)
              
              legend("topleft",legend=c("Data","Least-squares best-fit SIR model prediction"),col=c(1,2),lwd=3,bty="n",cex=0.7)
              
              ######################################################################## 
              # print the least squares info on the plot
              ######################################################################## 
              text(47,180,paste("When R0 = ",round(vR0[which.min(vleastsq)],2),sep=""),adj=0,cex=0.6)
              text(47,160,paste(" and t0  = ",vt0[which.min(vleastsq)],",",sep=""),adj=0,cex=0.6)
              text(47,130,paste(" the least-squares \n statistic is = ",round(min(vleastsq),1),sep=""),adj=0,cex=0.6)
            }
         }
   
   
      } # end check that the length of the incidence vector matches the length of the data vector
   } # end check that the incidence vector actually has some data in it
} # end loop over iterations




########################################################################
# write the results out to a file for later reference
########################################################################
zdat = data.frame(R0=vR0,t0=vt0,leastsq=vleastsq)
write.table(zdat,"results_of_least_squares_fit_to_midwest_flu_data.txt",row.names=F)

