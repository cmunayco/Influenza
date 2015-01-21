##################################################################################
# An R script to solve ODE's of an SIR model 
# http://www.sherrytowers.com/sir_func.R
#
# Author: Sherry Towers
#         admin@sherrytowers.com
#
# Created: Dec 1st, 2012
# Copyright Sherry Towers, 2012
#
# This script is not guaranteed to be free of bugs and/or errors.
#
# This script can be freely used and shared as long as the author and
# copyright information in this header remain intact.
##################################################################################

require("deSolve")

##################################################################################
##################################################################################
# this is a function which, given a value of S,I and R at time t
# calculates the time derivatives of S I and R
# vparameters contains the parameters of the model, like the
# recovery period, gamma, and the transmission rate, beta
#
# this function gets passed to the ODEsolve package
##################################################################################
SIRfunc=function(t, x, vparameters){
   S = x[1]  # the value of S at time t
   I = x[2]  # the value of I at time t
   R = x[3]  # the value of R at time t
   if (I<0) I=0 # this is a cross check to ensure that we always have sensical values of I

   with(as.list(vparameters),{
      npop = S+I+R   # the population size is always S+I+R because there are no births or deaths in the model
      dS = -beta*S*I/npop            # the derivative of S wrt time
      dI = +beta*S*I/npop - gamma*I  # the derivative of I wrt time
      dR = +gamma*I                  # the derivative of R wrt time
      out = c(dS,dI,dR)
      list(out)
   })
}

##################################################################################
##################################################################################
# this is a function which, given a value of S,I and R at time t
# calculates the time derivatives of S I and R
# This model includes births and deaths (with equal rate mu)
#
# vparameters contains the parameters of the model, like the
# recovery period, gamma, and the transmission rate, beta
#
# this function gets passed to the ODEsolve package
##################################################################################
SIRfunc_with_demographics=function(t, x, vparameters){
   S = x[1]  # the value of S at time t
   I = x[2]  # the value of I at time t
   R = x[3]  # the value of R at time t
   if (I<0) I=0 # this is a cross check to ensure that we always have sensical values of I

   with(as.list(vparameters),{
      npop = S+I+R   # the population size is always S+I+R because there are no births or deaths in the model
      dS = -beta*S*I/npop - mu*S + npop*mu  # the derivative of S wrt time
      dI = +beta*S*I/npop - gamma*I - mu*I  # the derivative of I wrt time
      dR = +gamma*I - mu*R                  # the derivative of R wrt time
      out = c(dS,dI,dR)
      list(out)
   })
}

##################################################################################
##################################################################################
# This function calculations the derivatives of an SIR model with vaccination
# Rvac is the vaccinated (and now immune) compartment
# The vaccination begins at time_vaccination_begins, and ends at 
# time_vaccination_ends
#
# vparameters contains the parameters of the model, like the
# recovery period, gamma, and the transmission rate, beta
# rho is the vaccination rate (fraction of the population that
# can be vaccinated per unit time.
#
# this function gets passed to the ODEsolve package
##################################################################################
SIRfunc_with_vaccination=function(t, x, vparameters){
   S    = x[1]  # the value of S at time t
   I    = x[2]  # the value of I at time t
   R    = x[3]  # the value of R at time t
   Rvac = x[4]  # the value of Rvac at time t

   if (S<0) S=0 # this is a cross check to ensure that we always have sensical values of S
   if (I<0) I=0 # this is a cross check to ensure that we always have sensical values of I

   with(as.list(vparameters),{
      npop = S+I+R+Rvac   # the population size is always S+I+R because there are no births or deaths in the model
      dS    = -beta*S*I/npop            # the derivative of S wrt time
      dI    = +beta*S*I/npop - gamma*I  # the derivative of I wrt time
      dR    = +gamma*I                  # the derivative of R wrt time
      dRvac = 0
      if (t>=time_vaccination_begins&t<=time_vaccination_ends){
         dS    = dS - rho*S
         dRvac = +rho*S
      }
      out = c(dS,dI,dR,dRvac)
      list(out)
   })
}



