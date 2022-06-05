# ===========================================================================
# Function to estimate the net oxygen production rate from oxygen time series 
# in a river
# ===========================================================================

# Estimate the net oxygen production rate from oxygen time series in a river
# as described in 
#
#   Reichert, P., Uehlinger, U. and Acuna, V.
#   Estimating stream metabolism from oxygen concentrations:
#   The Effect of Spatial Inhomogeneity
#   Submitted manuscript
#
# Last modification: Dec. 20, 2008


# load library for smoothing of value and derivative:
# ---------------------------------------------------

source("smooth.r")


# function for calculating oxygen saturation concentration:
# ---------------------------------------------------------

# Call:    O2.sat(Temp,p=1013.25,method="")
# -----
#
# Input:
# ------
#
# Temp     (vector of) water temperature in degC
# p        (vector of) air pressure in hPa
# method   optional string either "Buehrer" or ""
#
# Output:
# -------
#
# (vector of) dissolved oxygen saturation concentration

O2.sat <- function(Temp,p=1013.25,method="")
{
   if ( method == "Buehrer" )
   {
      # Bührer, H., Computerprogramm zur Bekanntgabe aktueller Seedaten, 
      # Schweiz. Z. Hydrol. 37, 332-346, 1975:
      O2.sat <- (14.60307-0.4021469*Temp+0.00768703*Temp^2-0.0000692575*Temp^3)*p/1013.25
   }
   else
   {
      # unknown reference:
      O2.sat <- exp(7.7117-1.31403*log(Temp+45.93))*p/1013.25 
   }

   return(O2.sat)
}


# function for calculating net production from oxygen and temp. time series:
# --------------------------------------------------------------------------

# Call:    calc_np(t,sigma,
# -----            t.p,p,
#                  t.dn,Temp.dn,O2.dn,
#                  K,d=1,tau=NA,
#                  t.up=NA,Temp.up=NA,O2.up=NA,
#                  method="noderiv")
#
# Input:
# ------
#
# t        (vector of) time at which to estimate net production (time unit)
# sigma    standard deviation of normal weights for smooting (see smooth)
# t.p      (vector of) time at which air pressure is measured (time unit)
# p        (vector of) air pressure measurements (hPa)
# t.dn     vector of time points at which downstream measurements are taken
#          (time unit)
# Temp.dn  vector of downstream temperature measurements (degC)
# O2.dn    vector of downstream dissolved oxygen measurements (g/m3)
# K        reaeration coefficient of the river reach (1/time unit)
# d        river depth (optional; if not given net production is calculated
#          as gO2/m3)
# tau      transport time from upstream to downstream measurement sites
#          (time unit; only needed for two station technique)
# t.up     vector of time points at which upstream measurements are taken
#          (time unit; only needed for two station technique)
# Temp.up  vector of upstream temperature measurements (degC; only needed for 
#          two station technique)
# O2.up    vector of up stream dissolved oxygen measurements (g/m3; only 
#          needed for two station technique)
#
# Output:
# -------
#
# List of:
# t        (vector of) time at which net production is estimated (time unit)
# np       (vector of) net production estimates as g/m2/(time unit) or
#          (if depth is not provided) as g/m3/(time unit)

calc_np <- function(t,sigma,
                    t.p,p,
                    t.dn,Temp.dn,O2.dn,
                    K,d=1,tau=NA,
                    t.up=NA,Temp.up=NA,O2.up=NA,
                    method="noderiv")
{
   smooth.p    <- smooth(x=t.p,y=p,sigma=sigma,newx=t)
   if ( method != "noderiv" )
   {
      smooth.Temp.dn <- smooth(x=t.dn,y=Temp.dn,sigma=sigma,newx=t)
      smooth.O2.dn   <- smooth(x=t.dn,y=O2.dn ,sigma=sigma,newx=t)

      if ( !is.na(t.up[1])    & !is.na(O2.up[1]) & 
           !is.na(Temp.up[1]) & !is.na(tau)        )
      {
         smooth.Temp.up <- smooth(x=t.up,y=Temp.up,sigma=sigma,newx=t)
         smooth.O2.up   <- smooth(x=t.up,y=O2.up,sigma=sigma,newx=t-tau)
         O2.sat         <- O2.sat(0.5*(smooth.Temp.up$y+smooth.Temp.dn$y),
                                  smooth.p$y)
         a              <- exp(-K*tau)
         b              <- 1/(1-a)
         c              <- 1 + (K*tau-1)*a
         np             <- d*(K*((smooth.O2.dn$y-a*smooth.O2.up$y)*b-O2.sat)
                              + (smooth.O2.dn$ydot-a*smooth.O2.up$ydot)*b*b*c)
      }
      else
      {
         O2.sat         <- O2.sat(smooth.Temp.dn$y,smooth.p$y)
         np             <- d*(K*(smooth.O2.dn$y-O2.sat) + smooth.O2.dn$ydot)
      }
   }
   else
   {
      if ( !is.na(t.up[1])    & !is.na(O2.up[1]) & 
           !is.na(Temp.up[1]) & !is.na(tau)        )
      {
         t.dn.eval   <- t + 1/K - tau*exp(-K*tau)/(1-exp(-K*tau))
         smooth.Temp.dn <- smooth(x=t.dn,y=Temp.dn,sigma=sigma,newx=t)
         smooth.Temp.up <- smooth(x=t.up,y=Temp.up,sigma=sigma,newx=t)
         smooth.O2.dn   <- smooth(x=t.dn,y=O2.dn,sigma=sigma,newx=t.dn.eval)
         smooth.O2.up   <- smooth(x=t.up,y=O2.up,sigma=sigma,newx=t.dn.eval-tau)
         a              <- exp(-K*tau)
         b              <- 1/(1-a)
         O2.sat.val     <- O2.sat(0.5*(smooth.Temp.up$y+smooth.Temp.dn$y),
                                  smooth.p$y)
         np             <- d*K*((smooth.O2.dn$y-a*smooth.O2.up$y)*b-O2.sat.val)
      }
      else
      {
         t.dn.eval      <- t + 1/K
         smooth.Temp.dn <- smooth(x=t.dn,y=Temp.dn,sigma=sigma,newx=t)
         smooth.O2.dn   <- smooth(x=t.dn,y=O2.dn,sigma=sigma,newx=t.dn.eval)
         O2.sat.val     <- O2.sat(smooth.Temp.dn$y,smooth.p$y)
         np             <- d*K*(smooth.O2.dn$y-O2.sat.val)
      }
   }
   return(list(t=t,np=np))
}

