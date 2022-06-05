##############################################################################
#                                                                            #
# smooth                                                                     #
# ------                                                                     #
#                                                                            #
# Function for smoothing observed data and derivatives by local quadratic    #
# regression implemented by Gaussian weighting of the data points.           #
# The smoothing parameter is the standard deviation of the Gaussian weights. #
#                                                                            #
# Peter Reichert 05.09.2008 , last modification 16.12.2008                   #
#                                                                            #
##############################################################################

# Call:   smooth(x,y,sigma,newx=NA)
# -----
#
# Input:
# ------
#
# x       vector of x-coordinates of data points to be smoothed
# y       vector of y-coordinates of data points to be smoothed
#         (x and y must be of the same length)
# sigma   standard deviation of Gaussian distribution used as weights for
#         local quadratic regression
# newx    optional vector of x-coordinates at which smoothed results and
#         derivatives are to be calculated (if not specified, results
#         are provided at the same locations as there is data available)
#
# Output:
# -------
#
# List of
# x       x-coordinates at which smoothed results and derivatives are available
# y       smoothed results at the locations x
# ydot    derivatives of smoothed results at the locations x

smooth <- function(x,y,sigma,newx=NA)
{
   ind <- !is.na(y)
   x.loc <- x[ind]
   y.loc <- y[ind]

   if ( length(x.loc) < 1 )
   {
      stop("*** error in smooth: length of x is zero ***")
   }
   if ( length(x.loc) != length(y.loc) ) 
   {
      stop("*** error in smooth: length of x and y different ***")
   }
   if ( ! sigma > 0 ) 
   {
      stop("*** error in smooth: sigma is not positive ***")
   }
   if ( is.na(newx[1]) ) newx <- x.loc
   n <- length(newx)
   ysmooth <- rep(NA,n)
   ydot    <- rep(NA,n)
   for ( i in 1:n )
   {
      ind <- x.loc>newx[i]-2*sigma & x.loc<newx[i]+2*sigma
      if ( sum(ifelse(x,1,0)) > 0 )
      {
         ind <- x.loc>newx[i]-5*sigma & x.loc<newx[i]+5*sigma
         x1  <- x.loc[ind]-newx[i]
         x2  <- (x.loc[ind]-newx[i])^2
         if ( length(x1) == 1 )  # use value
         {
            ysmooth[i] <- y.loc[ind]
         }
         else
         {
            if ( length(x1) == 2 )  # use weighted mean
            {
               weights <- dnorm(x.loc[ind],mean=newx[i],sd=sigma)
               weights <- weights/sum(weights)
               ysmooth[i] <- weights[1]*y.loc[ind][1] + weights[2]*y.loc[ind][2]
            }
            else
            {
               if ( length(x1) == 3 )  # use local linear regression
               {
                  res.lm     <- lm(y.loc[ind] ~ x1,
                                   weights=dnorm(x.loc[ind],mean=newx[i],
                                                            sd=sigma))
                  ysmooth[i] <- coef(res.lm)[1]
                  ydot[i]    <- coef(res.lm)[2]
               }
               else  # use local quadratic regression
               {
                  res.lm     <- lm(y.loc[ind] ~ x1 + x2,
                                   weights=dnorm(x.loc[ind],mean=newx[i],
                                                            sd=sigma))
                  ysmooth[i] <- coef(res.lm)[1]
                  ydot[i]    <- coef(res.lm)[2]
               }
            }
         }
      }
   }
   return(list(x=newx,y=ysmooth,ydot=ydot))
}

############################################################################
