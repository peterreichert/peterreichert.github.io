# ============================================================================
# Extimate the net oxygen production rate from oxygen time series in the river
# Luteren, Switzerland
# ============================================================================
#
# Estimation is based on the methods described in
#
#   Reichert, P., Uehlinger, U. and Acuna, V.
#   Estimating stream metabolism from oxygen concentrations:
#   The Effect of Spatial Inhomogeneity
#   Submitted manuscript
#
# and
#
#   Acuna, V., Uehlinger, U. and Reichert, P.
#   Whole-stream metabolism estimation with open-channel techniques:
#   effects of the spatial heterogeneity and comparison of the single
#   and two-station techniques
#   Submitted manuscript
#
# This script is based on the function "calc_np" implemented in the file
# "calc_np.r". This function uses the smoothing function "smooth" implemented
# in the file "smooth.r".
#
# Last modification: Dec. 20, 2008


# global parameters:
# ==================

sigma         <-    0.5
nout          <-   20
sampsize.mc   <- 3000
sampsize.sens <- 1000
conf          <- c(0.05,0.95)

min.O2        <-    8.0
max.O2        <-   10.5
min.np        <-   -0.8
max.np        <-    0.0
min.np.unc    <-   -1.2
max.np.unc    <-    0.2


# load library for estimating the net oxygen production rate (np):
# ================================================================

source("calc_np.r")


# load library for sensitivity analysis:
# ======================================

if ( !require("sensitivity") ) { install.packages("sensitivity"); library(sensitivity) }


# read data and convert time to hours:
# ====================================

data        <- read.table("calc_np_luteren_o2.dat",header=TRUE,sep="\t")
datetime    <- strptime(data$datetime,format="%d/%m/%Y %H:%M")
data$t      <- as.numeric(datetime-strptime(substring(data$datetime[1],1,10),
                          format="%d/%m/%Y"))
par         <- read.table("calc_np_luteren_par.dat",header=T,sep="\t")

par$Ktau    <- par$K*par$tau
par$sd_Ktau <- sqrt(par$Ktau^2*((par$sd_K/par$K)^2+(par$sd_tau/par$tau)^2))

tout <- min(data$t)+(max(data$t)-min(data$t))*(0:(nout-1))/(nout-1)


# calculate best estimates:
# =========================

methods <- c("deriv","noderiv")

np.est <- list()
for ( reach in 1:4 )
{
   tau <- par$tau[reach]/60
   K   <- par$K[reach]*60
   d   <- par$d[reach]
   np.est[[reach]] <- list()
   for ( technique in 1:2 )
   {
      np.est[[reach]][[technique]] <- list()
      for ( version in 1:2 )
      {
         if ( technique == 1 )
         {
            res <- calc_np(t=tout,sigma=sigma,
                           t.p=data$t,p=data$AP,
                           t.dn=data$t,
                           Temp.dn=data[,paste("T",reach+1,sep="")],
                           O2.dn=data[,paste("C",reach+1,sep="")],
                           K=K,d=d,
                           method=methods[version])
         }
         else
         {
            res <- calc_np(t=tout,sigma=sigma,
                           t.p=data$t,p=data$AP,
                           t.dn=data$t,
                           Temp.dn=data[,paste("T",reach+1,sep="")],
                           O2.dn=data[,paste("C",reach+1,sep="")],
                           K=K,d=d,tau=tau,
                           t.up=data$t,
                           Temp.up=data[,paste("T",reach,sep="")],
                           O2.up=data[,paste("C",reach,sep="")],
                           method=methods[version])
         }
         np.est[[reach]][[technique]][[version]] <- res$np
      }
   }
}


# plot best estimates:
# ====================

# plot comparisons of different techniques:
# -----------------------------------------

for ( technique in 1:2 )
{
   for ( version in 1:2 )
   {
      pdf(paste("luteren_c_np_",technique,"_",version,".pdf",sep=""),
          width=6,height=6)
      par.def <- par(no.readonly=TRUE)
      par(mfrow=c(2,1),mar=c(4.5,4.5,4,2))  # c(bottom, left, top, right)

      # plot concentrations:

      plot(numeric(0),numeric(0),type="n",
           xlim=range(data$t),ylim=c(min.O2,max.O2),
           xlab="time [h]",ylab=expression(C*" "*group("[",gO[2]/m^3,"]")),
           main="Measured oxygen concentrations and calculated saturation")
      abline(v=(1:4)*24)
      Temp <- (data$T1+data$T2+data$T3+data$T4+data$T5)/5
      p    <- data$AP
      C.sat <- O2.sat(Temp,p)
      lines(data$t,C.sat,lty="dotted",lwd=2)
      for ( i in 1:5 ) lines(data$t,data[,paste("C",i  ,sep="")],lty=i,lwd=2)
      legend(x=min(data$t)+0.86*(max(data$t)-min(data$t)),
             y=min.O2+0.55*(max.O2-min.O2),
             legend=paste("station",1:5),
             lty=1:5,lwd=2,cex=0.6)

      # plot net production:

      plot(numeric(0),numeric(0),type="n",
           xlim=range(data$t),ylim=c(min.np,max.np),
           xlab="time [h]",ylab=expression(d%.%NP*" "*group("[",gO[2]/m^2/h,"]")),
           main=paste("Net Production,",technique,"station(s), method:",version))
      #abline(h=0)
      abline(v=(1:4)*24)
      for ( reach in 1:4 ) lines(tout,np.est[[reach]][[technique]][[version]],lty=reach,lwd=2)
      legend(x=min(data$t)+0.87*(max(data$t)-min(data$t)),
             y=min.np+0.45*(max.np-min.np),
             legend=paste("reach",1:4),
             lty=1:5,lwd=2,cex=0.6)

      par(par.def)
      dev.off()
   }
}

# selected best estimates plots:
# ------------------------------

pdf("luteren_meas.pdf",width=6,height=3)
par.def <- par(no.readonly=TRUE)
par(mfrow=c(1,1),mar=c(4.5,5,1,2))  # c(bottom, left, top, right)
plot(numeric(0),numeric(0),type="n",
     xlim=range(data$t),ylim=c(min.O2,max.O2),
     xlab="time [h]",ylab=expression(C*", "*C[sat]*" "*group("[",gO[2]/m^3,"]")))
abline(v=(1:4)*24)
p <- data$AP
for ( i in 1:5 ) 
{
   Temp <- data[,paste("T",i,sep="")]
   C.sat <- O2.sat(Temp,p)
   lines(data$t,C.sat,lty=i,lwd=2)
   lines(data$t,data[,paste("C",i  ,sep="")],lty=i,lwd=2)
}
legend(x=min(data$t)+0.86*(max(data$t)-min(data$t)),
       y=min.O2+0.35*(max.O2-min.O2),
       legend=paste("station",1:5),
       lty=1:5,lwd=2,cex=0.6)
par(par.def)
dev.off()

pdf("luteren_np_2sta.pdf",width=6,height=3)
par.def <- par(no.readonly=TRUE)
par(mfrow=c(1,1),mar=c(4.5,5,1,2))  # c(bottom, left, top, right)
plot(numeric(0),numeric(0),type="n",
     xlim=range(data$t),ylim=c(min.np,max.np),
     xlab="time [h]",ylab=expression(d%.%NP*" "*group("[",gO[2]/m^2/h,"]")))
#abline(h=0)
abline(v=(1:4)*24)
for ( reach in 1:4 ) lines(tout,np.est[[reach]][[2]][[2]],lty=reach,lwd=2)
legend(x=min(data$t)+0.87*(max(data$t)-min(data$t)),
       y=min.np+0.30*(max.np-min.np),
       legend=paste("reach",1:4),
       lty=1:5,lwd=2,cex=0.6)
par(par.def)
dev.off()

pdf("luteren_np_1sta.pdf",width=6,height=3)
par.def <- par(no.readonly=TRUE)
par(mfrow=c(1,1),mar=c(4.5,5,1,2))  # c(bottom, left, top, right)
plot(numeric(0),numeric(0),type="n",
     xlim=range(data$t),ylim=c(min.np,max.np),
     xlab="time [h]",ylab=expression(d%.%NP*" "*group("[",gO[2]/m^2/h,"]")))
#abline(h=0)
abline(v=(1:4)*24)
for ( reach in 1:4 ) lines(tout,np.est[[reach]][[1]][[2]],lty=reach,lwd=2)
legend(x=min(data$t)+0.87*(max(data$t)-min(data$t)),
       y=min.np+0.30*(max.np-min.np),
       legend=paste("reach",1:4),
       lty=1:5,lwd=2,cex=0.6)
par(par.def)
dev.off()


# Monte Carlo simulation:
# =======================

qlnorm.bounded <- function(p,meanlog=0,sdlog=1,lower.tail=TRUE,log.p=FALSE)
{
   if ( sum(ifelse(is.na(p),1,0)) > 0 ) stop("qlnorm.bounded: p contains NA(s)")
   p <- 0.01 + p*0.98   # cut tails
   res <- 
      qlnorm(p,meanlog=meanlog,sdlog=sdlog,lower.tail=lower.tail,log.p=log.p)
   return(res)
}

np.samp <- list()
for ( reach in 1:4 )
{
   # draw parameter sample:
   # ----------------------

   par.mean <- c(K   = par$K[reach]*60,
                 d   = par$d[reach],
                 tau = par$tau[reach]/60)
   par.sd   <- c(K   = par$sd_K[reach]*60,
                 d   = par$sd_d[reach],
                 tau = par$sd_tau[reach]/60)
   par.sdlog   <- sqrt(log(1+par.sd^2/par.mean^2))
   par.meanlog <- log(par.mean) - 0.5*par.sdlog^2

   parsamp <- as.data.frame(matrix(nrow=sampsize.mc,ncol=length(par.mean)))
   names(parsamp) <- names(par.mean)
   for ( j in 1:length(par.mean) ) 
   {
      #parsamp[,j] <- rlnorm(n=sampsize.mc,meanlog=par.meanlog[j],sdlog=par.sdlog[j])
      parsamp[,j] <- runif(n=sampsize.mc)
      parsamp[,j] <- qlnorm.bounded(parsamp[,j],
                                    meanlog=par.meanlog[j],sdlog=par.sdlog[j])
   }
   write.table(parsamp,paste("luteren_parsamp_reach",reach,".dat",sep=""),
               row.names=F,sep="\t")

   # calculate results:
   # ------------------

   np.samp[[reach]] <- list()
   for ( technique in 1:2 )
   {
      np.samp[[reach]][[technique]] <- as.data.frame(matrix(nrow=sampsize.mc,ncol=length(tout)))
      names(np.samp[[reach]][[technique]]) <- tout
      for ( i in 1:sampsize.mc )
      {
         if ( technique == 1 )
         {
            res <- calc_np(t=tout,sigma=sigma,
                           t.p=data$t,p=data$AP,
                           t.dn=data$t,
                           Temp.dn=data[,paste("T",reach+1,sep="")],
                           O2.dn=data[,paste("C",reach+1,sep="")],
                           K=parsamp$K[i],d=parsamp$d[i])
         }
         else
         {
            res <- calc_np(t=tout,sigma=sigma,
                           t.p=data$t,p=data$AP,
                           t.dn=data$t,
                           Temp.dn=data[,paste("T",reach+1,sep="")],
                           O2.dn=data[,paste("C",reach+1,sep="")],
                           K=parsamp$K[i],d=parsamp$d[i],tau=parsamp$tau[i],
                           t.up=data$t,
                           Temp.up=data[,paste("T",reach,sep="")],
                           O2.up=data[,paste("C",reach,sep="")])
         }
         np.samp[[reach]][[technique]][i,] <- res$np
      }
      write.table(np.samp[[reach]][[technique]],
                  paste("luteren_ressamp_reach",reach,"_",technique,"sta.dat",sep=""),
                  row.names=F,sep="\t")
   }
}


# plot results:
# -------------

pdf("luteren_o2_np.pdf",width=6,height=6)

par.def <- par(no.readonly=TRUE)
par(mfrow=c(4,3),xaxs="i",yaxs="i",mar=c(2.5,4,2.5,1))  # c(bottom, left, top, right)

for ( reach in 1:4 )
{
   # plot oxygen concentration:
   
   plot(numeric(0),numeric(0),type="n",
        xlim=range(data$t),ylim=c(min.O2,max.O2),
        axes=FALSE,xlab="",ylab="")
   lines(c(min(data$t),rep(max(data$t),2),rep(min(data$t),2)),
         c(rep(min.O2,2),rep(max.O2,2),min.O2),lwd=1.5)
   axis(side=1,at=c(12,24,36,48,60,72,84),
        labels=c("June 24","","June 25","","June 26","","June 27"),cex.axis=0.8)
   axis(side=2)
   abline(v=(1:4)*24)
   mtext(expression(C*", "*C[sat]*" "*group("[",gO[2]/m^3,"]")),
         side=2,line=2.5,cex=0.55)
   mtext(paste("Reach ",reach,"; dissolved oxygen",sep=""),
         line=0.5,cex=0.6)
   lines(data$t,data[,paste("C",reach+1,sep="")],lty="solid",lwd=2)
   lines(data$t,data[,paste("C",reach  ,sep="")],lty="dashed",lwd=2)
   p <- data$AP
   Temp <- 0.5*(data[,paste("T",reach,sep="")]+data[,paste("T",reach,sep="")])
   C.sat <- O2.sat(Temp,p)
   lines(data$t,C.sat,lty="dotted",lwd=2)
   
   # plot net production:
   
   for ( technique in c(2,1) )
   {
      plot(numeric(0),numeric(0),type="n",
           xlim=range(data$t),ylim=c(min.np.unc,max.np.unc),
           axes=FALSE,xlab="",ylab="")
      lines(c(min(data$t),rep(max(data$t),2),rep(min(data$t),2)),
            c(rep(min.np.unc,2),rep(max.np.unc,2),min.np.unc),lwd=1.5)
      axis(side=1,at=c(12,24,36,48,60,72,84),
           labels=c("June 24","","June 25","","June 26","","June 27"),cex.axis=0.8)
      axis(side=2)
      mtext(expression(d%.%NP*" "*group("[",gO[2]/m^2/h,"]")),
            side=2,line=2.5,cex=0.55)
      mtext(paste("Reach ",reach,"; ",c("single station technique",
                  "two stations technique")[technique],sep=""),
            line=0.5,cex=0.6)
      abline(h=0)
      abline(v=(1:4)*24)
      lines(tout,np.est[[reach]][[technique]][[2]],lwd=2)
      lower <- numeric(length(tout))
      upper <- numeric(length(tout))
      for ( i in 1:length(tout) )
      {
         q <- quantile(np.samp[[reach]][[technique]][,i],conf)
         lower[i] <- q[1]
         upper[i] <- q[2]
      }
      lines(tout,lower,lty="dashed")
      lines(tout,upper,lty="dashed")
   }
}

par(par.def)

dev.off()


# Variance-based sensitivity analysis:
# ====================================

# calculate results:
# ------------------

calc_np_forsens <- function(parsamp,tout,sigma,reach,technique,data)
{
   res <- rep(NA,nrow(parsamp))
   for ( i in 1:nrow(parsamp) )
   {
      if ( technique == 1 )
      {
         res.np <- calc_np(t=tout,sigma=sigma,
                           t.p=data$t,p=data$AP,
                           t.dn=data$t,
                           Temp.dn=data[,paste("T",reach+1,sep="")],
                           O2.dn=data[,paste("C",reach+1,sep="")],
                           K=parsamp$K[i],d=parsamp$d[i])
      }
      else
      {
         res.np <- calc_np(t=tout,sigma=sigma,
                           t.p=data$t,p=data$AP,
                           t.dn=data$t,
                           Temp.dn=data[,paste("T",reach+1,sep="")],
                           O2.dn=data[,paste("C",reach+1,sep="")],
                           K=parsamp$K[i],d=parsamp$d[i],tau=parsamp$tau[i],
                           t.up=data$t,
                           Temp.up=data[,paste("T",reach,sep="")],
                           O2.up=data[,paste("C",reach,sep="")])
      }
      res[i] <- res.np$np
   }   
   return(res)
}

# test only:
#reach <- 1; technique <- 1
#calc_np_forsens(parsamp[1:10,],tout[1],sigma,reach,technique,data)
#calc_np_forsens(data.frame(K=0.01,d=0.1,tau=0.1),tout[1],sigma,reach,technique,data)

np.sens <- list()
for ( reach in 1:4 )
{
   # specify parameter distributions:
   # --------------------------------

   par.mean <- c(K   = par$K[reach]*60,
                 d   = par$d[reach],
                 tau = par$tau[reach]/60)
   par.sd   <- c(K   = par$sd_K[reach]*60,
                 d   = par$sd_d[reach],
                 tau = par$sd_tau[reach]/60)
   par.sdlog   <- sqrt(log(1+par.sd^2/par.mean^2))
   par.meanlog <- log(par.mean) - 0.5*par.sdlog^2

   # calculate results:
   # ------------------

   np.sens[[reach]] <- list()
   for ( technique in 1:2 )
   {
      np.sens[[reach]][[technique]] <- list()
      for ( i in 1:length(tout) )
      {
         np.sens[[reach]][[technique]][[i]] <-
            fast99(model     = calc_np_forsens,
                   factors   = c("K","d","tau"),
                   n         = sampsize.sens,
                   q         = "qlnorm.bounded",
                   q.arg     = list(list(meanlog=par.meanlog["K"],sdlog=par.sdlog["K"]),
                                    list(meanlog=par.meanlog["d"],sdlog=par.sdlog["d"]),
                                    list(meanlog=par.meanlog["tau"],sdlog=par.sdlog["tau"])),
                   tout      = tout[i],
                   sigma     = sigma,
                   reach     = reach,
                   technique = technique,
                   data      = data)
      }
   }
}


# plot results:
# -------------

max.y <- 1.0

pdf("luteren_sens.pdf",width=4.5,height=6)

par.def <- par(no.readonly=TRUE)
par(mfrow=c(4,2),xaxs="i",yaxs="i",mar=c(2.5,4,2.5,1))  # c(bottom, left, top, right)

for ( reach in 1:4 )
{
   for ( technique in c(2,1) )
   {
      parnames <- names(np.sens[[reach]][[technique]][[1]]$X)
      S1 <- as.data.frame(matrix(nrow=length(tout),ncol=length(parnames)))
      names(S1) <- parnames 
      for ( i in 1:length(tout) )
      {
         S1[i,] <- np.sens[[reach]][[technique]][[i]]$D1 /
                   np.sens[[reach]][[technique]][[i]]$V
      }
      plot(numeric(0),numeric(0),type="n",
           xlim=range(data$t),ylim=c(0,max.y),
           axes=FALSE,xlab="",ylab="")
      polygon(c(min(tout),max(tout),rev(tout),tout[1]),c(0,0,rev(S1$K),0),col=grey(0.85))
      polygon(c(tout,rev(tout),tout[1]),c(S1$K,rev(S1$K+S1$d),S1$K[1]),col=grey(0.55))
      polygon(c(tout,rev(tout),tout[1]),c(S1$K+S1$d,rev(S1$K+S1$d+S1$tau),S1$K[1]+S1$d[1]),col="black")
      lines(c(min(data$t),rep(max(data$t),2),rep(min(data$t),2)),
            c(0,0,max.y,max.y,0),lwd=1.5)
      axis(side=1,at=c(12,24,36,48,60,72,84),
           labels=c("June 24","","June 25","","June 26","","June 27"),cex.axis=0.8)
      axis(side=2)
      mtext(expression("cumulative "*S[1]),
            side=2,line=2.5,cex=0.55)
      mtext(paste("Reach ",reach,"; ",c("single station technique",
                  "two stations technique")[technique],sep=""),
            line=0.5,cex=0.6)
      abline(h=0)
      abline(v=(1:4)*24)
   }
}

par(par.def)

dev.off()
