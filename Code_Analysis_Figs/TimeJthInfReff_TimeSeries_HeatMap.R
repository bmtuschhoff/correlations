# read in all files for individual-based SIR model for each parameter combination for the level of heterogeneity in susceptibility (cv_s), heterogeneity in transmission (cv_t), and correlation (cor)
# determine mean, median, 95% CI for R-effective and time to the jth infection given a major outbreak
# plot time series results

library(dplyr)
library(ggplot2)
library(metR)
library(ggpubr)
library(egg)
library(cowplot)

# number of simulations
sims <- 200

# set working directory to get data from
setwd("C:/Users/bmt5507/Documents/Cor_HiST_Copula_I1_R0.8")

# function to read in datasets that are tab delimited
read_sav <- function(fileName) {read.table(fileName, header=T, sep="\t")}


# function to get time to jth infection adjusted to longest time series
get_timejthinf <- function(epi_sim, cv_s, max_length)
{
  # determine whether there was a major outbreak
  # set threshold for major outbreak based on level of HiS because more HiS = smaller outbreak
  if(cv_s == 0 || cv_s == 0.5)
    threshold <- 200
  if(cv_s == 1)
    threshold <- 100
  if(cv_s == 3)
    threshold <- 50
  
  # check if final epidemic size (FES) is larger than threshold for major epidemic
  FES <- epi_sim$I[nrow(epi_sim)] + epi_sim$R[nrow(epi_sim)]
  
  if(FES > threshold)  # major epidemic, set jth inf for each time step
  {
    i <- 1
    timejthinf <- rep(Inf, max_length)
    for(inf in (1:max(epi_sim$inf_num)))
    {
      jinf <- which(epi_sim$inf_num == inf)
      index_time <- jinf[1]
      timejthinf[i] <- epi_sim[index_time,]$time
      i <- i + 1
    } 
    
    return(timejthinf)
  }
  else  # no major epidemic, return all NAs
    return(rep(NA, max_length))
}



# function to get R effective throughout epidemic adjusted to longest time series
get_Reff <- function(epi_sim, cv_s, max_length)
{
  # determine whether there was a major outbreak
  # set threshold for major outbreak based on level of HiS because more HiS = smaller outbreak
  if(cv_s == 0 || cv_s == 0.5)
    threshold <- 200
  if(cv_s == 1)
    threshold <- 100
  if(cv_s == 3)
    threshold <- 50
  
  # check if final epidemic size (FES) is larger than threshold for major epidemic
  FES <- epi_sim$I[nrow(epi_sim)] + epi_sim$R[nrow(epi_sim)]
  
  if(FES > threshold)  # major epidemic, find Reff at each time step
  {
    i <- 1
    Reff <- rep(NA, max_length)
    for(time in seq(0, max(epi_sim$time)+0.1, 0.1))
    {
      lessthan_time <- which(epi_sim$time <= time)
      index_time <- lessthan_time[length(lessthan_time)]
      Reff[i] <- epi_sim[index_time,]$Reff
      i <- i + 1
    }  
    
    return(Reff)
  }
  else  # no major epidemic, return all NAs
    return(rep(NA, max_length))
}


# levels of heterogeneity in transmission, susceptibility, and correlation
cv_ts <- c(0.5, 1, 3)
cv_ss <- c(0.5, 1, 3)
cors <- c(-1,-0.5,0,0.5,1)

# make empty lists to store all data  
all_timestojthinf <- vector("list", length=length(cv_ts)*length(cv_ss)*length(cors))
all_Reffs <- vector("list", length=length(cv_ts)*length(cv_ss)*length(cors))


## get data for each level of HiS, HiT, and correlation for a specific R0
combo <- 1
for(cv_s in cv_ss)
{
  for(cv_t in cv_ts)
  {
    for(corr in cors)
    {
      # list of file names
      filePattern <- paste("results_his",cv_s,"hit",cv_t,"cor",corr,"g0.1p1000R00.8initI1s*", sep="")
      files <- list.files(pattern=filePattern)
      
      # make a list that contains all files
      read.all <- lapply(files, read_sav)
      
      # determine length of longest time series based on time
      max_time <- max(unlist(lapply(read.all, function(x) {x$time[nrow(x)]})))
      max_inf <- max(unlist(lapply(read.all, function(x) {x$inf_num[nrow(x)]})))
      if(max_inf < 50)  # make sure can use for time to 50th inf
        max_inf <- 50
      max_length_time <- length(seq(0,max_time+0.1,0.1))

      ## time to jth infection
      # make data frame containing time to jth infection for each time series of the sims
      temp <- data.frame(sapply(read.all, get_timejthinf, cv_s=cv_s, max_length=max_inf))

      # make sure all numeric
      temp <- data.frame(sapply(temp, as.numeric))

      # find mean, 95% CIs for each time point
      timejthinf_mean <- apply(temp, 1, mean, na.rm=T)
      timejthinf_median <- apply(temp, 1, median, na.rm=T)
      timejthinf_low <- apply(temp, 1, function(x){quantile(x, probs=0.025, na.rm=T)})
      timejthinf_high <- apply(temp, 1, function(x){quantile(x, probs=0.975, na.rm=T)})

      # make data frame containing time, mean, and 95% CIs
      element_name <- paste(cv_s,cv_t,corr, sep="_")
      timejthinf <- data.frame(jinf=seq(1,max_inf,1),
                               avg=timejthinf_mean,
                               median=timejthinf_median,
                               low=timejthinf_low,
                               high=timejthinf_high)

      # add data frame to list of data frames for each param combo (cv_s, cv_t, cor) and name
      all_timestojthinf[[combo]] <- timejthinf
      names(all_timestojthinf)[[combo]] <- element_name
      
      
      ## R effective
      # make data frame containing Reff for each time series of the sims
      temp <- data.frame(sapply(read.all, get_Reff, cv_s=cv_s, max_length=max_length_time))
      
      # make sure all numeric
      temp <- data.frame(sapply(temp, as.numeric))
      
      # find mean, 95% CIs for each time point
      Reff_mean <- apply(temp, 1, mean, na.rm=T)
      Reff_median <- apply(temp, 1, median, na.rm=T)
      Reff_low <- apply(temp, 1, function(x){quantile(x, probs=0.025, na.rm=T)})
      Reff_high <- apply(temp, 1, function(x){quantile(x, probs=0.975, na.rm=T)})
      
      # make data frame containing time, mean, and 95% CIs
      element_name <- paste(cv_s,cv_t,corr, sep="_")
      Reff <- data.frame(time=seq(0,max_time+0.1,0.1),
                         avg=Reff_mean,
                         median=Reff_median,
                         low=Reff_low,
                         high=Reff_high)
      
      # add data frame to list of data frames for each param combo (cv_s, cv_t, cor) and name
      all_Reffs[[combo]] <- Reff
      names(all_Reffs)[[combo]] <- element_name
      
      combo <- combo + 1
    }
  }
}


########################################################################################
## separately get data for homogeneous case (HiS=0, HiT=0)
# list of file names
filePatterns <- c()
filePatterns[1] <- paste("results_his0hit0cor0g0.1p1000R00.8initI1s*", sep="")
filePatterns[2] <- paste("results_his0.5hit0cor0g0.1p1000R00.8initI1s*", sep="")
filePatterns[3] <- paste("results_his1hit0cor0g0.1p1000R00.8initI1s*", sep="")
filePatterns[4] <- paste("results_his3hit0cor0g0.1p1000R00.8initI1s*", sep="")
filePatterns[5] <- paste("results_his0hit0.5cor0g0.1p1000R00.8initI1s*", sep="")
filePatterns[6] <- paste("results_his0hit1cor0g0.1p1000R00.8initI1s*", sep="")
filePatterns[7] <- paste("results_his0hit3cor0g0.1p1000R00.8initI1s*", sep="")

# levels of heterogeneity in transmission and susceptibility
cv_ss <- c(0,0.5,1,3,0,0,0)
cv_ts <- c(0,0,0,0,0.5,1,3)

# get data for each level of HiS and HiT
for(i in 1:7)
{
  files <- list.files(pattern=filePatterns[i])

  # make a list that contains all files
  read.all <- lapply(files, read_sav)

  # determine length of longest time series based on time
  max_time <- max(unlist(lapply(read.all, function(x) {x$time[nrow(x)]})))
  max_inf <- max(unlist(lapply(read.all, function(x) {x$inf_num[nrow(x)]})))
  if(max_inf < 50)
    max_inf <- 50
  max_length_time <- length(seq(0,max_time+0.1,0.1))

  ## time to jth infection
  # make data frame containing time to jth infection for each time series of the sims
  temp <- data.frame(sapply(read.all, get_timejthinf, cv_s=cv_ss[i], max_length=max_inf))

  # make sure all numeric
  temp <- data.frame(sapply(temp, as.numeric))

  # find mean, 95% CIs for each time point
  timejthinf_mean <- apply(temp, 1, mean, na.rm=T)
  timejthinf_median <- apply(temp, 1, median, na.rm=T)
  timejthinf_low <- apply(temp, 1, function(x){quantile(x, probs=0.025, na.rm=T)})
  timejthinf_high <- apply(temp, 1, function(x){quantile(x, probs=0.975, na.rm=T)})

  # make data frame containing time, mean, and 95% CIs
  element_name <- paste(cv_ss[i],cv_ts[i], sep="_")
  timejthinf <- data.frame(jinf=seq(1,max_inf,1),
                           avg=timejthinf_mean,
                           median=timejthinf_median,
                           low=timejthinf_low,
                           high=timejthinf_high)

  # add data frame to list of data frames for each param combo (cv_s, cv_t, cor) and name
  all_timestojthinf[[combo]] <- timejthinf
  names(all_timestojthinf)[[combo]] <- element_name


  ## R effective
  # make data frame containing Reff for each time series of the sims
  temp <- data.frame(sapply(read.all, get_Reff, cv_s=cv_ss[i], max_length=max_length_time))

  # make sure all numeric
  temp <- data.frame(sapply(temp, as.numeric))

  # find mean, 95% CIs for each time point
  Reff_mean <- apply(temp, 1, mean, na.rm=T)
  Reff_median <- apply(temp, 1, median, na.rm=T)
  Reff_low <- apply(temp, 1, function(x){quantile(x, probs=0.025, na.rm=T)})
  Reff_high <- apply(temp, 1, function(x){quantile(x, probs=0.975, na.rm=T)})

  # make data frame containing time, mean, and 95% CIs
  element_name <- paste(cv_ss[i],cv_ts[i], sep="_")
  Reff <- data.frame(time=seq(0,max_time+0.1,0.1),
                     avg=Reff_mean,
                     median=Reff_median,
                     low=Reff_low,
                     high=Reff_high)

  # add data frame to list of data frames for each param combo (cv_s, cv_t, cor) and name
  all_Reffs[[combo]] <- Reff
  names(all_Reffs)[[combo]] <- element_name
  
  combo <- combo + 1
}



################################################################################
############################# plot results #####################################
################################################################################

###################### plot time to jth infection #################################
layout(mat=matrix(c(1,2,3,
                    4,5,6,
                    7,8,9,
                    10,10,10), nrow=4, ncol=3, byrow=T),heights = c(0.4,0.4,0.4,0.05))


par(mai = c(0.55, 0.6, 0.1, 0.075))

# HiT=3
plot(all_timestojthinf$'0.5_3_-1'$median, all_timestojthinf$'0.5_3_-1'$jinf, type="l", col="red", lty="solid", 
     xlab="", ylab="Number of infections", xlim=c(0,200), ylim=c(0,1000), cex.axis=1.3, cex.lab=1.5, lwd=2)
#polygon(c(rev(all_timestojthinf$'0.5_3_-1'$median), all_timestojthinf$'0.5_3_-1'$median), c(rev(all_timestojthinf$'0.5_3_-1'$low), all_timestojthinf$'0.5_3_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0.5_3_0'$median, all_timestojthinf$'0.5_3_0'$jinf, type="l", col="yellow3", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'0.5_3_0'$median), all_timestojthinf$'0.5_3_0'$median), c(rev(all_timestojthinf$'0.5_3_0'$low), all_timestojthinf$'0.5_3_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0.5_3_1'$median, all_timestojthinf$'0.5_3_1'$jinf, type="l", col="blue", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'0.5_3_1'$median), all_timestojthinf$'0.5_3_1'$median), c(rev(all_timestojthinf$'0.5_3_1'$low), all_timestojthinf$'0.5_3_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0_0'$median, all_timestojthinf$'0_0'$jinf, type="l", col="black", lwd=2)
lines(all_timestojthinf$'0_3'$median, all_timestojthinf$'0_3'$jinf, type="l", col="orange", lwd=2)
lines(all_timestojthinf$'0.5_0'$median, all_timestojthinf$'0.5_0'$jinf, type="l", col="orchid", lwd=2)

plot(all_timestojthinf$'1_3_-1'$median, all_timestojthinf$'1_3_-1'$jinf, type="l", col="red", lty="solid", 
     xlab="", ylab="", xlim=c(0,200), ylim=c(0,1000), cex.axis=1.3, cex.lab=1.5, lwd=2)
#polygon(c(rev(all_timestojthinf$'1_3_-1'$median), all_timestojthinf$'1_3_-1'$median), c(rev(all_timestojthinf$'1_3_-1'$low), all_timestojthinf$'1_3_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'1_3_0'$median, all_timestojthinf$'1_3_0'$jinf, type="l", col="yellow3", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'1_3_0'$median), all_timestojthinf$'1_3_0'$median), c(rev(all_timestojthinf$'1_3_0'$low), all_timestojthinf$'1_3_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'1_3_1'$median, all_timestojthinf$'1_3_1'$jinf, type="l", col="blue", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'1_3_1'$median), all_timestojthinf$'1_3_1'$median), c(rev(all_timestojthinf$'1_3_1'$low), all_timestojthinf$'1_3_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0_0'$median, all_timestojthinf$'0_0'$jinf, type="l", col="black", lwd=2)
lines(all_timestojthinf$'0_3'$median, all_timestojthinf$'0_3'$jinf, type="l", col="orange", lwd=2)
lines(all_timestojthinf$'1_0'$median, all_timestojthinf$'1_0'$jinf, type="l", col="orchid", lwd=2)

plot(all_timestojthinf$'3_3_-1'$median, all_timestojthinf$'3_3_-1'$jinf, type="l", col="red", 
     xlab="Time", ylab="", xlim=c(0,200), ylim=c(0,1000), cex.axis=1.3, cex.lab=1.5, lwd=2)
#polygon(c(rev(all_timestojthinf$'3_3_-1'$median), all_timestojthinf$'3_3_-1'$median), c(rev(all_timestojthinf$'3_3_-1'$low), all_timestojthinf$'3_3_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'3_3_0'$median, all_timestojthinf$'3_3_0'$jinf, type="l", col="yellow3", lwd=2)
#polygon(c(rev(all_timestojthinf$'3_3_0'$median), all_timestojthinf$'3_3_0'$median), c(rev(all_timestojthinf$'3_3_0'$low), all_timestojthinf$'3_3_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'3_3_1'$median, all_timestojthinf$'3_3_1'$jinf, type="l", col="blue", lwd=2)
#polygon(c(rev(all_timestojthinf$'3_3_1'$median), all_timestojthinf$'3_3_1'$median), c(rev(all_timestojthinf$'3_3_1'$low), all_timestojthinf$'3_3_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0_0'$median, all_timestojthinf$'0_0'$jinf, type="l", col="black", lwd=2)
lines(all_timestojthinf$'0_3'$median, all_timestojthinf$'0_3'$jinf, type="l", col="orange", lwd=2)
lines(all_timestojthinf$'3_0'$median, all_timestojthinf$'3_0'$jinf, type="l", col="orchid", lwd=2)


# HiT=1
plot(all_timestojthinf$'0.5_1_-1'$median, all_timestojthinf$'0.5_1_-1'$jinf, type="l", col="red", lty="solid", 
     xlab="", ylab="Number of infections", xlim=c(0,200), ylim=c(0,1000), cex.axis=1.3, cex.lab=1.5, lwd=2)
#polygon(c(rev(all_timestojthinf$'0.5_1_-1'$median), all_timestojthinf$'0.5_1_-1'$median), c(rev(all_timestojthinf$'0.5_1_-1'$low), all_timestojthinf$'0.5_1_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0.5_1_0'$median, all_timestojthinf$'0.5_1_0'$jinf, type="l", col="yellow3", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'0.5_1_0'$median), all_timestojthinf$'0.5_1_0'$median), c(rev(all_timestojthinf$'0.5_1_0'$low), all_timestojthinf$'0.5_1_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0.5_1_1'$median, all_timestojthinf$'0.5_1_1'$jinf, type="l", col="blue", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'0.5_1_1'$median), all_timestojthinf$'0.5_1_1'$median), c(rev(all_timestojthinf$'0.5_1_1'$low), all_timestojthinf$'0.5_1_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0_0'$median, all_timestojthinf$'0_0'$jinf, type="l", col="black", lwd=2)
lines(all_timestojthinf$'0_1'$median, all_timestojthinf$'0_1'$jinf, type="l", col="orange", lwd=2)
lines(all_timestojthinf$'0.5_0'$median, all_timestojthinf$'0.5_0'$jinf, type="l", col="orchid", lwd=2)

plot(all_timestojthinf$'1_1_-1'$median, all_timestojthinf$'1_1_-1'$jinf, type="l", col="red", lty="solid", 
     xlab="", ylab="", xlim=c(0,200), ylim=c(0,1000), cex.axis=1.3, cex.lab=1.5, lwd=2)
#polygon(c(rev(all_timestojthinf$'1_1_-1'$median), all_timestojthinf$'1_1_-1'$median), c(rev(all_timestojthinf$'1_1_-1'$low), all_timestojthinf$'1_1_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
#lines(all_timestojthinf$'1_1_-0.5'$median, all_timestojthinf$'1_1_-0.5'$jinf, type="l", col="red", lty="solid", lwd=2)
lines(all_timestojthinf$'1_1_0'$median, all_timestojthinf$'1_1_0'$jinf, type="l", col="yellow3", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'1_1_0'$median), all_timestojthinf$'1_1_0'$median), c(rev(all_timestojthinf$'1_1_0'$low), all_timestojthinf$'1_1_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'1_1_1'$median, all_timestojthinf$'1_1_1'$jinf, type="l", col="blue", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'1_1_1'$median), all_timestojthinf$'1_1_1'$median), c(rev(all_timestojthinf$'1_1_1'$low), all_timestojthinf$'1_1_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0_0'$median, all_timestojthinf$'0_0'$jinf, type="l", col="black", lwd=2)
lines(all_timestojthinf$'0_1'$median, all_timestojthinf$'0_1'$jinf, type="l", col="orange", lwd=2)
lines(all_timestojthinf$'1_0'$median, all_timestojthinf$'1_0'$jinf, type="l", col="orchid", lwd=2)

plot(all_timestojthinf$'3_1_-1'$median, all_timestojthinf$'3_1_-1'$jinf, type="l", col="red", 
     xlab="Time", ylab="", xlim=c(0,200), ylim=c(0,1000), cex.axis=1.3, cex.lab=1.5, lwd=2)
#polygon(c(rev(all_timestojthinf$'3_1_-1'$median), all_timestojthinf$'3_1_-1'$median), c(rev(all_timestojthinf$'3_1_-1'$low), all_timestojthinf$'3_1_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'3_1_0'$median, all_timestojthinf$'3_1_0'$jinf, type="l", col="yellow3", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'3_1_0'$median), all_timestojthinf$'3_1_0'$median), c(rev(all_timestojthinf$'3_1_0'$low), all_timestojthinf$'3_1_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'3_1_1'$median, all_timestojthinf$'3_1_1'$jinf, type="l", col="blue", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'3_1_1'$median), all_timestojthinf$'3_1_1'$median), c(rev(all_timestojthinf$'3_1_1'$low), all_timestojthinf$'3_1_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0_0'$median, all_timestojthinf$'0_0'$jinf, type="l", col="black", lwd=2)
lines(all_timestojthinf$'0_1'$median, all_timestojthinf$'0_1'$jinf, type="l", col="orange", lwd=2)
lines(all_timestojthinf$'3_0'$median, all_timestojthinf$'3_0'$jinf, type="l", col="orchid", lwd=2)


# HiT=0.5
plot(all_timestojthinf$'0.5_0.5_-1'$median, all_timestojthinf$'0.5_0.5_-1'$jinf, type="l", col="red", 
     xlab="Time", ylab="Number of infections", xlim=c(0,200), ylim=c(0,1000), cex.axis=1.3, cex.lab=1.5, lwd=2)
#polygon(c(rev(all_timestojthinf$'0.5_0.5_-1'$median), all_timestojthinf$'0.5_0.5_-1'$median), c(rev(all_timestojthinf$'0.5_0.5_-1'$low), all_timestojthinf$'0.5_0.5_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0.5_0.5_0'$median, all_timestojthinf$'0.5_0.5_0'$jinf, type="l", col="yellow3", lwd=2)
#polygon(c(rev(all_timestojthinf$'0.5_0.5_0'$median), all_timestojthinf$'0.5_0.5_0'$median), c(rev(all_timestojthinf$'0.5_0.5_0'$low), all_timestojthinf$'0.5_0.5_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0.5_0.5_1'$median, all_timestojthinf$'0.5_0.5_1'$jinf, type="l", col="blue", lwd=2)
#polygon(c(rev(all_timestojthinf$'0.5_0.5_1'$median), all_timestojthinf$'0.5_0.5_1'$median), c(rev(all_timestojthinf$'0.5_0.5_1'$low), all_timestojthinf$'0.5_0.5_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0_0'$median, all_timestojthinf$'0_0'$jinf, type="l", col="black", lwd=2)
lines(all_timestojthinf$'0_0.5'$median, all_timestojthinf$'0_0.5'$jinf, type="l", col="orange", lwd=2)
lines(all_timestojthinf$'0.5_0'$median, all_timestojthinf$'0.5_0'$jinf, type="l", col="orchid", lwd=2)

plot(all_timestojthinf$'1_0.5_-1'$median, all_timestojthinf$'1_0.5_-1'$jinf, type="l", col="red", lty="solid", 
     xlab="Time", ylab="", xlim=c(0,200), ylim=c(0,1000), cex.axis=1.3, cex.lab=1.5, lwd=2)
#polygon(c(rev(all_timestojthinf$'1_0.5_-1'$median), all_timestojthinf$'1_0.5_-1'$median), c(rev(all_timestojthinf$'1_0.5_-1'$low), all_timestojthinf$'1_0.5_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'1_0.5_0'$median, all_timestojthinf$'1_0.5_0'$jinf, type="l", col="yellow3", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'1_0.5_0'$median), all_timestojthinf$'1_0.5_0'$median), c(rev(all_timestojthinf$'1_0.5_0'$low), all_timestojthinf$'1_0.5_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'1_0.5_1'$median, all_timestojthinf$'1_0.5_1'$jinf, type="l", col="blue", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'1_0.5_1'$median), all_timestojthinf$'1_0.5_1'$median), c(rev(all_timestojthinf$'1_0.5_1'$low), all_timestojthinf$'1_0.5_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0_0'$median, all_timestojthinf$'0_0'$jinf, type="l", col="black", lwd=2)
lines(all_timestojthinf$'0_0.5'$median, all_timestojthinf$'0_0.5'$jinf, type="l", col="orange", lwd=2)
lines(all_timestojthinf$'1_0'$median, all_timestojthinf$'1_0'$jinf, type="l", col="orchid", lwd=2)

plot(all_timestojthinf$'3_0.5_-1'$median, all_timestojthinf$'3_0.5_-1'$jinf, type="l", col="red", lty="solid", 
     xlab="Time", ylab="", xlim=c(0,200), ylim=c(0,1000), cex.axis=1.3, cex.lab=1.5, lwd=2)
#polygon(c(rev(all_timestojthinf$'3_0.5_-1'$median), all_timestojthinf$'3_0.5_-1'$median), c(rev(all_timestojthinf$'3_0.5_-1'$low), all_timestojthinf$'3_0.5_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'3_0.5_0'$median, all_timestojthinf$'3_0.5_0'$jinf, type="l", col="yellow3", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'3_0.5_0'$median), all_timestojthinf$'3_0.5_0'$median), c(rev(all_timestojthinf$'3_0.5_0'$low), all_timestojthinf$'3_0.5_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'3_0.5_1'$median, all_timestojthinf$'3_0.5_1'$jinf, type="l", col="blue", lty="solid", lwd=2)
#polygon(c(rev(all_timestojthinf$'3_0.5_1'$median), all_timestojthinf$'3_0.5_1'$median), c(rev(all_timestojthinf$'3_0.5_1'$low), all_timestojthinf$'3_0.5_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_timestojthinf$'0_0'$median, all_timestojthinf$'0_0'$jinf, type="l", col="black", lwd=2)
lines(all_timestojthinf$'0_0.5'$median, all_timestojthinf$'0_0.5'$jinf, type="l", col="orange", lwd=2)
lines(all_timestojthinf$'3_0'$median, all_timestojthinf$'3_0'$jinf, type="l", col="orchid", lwd=2)


par(mai=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top", bty="n", inset=0, horiz=T, legend=c("Hom","HiT only", "HiS only","Neg cor","No cor","Pos cor"),
       col=c("black","orange","orchid","red","yellow3","blue"),
       lty="solid", lwd=2, cex=1.5)

par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))

# 10 x 9


###################### plot R effective over time #################################
par(mfrow=c(3,3))

# HiS=0.5, HiT change
plot(all_Reffs$'0.5_0.5_-1'$time, all_Reffs$'0.5_0.5_-1'$avg, type="l", col="red", 
     xlab="Time", ylab="R effective", ylim=c(0,6))
#polygon(c(rev(all_Reffs$'0.5_0.5_-1'$time), all_Reffs$'0.5_0.5_-1'$time), c(rev(all_Reffs$'0.5_0.5_-1'$low), all_Reffs$'0.5_0.5_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0.5_0.5_0'$time, all_Reffs$'0.5_0.5_0'$avg, type="l", col="yellow3")
#polygon(c(rev(all_Reffs$'0.5_0.5_0'$time), all_Reffs$'0.5_0.5_0'$time), c(rev(all_Reffs$'0.5_0.5_0'$low), all_Reffs$'0.5_0.5_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0.5_0.5_1'$time, all_Reffs$'0.5_0.5_1'$avg, type="l", col="blue")
#polygon(c(rev(all_Reffs$'0.5_0.5_1'$time), all_Reffs$'0.5_0.5_1'$time), c(rev(all_Reffs$'0.5_0.5_1'$low), all_Reffs$'0.5_0.5_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0_0'$time, all_Reffs$'0_0'$avg, type="l", col="black")

plot(all_Reffs$'0.5_1_-1'$time, all_Reffs$'0.5_1_-1'$avg, type="l", col="red", lty="solid", 
      xlab="Time", ylab="R effective", ylim=c(0,6))
#polygon(c(rev(all_Reffs$'0.5_1_-1'$time), all_Reffs$'0.5_1_-1'$time), c(rev(all_Reffs$'0.5_1_-1'$low), all_Reffs$'0.5_1_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0.5_1_0'$time, all_Reffs$'0.5_1_0'$avg, type="l", col="yellow3", lty="solid")
#polygon(c(rev(all_Reffs$'0.5_1_0'$time), all_Reffs$'0.5_1_0'$time), c(rev(all_Reffs$'0.5_1_0'$low), all_Reffs$'0.5_1_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0.5_1_1'$time, all_Reffs$'0.5_1_1'$avg, type="l", col="blue", lty="solid")
#polygon(c(rev(all_Reffs$'0.5_1_1'$time), all_Reffs$'0.5_1_1'$time), c(rev(all_Reffs$'0.5_1_1'$low), all_Reffs$'0.5_1_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0_0'$time, all_Reffs$'0_0'$avg, type="l", col="black")

plot(all_Reffs$'0.5_3_-1'$time, all_Reffs$'0.5_3_-1'$avg, type="l", col="red", lty="solid", 
      xlab="Time", ylab="R effective", ylim=c(0,6))
#polygon(c(rev(all_Reffs$'0.5_3_-1'$time), all_Reffs$'0.5_3_-1'$time), c(rev(all_Reffs$'0.5_3_-1'$low), all_Reffs$'0.5_3_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0.5_3_0'$time, all_Reffs$'0.5_3_0'$avg, type="l", col="yellow3", lty="solid")
#polygon(c(rev(all_Reffs$'0.5_3_0'$time), all_Reffs$'0.5_3_0'$time), c(rev(all_Reffs$'0.5_3_0'$low), all_Reffs$'0.5_3_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0.5_3_1'$time, all_Reffs$'0.5_3_1'$avg, type="l", col="blue", lty="solid")
#polygon(c(rev(all_Reffs$'0.5_3_1'$time), all_Reffs$'0.5_3_1'$time), c(rev(all_Reffs$'0.5_3_1'$low), all_Reffs$'0.5_3_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0_0'$time, all_Reffs$'0_0'$avg, type="l", col="black")


# HiS=1
plot(all_Reffs$'1_0.5_-1'$time, all_Reffs$'1_0.5_-1'$avg, type="l", col="red", lty="solid", 
      xlab="Time", ylab="R effective", ylim=c(0,6))
#polygon(c(rev(all_Reffs$'1_0.5_-1'$time), all_Reffs$'1_0.5_-1'$time), c(rev(all_Reffs$'1_0.5_-1'$low), all_Reffs$'1_0.5_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'1_0.5_0'$time, all_Reffs$'1_0.5_0'$avg, type="l", col="yellow3", lty="solid")
#polygon(c(rev(all_Reffs$'1_0.5_0'$time), all_Reffs$'1_0.5_0'$time), c(rev(all_Reffs$'1_0.5_0'$low), all_Reffs$'1_0.5_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'1_0.5_1'$time, all_Reffs$'1_0.5_1'$avg, type="l", col="blue", lty="solid")
#polygon(c(rev(all_Reffs$'1_0.5_1'$time), all_Reffs$'1_0.5_1'$time), c(rev(all_Reffs$'1_0.5_1'$low), all_Reffs$'1_0.5_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0_0'$time, all_Reffs$'0_0'$avg, type="l", col="black")

plot(all_Reffs$'1_1_-1'$time, all_Reffs$'1_1_-1'$avg, type="l", col="red", lty="solid", 
     xlab="Time", ylab="R effective", ylim=c(0,6))
#polygon(c(rev(all_Reffs$'1_1_-1'$time), all_Reffs$'1_1_-1'$time), c(rev(all_Reffs$'1_1_-1'$low), all_Reffs$'1_1_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'1_1_0'$time, all_Reffs$'1_1_0'$avg, type="l", col="yellow3", lty="solid")
#polygon(c(rev(all_Reffs$'1_1_0'$time), all_Reffs$'1_1_0'$time), c(rev(all_Reffs$'1_1_0'$low), all_Reffs$'1_1_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'1_1_1'$time, all_Reffs$'1_1_1'$avg, type="l", col="blue", lty="solid")
#polygon(c(rev(all_Reffs$'1_1_1'$time), all_Reffs$'1_1_1'$time), c(rev(all_Reffs$'1_1_1'$low), all_Reffs$'1_1_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0_0'$time, all_Reffs$'0_0'$avg, type="l", col="black")

plot(all_Reffs$'1_3_-1'$time, all_Reffs$'1_3_-1'$avg, type="l", col="red", lty="solid", 
     xlab="Time", ylab="R effective", ylim=c(0,6))
#polygon(c(rev(all_Reffs$'1_3_-1'$time), all_Reffs$'1_3_-1'$time), c(rev(all_Reffs$'1_3_-1'$low), all_Reffs$'1_3_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'1_3_0'$time, all_Reffs$'1_3_0'$avg, type="l", col="yellow3", lty="solid")
#polygon(c(rev(all_Reffs$'1_3_0'$time), all_Reffs$'1_3_0'$time), c(rev(all_Reffs$'1_3_0'$low), all_Reffs$'1_3_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'1_3_1'$time, all_Reffs$'1_3_1'$avg, type="l", col="blue", lty="solid")
#polygon(c(rev(all_Reffs$'1_3_1'$time), all_Reffs$'1_3_1'$time), c(rev(all_Reffs$'1_3_1'$low), all_Reffs$'1_3_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0_0'$time, all_Reffs$'0_0'$avg, type="l", col="black")


# HiS=3
plot(all_Reffs$'3_0.5_-1'$time, all_Reffs$'3_0.5_-1'$avg, type="l", col="red", lty="solid", 
      xlab="Time", ylab="R effective", ylim=c(0,6))
#polygon(c(rev(all_Reffs$'3_0.5_-1'$time), all_Reffs$'3_0.5_-1'$time), c(rev(all_Reffs$'3_0.5_-1'$low), all_Reffs$'3_0.5_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'3_0.5_0'$time, all_Reffs$'3_0.5_0'$avg, type="l", col="yellow3", lty="solid")
#polygon(c(rev(all_Reffs$'3_0.5_0'$time), all_Reffs$'3_0.5_0'$time), c(rev(all_Reffs$'3_0.5_0'$low), all_Reffs$'3_0.5_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'3_0.5_1'$time, all_Reffs$'3_0.5_1'$avg, type="l", col="blue", lty="solid")
#polygon(c(rev(all_Reffs$'3_0.5_1'$time), all_Reffs$'3_0.5_1'$time), c(rev(all_Reffs$'3_0.5_1'$low), all_Reffs$'3_0.5_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0_0'$time, all_Reffs$'0_0'$avg, type="l", col="black")

plot(all_Reffs$'3_1_-1'$time, all_Reffs$'3_1_-1'$avg, type="l", col="red", 
     xlab="Time", ylab="R effective", ylim=c(0,6))
#polygon(c(rev(all_Reffs$'3_1_-1'$time), all_Reffs$'3_1_-1'$time), c(rev(all_Reffs$'3_1_-1'$low), all_Reffs$'3_1_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'3_1_0'$time, all_Reffs$'3_1_0'$avg, type="l", col="yellow3", lty="solid")
#polygon(c(rev(all_Reffs$'3_1_0'$time), all_Reffs$'3_1_0'$time), c(rev(all_Reffs$'3_1_0'$low), all_Reffs$'3_1_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'3_1_1'$time, all_Reffs$'3_1_1'$avg, type="l", col="blue", lty="solid")
#polygon(c(rev(all_Reffs$'3_1_1'$time), all_Reffs$'3_1_1'$time), c(rev(all_Reffs$'3_1_1'$low), all_Reffs$'3_1_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0_0'$time, all_Reffs$'0_0'$avg, type="l", col="black")

plot(all_Reffs$'3_3_-1'$time, all_Reffs$'3_3_-1'$avg, type="l", col="red", 
     xlab="Time", ylab="R effective", ylim=c(0,6))
#polygon(c(rev(all_Reffs$'3_3_-1'$time), all_Reffs$'3_3_-1'$time), c(rev(all_Reffs$'3_3_-1'$low), all_Reffs$'3_3_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'3_3_0'$time, all_Reffs$'3_3_0'$avg, type="l", col="yellow3")
#polygon(c(rev(all_Reffs$'3_3_0'$time), all_Reffs$'3_3_0'$time), c(rev(all_Reffs$'3_3_0'$low), all_Reffs$'3_3_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'3_3_1'$time, all_Reffs$'3_3_1'$avg, type="l", col="blue")
#polygon(c(rev(all_Reffs$'3_3_1'$time), all_Reffs$'3_3_1'$time), c(rev(all_Reffs$'3_3_1'$low), all_Reffs$'3_3_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Reffs$'0_0'$time, all_Reffs$'0_0'$avg, type="l", col="black")

par(mfrow=c(1,1))


###########################################################################################
#### plot phase planes R effective vs Susceptibles ########



layout(mat=matrix(c(1,2,3,
                    4,5,6,
                    7,8,9,
                    10,10,10), nrow=4, ncol=3, byrow=T),heights = c(0.4,0.4,0.4,0.05))


par(mai = c(0.55, 0.6, 0.07, 0.075))


# HiT=3
plot(all_Sseries$'0_0'$avg, all_Reffs$'0_0'$avg, type="l", col="black",lwd=2,
     xlab="", ylab=expression(paste(italic(R)[italic(e)])), xlim=c(1000,0), ylim=c(0,10), cex.axis=1.3, cex.lab=1.5)
#polygon(c(rev(all_Sseries$'0_0'$avg), all_Sseries$'0_0'$avg), c(rev(all_Reffs$'0_0'$low), all_Reffs$'0_0'$high),
#        col=adjustcolor("black", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0_0'$avg[seq(1,max(all_Sseries$'0_0'$time/0.1),100)], all_Reffs$'0_0'$avg[seq(1,max(all_Reffs$'0_0'$time/0.1),100)], labels=all_Reffs$'0_0'$time[seq(1,max(all_Reffs$'0_0'$time/0.1),100)]/10, col="black", cex=2)
lines(all_Sseries$'0_3'$avg, all_Reffs$'0_3'$avg, col="orange", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'0_3'$avg), all_Sseries$'0_3'$avg), c(rev(all_Reffs$'0_3'$low), all_Reffs$'0_3'$high),
#        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0_3'$avg[seq(1,max(all_Sseries$'0_3'$time/0.1),100)], all_Reffs$'0_3'$avg[seq(1,max(all_Reffs$'0_3'$time/0.1),100)], labels=all_Reffs$'0_3'$time[seq(1,max(all_Reffs$'0_3'$time/0.1),100)]/10, col="orange", cex=2)
lines(all_Sseries$'0.5_0'$avg, all_Reffs$'0.5_0'$avg, col="orchid", lty="solid",lwd=2)
text(all_Sseries$'0.5_0'$avg[seq(1,max(all_Sseries$'0.5_0'$time/0.1),100)], all_Reffs$'0.5_0'$avg[seq(1,max(all_Reffs$'0.5_0'$time/0.1),100)], labels=all_Reffs$'0.5_0'$time[seq(1,max(all_Reffs$'0.5_0'$time/0.1),100)]/10, col="orchid", cex=2)
#polygon(c(rev(all_Sseries$'0.5_0'$avg), all_Sseries$'0.5_0'$avg), c(rev(all_Reffs$'0.5_0'$low), all_Reffs$'0.5_0'$high),
#        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'0.5_3_0'$avg, all_Reffs$'0.5_3_0'$avg, col="yellow3", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'0.5_3_0'$avg), all_Sseries$'0.5_3_0'$avg), c(rev(all_Reffs$'0.5_3_0'$low), all_Reffs$'0.5_3_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0.5_3_0'$avg[seq(1,max(all_Sseries$'0.5_3_0'$time/0.1),100)], all_Reffs$'0.5_3_0'$avg[seq(1,max(all_Reffs$'0.5_3_0'$time/0.1),100)], labels=all_Reffs$'0.5_3_0'$time[seq(1,max(all_Reffs$'0.5_3_0'$time/0.1),100)]/10, col="yellow3", cex=2)
lines(all_Sseries$'0.5_3_-1'$avg, all_Reffs$'0.5_3_-1'$avg, col="red",lwd=2)
#polygon(c(rev(all_Sseries$'0.5_3_-1'$avg), all_Sseries$'0.5_3_-1'$avg), c(rev(all_Reffs$'0.5_3_-1'$low), all_Reffs$'0.5_3_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0.5_3_-1'$avg[seq(1,max(all_Sseries$'0.5_3_-1'$time/0.1),100)], all_Reffs$'0.5_3_-1'$avg[seq(1,max(all_Reffs$'0.5_3_-1'$time/0.1),100)], labels=all_Reffs$'0.5_3_-1'$time[seq(1,max(all_Reffs$'0.5_3_-1'$time/0.1),100)]/10, col="red", cex=2)
lines(all_Sseries$'0.5_3_1'$avg, all_Reffs$'0.5_3_1'$avg, type="l", col="blue", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'0.5_3_1'$avg), all_Sseries$'0.5_3_1'$avg), c(rev(all_Reffs$'0.5_3_1'$low), all_Reffs$'0.5_3_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0.5_3_1'$avg[seq(1,max(all_Sseries$'0.5_3_1'$time/0.1),100)], all_Reffs$'0.5_3_1'$avg[seq(1,max(all_Reffs$'0.5_3_1'$time/0.1),100)], labels=all_Reffs$'0.5_3_1'$time[seq(1,max(all_Reffs$'0.5_3_1'$time/0.1),100)]/10, col="blue", cex=2)
abline(h=1, col="gray45", lty="dashed")

plot(all_Sseries$'0_0'$avg, all_Reffs$'0_0'$avg, type="l", col="black",lwd=2,
     xlab="", ylab="", xlim=c(1000,0), ylim=c(0,10), cex.axis=1.3, cex.lab=1.5)
#polygon(c(rev(all_Sseries$'0_0'$avg), all_Sseries$'0_0'$avg), c(rev(all_Reffs$'0_0'$low), all_Reffs$'0_0'$high),
#        col=adjustcolor("black", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0_0'$avg[seq(1,max(all_Sseries$'0_0'$time/0.1),100)], all_Reffs$'0_0'$avg[seq(1,max(all_Reffs$'0_0'$time/0.1),100)], labels=all_Reffs$'0_0'$time[seq(1,max(all_Reffs$'0_0'$time/0.1),100)]/10, col="black", cex=2)
lines(all_Sseries$'0_3'$avg, all_Reffs$'0_3'$avg, col="orange", lty="solid",lwd=2)
text(all_Sseries$'0_3'$avg[seq(1,max(all_Sseries$'0_3'$time/0.1),100)], all_Reffs$'0_3'$avg[seq(1,max(all_Reffs$'0_3'$time/0.1),100)], labels=all_Reffs$'0_3'$time[seq(1,max(all_Reffs$'0_3'$time/0.1),100)]/10, col="orange", cex=2)
#polygon(c(rev(all_Sseries$'0_3'$avg), all_Sseries$'0_3'$avg), c(rev(all_Reffs$'0_3'$low), all_Reffs$'0_3'$high),
#        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'1_0'$avg, all_Reffs$'1_0'$avg, col="orchid", lty="solid",lwd=2)
text(all_Sseries$'1_0'$avg[seq(1,max(all_Sseries$'1_0'$time/0.1),100)], all_Reffs$'1_0'$avg[seq(1,max(all_Reffs$'1_0'$time/0.1),100)], labels=all_Reffs$'1_0'$time[seq(1,max(all_Reffs$'1_0'$time/0.1),100)]/10, col="orchid", cex=2)
#polygon(c(rev(all_Sseries$'1_0'$avg), all_Sseries$'1_0'$avg), c(rev(all_Reffs$'1_0'$low), all_Reffs$'1_0'$high),
#        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'1_3_0'$avg, all_Reffs$'1_3_0'$avg, col="yellow3", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'1_3_0'$avg), all_Sseries$'1_3_0'$avg), c(rev(all_Reffs$'1_3_0'$low), all_Reffs$'1_3_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
text(all_Sseries$'1_3_0'$avg[seq(1,max(all_Sseries$'1_3_0'$time/0.1),100)], all_Reffs$'1_3_0'$avg[seq(1,max(all_Reffs$'1_3_0'$time/0.1),100)], labels=all_Reffs$'1_3_0'$time[seq(1,max(all_Reffs$'1_3_0'$time/0.1),100)]/10, col="yellow3", cex=2)
lines(all_Sseries$'1_3_-1'$avg, all_Reffs$'1_3_-1'$avg, col="red",lwd=2)
#polygon(c(rev(all_Sseries$'1_3_-1'$avg), all_Sseries$'1_3_-1'$avg), c(rev(all_Reffs$'1_3_-1'$low), all_Reffs$'1_3_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
text(all_Sseries$'1_3_-1'$avg[seq(1,max(all_Sseries$'1_3_-1'$time/0.1),100)], all_Reffs$'1_3_-1'$avg[seq(1,max(all_Reffs$'1_3_-1'$time/0.1),100)], labels=all_Reffs$'1_3_-1'$time[seq(1,max(all_Reffs$'1_3_-1'$time/0.1),100)]/10, col="red", cex=2)
lines(all_Sseries$'1_3_1'$avg, all_Reffs$'1_3_1'$avg, type="l", col="blue", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'1_3_1'$avg), all_Sseries$'1_3_1'$avg), c(rev(all_Reffs$'1_3_1'$low), all_Reffs$'1_3_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
text(all_Sseries$'1_3_1'$avg[seq(1,max(all_Sseries$'1_3_1'$time/0.1),100)], all_Reffs$'1_3_1'$avg[seq(1,max(all_Reffs$'1_3_1'$time/0.1),100)], labels=all_Reffs$'1_3_1'$time[seq(1,max(all_Reffs$'1_3_1'$time/0.1),100)]/10, col="blue", cex=2)
abline(h=1, col="gray75", lty="dashed")

plot(all_Sseries$'0_0'$avg, all_Reffs$'0_0'$avg, type="l", col="black",lwd=2,
     xlab="", ylab="", xlim=c(1000,0), ylim=c(0,10), cex.axis=1.3, cex.lab=1.5)
#polygon(c(rev(all_Sseries$'0_0'$avg), all_Sseries$'0_0'$avg), c(rev(all_Reffs$'0_0'$low), all_Reffs$'0_0'$high),
#        col=adjustcolor("black", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0_0'$avg[seq(1,max(all_Sseries$'0_0'$time/0.1),100)], all_Reffs$'0_0'$avg[seq(1,max(all_Reffs$'0_0'$time/0.1),100)], labels=all_Reffs$'0_0'$time[seq(1,max(all_Reffs$'0_0'$time/0.1),100)]/10, col="black", cex=2)
lines(all_Sseries$'0_3'$avg, all_Reffs$'0_3'$avg, col="orange", lty="solid",lwd=2)
text(all_Sseries$'0_3'$avg[seq(1,max(all_Sseries$'0_3'$time/0.1),100)], all_Reffs$'0_3'$avg[seq(1,max(all_Reffs$'0_3'$time/0.1),100)], labels=all_Reffs$'0_3'$time[seq(1,max(all_Reffs$'0_3'$time/0.1),100)]/10, col="orange", cex=2)
#polygon(c(rev(all_Sseries$'0_3'$avg), all_Sseries$'0_3'$avg), c(rev(all_Reffs$'0_3'$low), all_Reffs$'0_3'$high),
#        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_0'$avg, all_Reffs$'3_0'$avg, col="orchid", lty="solid",lwd=2)
text(all_Sseries$'3_0'$avg[seq(1,max(all_Sseries$'3_0'$time/0.1),100)], all_Reffs$'3_0'$avg[seq(1,max(all_Reffs$'3_0'$time/0.1),100)], labels=all_Reffs$'3_0'$time[seq(1,max(all_Reffs$'3_0'$time/0.1),100)]/10, col="orchid", cex=2)
#polygon(c(rev(all_Sseries$'3_0'$avg), all_Sseries$'3_0'$avg), c(rev(all_Reffs$'3_0'$low), all_Reffs$'3_0'$high),
#        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_3_0'$avg, all_Reffs$'3_3_0'$avg, col="yellow3", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'3_3_0'$avg), all_Sseries$'3_3_0'$avg), c(rev(all_Reffs$'3_3_0'$low), all_Reffs$'3_3_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
text(all_Sseries$'3_3_0'$avg[seq(1,max(all_Sseries$'3_3_0'$time/0.1),100)], all_Reffs$'3_3_0'$avg[seq(1,max(all_Reffs$'3_3_0'$time/0.1),100)], labels=all_Reffs$'3_3_0'$time[seq(1,max(all_Reffs$'3_3_0'$time/0.1),100)]/10, col="yellow3", cex=2)
lines(all_Sseries$'3_3_-1'$avg, all_Reffs$'3_3_-1'$avg, col="red",lwd=2)
#polygon(c(rev(all_Sseries$'3_3_-1'$avg), all_Sseries$'3_3_-1'$avg), c(rev(all_Reffs$'3_3_-1'$low), all_Reffs$'3_3_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
text(all_Sseries$'3_3_-1'$avg[seq(1,max(all_Sseries$'3_3_-1'$time/0.1),100)], all_Reffs$'3_3_-1'$avg[seq(1,max(all_Reffs$'3_3_-1'$time/0.1),100)], labels=all_Reffs$'3_3_-1'$time[seq(1,max(all_Reffs$'3_3_-1'$time/0.1),100)]/10, col="red", cex=2)
lines(all_Sseries$'3_3_1'$avg, all_Reffs$'3_3_1'$avg, type="l", col="blue", lty="solid", lwd=2)
#polygon(c(rev(all_Sseries$'3_3_1'$avg), all_Sseries$'3_3_1'$avg), c(rev(all_Reffs$'3_3_1'$low), all_Reffs$'3_3_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
text(all_Sseries$'3_3_1'$avg[seq(1,max(all_Sseries$'3_3_1'$time/0.1),100)], all_Reffs$'3_3_1'$avg[seq(1,max(all_Reffs$'3_3_1'$time/0.1),100)], labels=all_Reffs$'3_3_1'$time[seq(1,max(all_Reffs$'3_3_1'$time/0.1),100)]/10, col="blue", cex=2)
abline(h=1, col="gray75", lty="dashed")



# HiT=1
plot(all_Sseries$'0_0'$avg, all_Reffs$'0_0'$avg, type="l", col="black",lwd=2,
     xlab="", ylab=expression(paste(italic(R)[italic(e)])), xlim=c(1000,0), ylim=c(0,10), cex.axis=1.3, cex.lab=1.5)
#polygon(c(rev(all_Sseries$'0_0'$avg), all_Sseries$'0_0'$avg), c(rev(all_Reffs$'0_0'$low), all_Reffs$'0_0'$high),
#        col=adjustcolor("black", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0_0'$avg[seq(1,max(all_Sseries$'0_0'$time/0.1),100)], all_Reffs$'0_0'$avg[seq(1,max(all_Reffs$'0_0'$time/0.1),100)], labels=all_Reffs$'0_0'$time[seq(1,max(all_Reffs$'0_0'$time/0.1),100)]/10, col="black", cex=2)
lines(all_Sseries$'0_1'$avg, all_Reffs$'0_1'$avg, col="orange", lty="solid",lwd=2)
text(all_Sseries$'0_1'$avg[seq(1,max(all_Sseries$'0_1'$time/0.1),100)], all_Reffs$'0_1'$avg[seq(1,max(all_Reffs$'0_1'$time/0.1),100)], labels=all_Reffs$'0_1'$time[seq(1,max(all_Reffs$'0_1'$time/0.1),100)]/10, col="orange", cex=2)
#polygon(c(rev(all_Sseries$'0_1'$avg), all_Sseries$'0_1'$avg), c(rev(all_Reffs$'0_1'$low), all_Reffs$'0_1'$high),
#        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'0.5_0'$avg, all_Reffs$'0.5_0'$avg, col="orchid", lty="solid",lwd=2)
text(all_Sseries$'0.5_0'$avg[seq(1,max(all_Sseries$'0.5_0'$time/0.1),100)], all_Reffs$'0.5_0'$avg[seq(1,max(all_Reffs$'0.5_0'$time/0.1),100)], labels=all_Reffs$'0.5_0'$time[seq(1,max(all_Reffs$'0.5_0'$time/0.1),100)]/10, col="orchid", cex=2)
#polygon(c(rev(all_Sseries$'0.5_0'$avg), all_Sseries$'0.5_0'$avg), c(rev(all_Reffs$'0.5_0'$low), all_Reffs$'0.5_0'$high),
#        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'0.5_1_0'$avg, all_Reffs$'0.5_1_0'$avg, col="yellow3", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'0.5_1_0'$avg), all_Sseries$'0.5_1_0'$avg), c(rev(all_Reffs$'0.5_1_0'$low), all_Reffs$'0.5_1_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0.5_1_0'$avg[seq(1,max(all_Sseries$'0.5_1_0'$time/0.1),100)], all_Reffs$'0.5_1_0'$avg[seq(1,max(all_Reffs$'0.5_1_0'$time/0.1),100)], all_Reffs$'0.5_1_0'$time[seq(1,max(all_Reffs$'0.5_1_0'$time/0.1),100)]/10, col="yellow3", cex=2)
lines(all_Sseries$'0.5_1_-1'$avg, all_Reffs$'0.5_1_-1'$avg, col="red",lwd=2)
#polygon(c(rev(all_Sseries$'0.5_1_-1'$avg), all_Sseries$'0.5_1_-1'$avg), c(rev(all_Reffs$'0.5_1_-1'$low), all_Reffs$'0.5_1_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0.5_1_-1'$avg[seq(1,max(all_Sseries$'0.5_1_-1'$time/0.1),100)], all_Reffs$'0.5_1_-1'$avg[seq(1,max(all_Reffs$'0.5_1_-1'$time/0.1),100)], labels=all_Reffs$'0.5_1_-1'$time[seq(1,max(all_Reffs$'0.5_1_-1'$time/0.1),100)]/10, col="red", cex=2)
lines(all_Sseries$'0.5_1_1'$avg, all_Reffs$'0.5_1_1'$avg, type="l", col="blue", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'0.5_1_1'$avg), all_Sseries$'0.5_1_1'$avg), c(rev(all_Reffs$'0.5_1_1'$low), all_Reffs$'0.5_1_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0.5_1_1'$avg[seq(1,max(all_Sseries$'0.5_1_1'$time/0.1),100)], all_Reffs$'0.5_1_1'$avg[seq(1,max(all_Reffs$'0.5_1_1'$time/0.1),100)], labels=all_Reffs$'0.5_1_1'$time[seq(1,max(all_Reffs$'0.5_1_1'$time/0.1),100)]/10, col="blue", cex=2)
abline(h=1, col="gray75", lty="dashed")

plot(all_Sseries$'0_0'$avg, all_Reffs$'0_0'$avg, type="l", col="black",lwd=2,
     xlab="", ylab="", xlim=c(1000,0), ylim=c(0,10), cex.axis=1.3, cex.lab=1.5)
#polygon(c(rev(all_Sseries$'0_0'$avg), all_Sseries$'0_0'$avg), c(rev(all_Reffs$'0_0'$low), all_Reffs$'0_0'$high),
#        col=adjustcolor("black", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0_0'$avg[seq(1,max(all_Sseries$'0_0'$time/0.1),100)], all_Reffs$'0_0'$avg[seq(1,max(all_Reffs$'0_0'$time/0.1),100)], labels=all_Reffs$'0_0'$time[seq(1,max(all_Reffs$'0_0'$time/0.1),100)]/10, col="black", cex=2)
lines(all_Sseries$'0_1'$avg, all_Reffs$'0_1'$avg, col="orange", lty="solid",lwd=2)
text(all_Sseries$'0_1'$avg[seq(1,max(all_Sseries$'0_1'$time/0.1),100)], all_Reffs$'0_1'$avg[seq(1,max(all_Reffs$'0_1'$time/0.1),100)], labels=all_Reffs$'0_1'$time[seq(1,max(all_Reffs$'0_1'$time/0.1),100)]/10, col="orange", cex=2)
#polygon(c(rev(all_Sseries$'0_1'$avg), all_Sseries$'0_1'$avg), c(rev(all_Reffs$'0_1'$low), all_Reffs$'0_1'$high),
#        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'1_0'$avg, all_Reffs$'1_0'$avg, col="orchid", lty="solid",lwd=2)
text(all_Sseries$'1_0'$avg[seq(1,max(all_Sseries$'1_0'$time/0.1),100)], all_Reffs$'1_0'$avg[seq(1,max(all_Reffs$'1_0'$time/0.1),100)], labels=all_Reffs$'1_0'$time[seq(1,max(all_Reffs$'1_0'$time/0.1),100)]/10, col="orchid", cex=2)
#polygon(c(rev(all_Sseries$'1_0'$avg), all_Sseries$'1_0'$avg), c(rev(all_Reffs$'1_0'$low), all_Reffs$'1_0'$high),
#        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'1_1_0'$avg, all_Reffs$'1_1_0'$avg, col="yellow3", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'1_1_0'$avg), all_Sseries$'1_1_0'$avg), c(rev(all_Reffs$'1_1_0'$low), all_Reffs$'1_1_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
text(all_Sseries$'1_1_0'$avg[seq(1,max(all_Sseries$'1_1_0'$time/0.1),100)], all_Reffs$'1_1_0'$avg[seq(1,max(all_Reffs$'1_1_0'$time/0.1),100)], labels=all_Reffs$'1_1_0'$time[seq(1,max(all_Reffs$'1_1_0'$time/0.1),100)]/10, col="yellow3", cex=2)
lines(all_Sseries$'1_1_-1'$avg, all_Reffs$'1_1_-1'$avg, col="red",lwd=2)
#polygon(c(rev(all_Sseries$'1_1_-1'$avg), all_Sseries$'1_1_-1'$avg), c(rev(all_Reffs$'1_1_-1'$low), all_Reffs$'1_1_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
text(all_Sseries$'1_1_-1'$avg[seq(1,max(all_Sseries$'1_1_-1'$time/0.1),100)], all_Reffs$'1_1_-1'$avg[seq(1,max(all_Reffs$'1_1_-1'$time/0.1),100)], labels=all_Reffs$'1_1_-1'$time[seq(1,max(all_Reffs$'1_1_-1'$time/0.1),100)]/10, col="red", cex=2)
lines(all_Sseries$'1_1_1'$avg, all_Reffs$'1_1_1'$avg, type="l", col="blue", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'1_1_1'$avg), all_Sseries$'1_1_1'$avg), c(rev(all_Reffs$'1_1_1'$low), all_Reffs$'1_1_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
text(all_Sseries$'1_1_1'$avg[seq(1,max(all_Sseries$'1_1_1'$time/0.1),100)], all_Reffs$'1_1_1'$avg[seq(1,max(all_Reffs$'1_1_1'$time/0.1),100)], labels=all_Reffs$'1_1_1'$time[seq(1,max(all_Reffs$'1_1_1'$time/0.1),100)]/10, col="blue", cex=2)
abline(h=1, col="gray75", lty="dashed")

plot(all_Sseries$'0_0'$avg, all_Reffs$'0_0'$avg, type="l", col="black",lwd=2,
     xlab="", ylab="", xlim=c(1000,0), ylim=c(0,10), cex.axis=1.3, cex.lab=1.5)
#polygon(c(rev(all_Sseries$'0_0'$avg), all_Sseries$'0_0'$avg), c(rev(all_Reffs$'0_0'$low), all_Reffs$'0_0'$high),
#        col=adjustcolor("black", alpha.f=0.2) , lty = 0)
#text(all_Sseries$'0_0'$avg[seq(1,max(all_Sseries$'0_0'$time/0.1),100)], all_Reffs$'0_0'$avg[seq(1,max(all_Reffs$'0_0'$time/0.1),100)], col="black", cex=2)
text(all_Sseries$'0_0'$avg[seq(1,max(all_Sseries$'0_0'$time/0.1),100)], all_Reffs$'0_0'$avg[seq(1,max(all_Reffs$'0_0'$time/0.1),100)], labels=all_Reffs$'0_0'$time[seq(1,max(all_Reffs$'0_0'$time/0.1),100)]/10, col="black", cex=2)
lines(all_Sseries$'0_1'$avg, all_Reffs$'0_1'$avg, col="orange", lty="solid",lwd=2)
text(all_Sseries$'0_1'$avg[seq(1,max(all_Sseries$'0_1'$time/0.1),100)], all_Reffs$'0_1'$avg[seq(1,max(all_Reffs$'0_1'$time/0.1),100)], labels=all_Reffs$'0_1'$time[seq(1,max(all_Reffs$'0_1'$time/0.1),100)]/10, col="orange", cex=2)
#polygon(c(rev(all_Sseries$'0_1'$avg), all_Sseries$'0_1'$avg), c(rev(all_Reffs$'0_1'$low), all_Reffs$'0_1'$high),
#        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_0'$avg, all_Reffs$'3_0'$avg, col="orchid", lty="solid",lwd=2)
text(all_Sseries$'3_0'$avg[seq(1,max(all_Sseries$'3_0'$time/0.1),100)], all_Reffs$'3_0'$avg[seq(1,max(all_Reffs$'3_0'$time/0.1),100)], labels=all_Reffs$'3_0'$time[seq(1,max(all_Reffs$'3_0'$time/0.1),100)]/10, col="orchid", cex=2)
#polygon(c(rev(all_Sseries$'3_0'$avg), all_Sseries$'3_0'$avg), c(rev(all_Reffs$'3_0'$low), all_Reffs$'3_0'$high),
#        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_1_0'$avg, all_Reffs$'3_1_0'$avg, col="yellow3", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'3_1_0'$avg), all_Sseries$'3_1_0'$avg), c(rev(all_Reffs$'3_1_0'$low), all_Reffs$'3_1_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
text(all_Sseries$'3_1_0'$avg[seq(1,max(all_Sseries$'3_1_0'$time/0.1),100)], all_Reffs$'3_1_0'$avg[seq(1,max(all_Reffs$'3_1_0'$time/0.1),100)], labels=all_Reffs$'3_1_0'$time[seq(1,max(all_Reffs$'3_1_0'$time/0.1),100)]/10, col="yellow3", cex=2)
lines(all_Sseries$'3_1_-1'$avg, all_Reffs$'3_1_-1'$avg, col="red",lwd=2)
#polygon(c(rev(all_Sseries$'3_1_-1'$avg), all_Sseries$'3_1_-1'$avg), c(rev(all_Reffs$'3_1_-1'$low), all_Reffs$'3_1_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
text(all_Sseries$'3_1_-1'$avg[seq(1,max(all_Sseries$'3_1_-1'$time/0.1),100)], all_Reffs$'3_1_-1'$avg[seq(1,max(all_Reffs$'3_1_-1'$time/0.1),100)], labels=all_Reffs$'3_1_-1'$time[seq(1,max(all_Reffs$'3_1_-1'$time/0.1),100)]/10, col="red", cex=2)
lines(all_Sseries$'3_1_1'$avg, all_Reffs$'3_1_1'$avg, type="l", col="blue", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'3_1_1'$avg), all_Sseries$'3_1_1'$avg), c(rev(all_Reffs$'3_1_1'$low), all_Reffs$'3_1_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
text(all_Sseries$'3_1_1'$avg[seq(1,max(all_Sseries$'3_1_1'$time/0.1),100)], all_Reffs$'3_1_1'$avg[seq(1,max(all_Reffs$'3_1_1'$time/0.1),100)], labels=all_Reffs$'3_1_1'$time[seq(1,max(all_Reffs$'3_1_1'$time/0.1),100)]/10, col="blue", cex=2)
abline(h=1, col="gray75", lty="dashed")


# HiT=0.5
plot(all_Sseries$'0_0'$avg, all_Reffs$'0_0'$avg, type="l", col="black",lwd=2,
     xlab="Number of susceptible individuals", ylab=expression(paste(italic(R)[italic(e)])), xlim=c(1000,0), ylim=c(0,10), cex.axis=1.3, cex.lab=1.5)
#polygon(c(rev(all_Sseries$'0_0'$avg), all_Sseries$'0_0'$avg), c(rev(all_Reffs$'0_0'$low), all_Reffs$'0_0'$high),
#        col=adjustcolor("black", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0_0'$avg[seq(1,max(all_Sseries$'0_0'$time/0.1),100)], all_Reffs$'0_0'$avg[seq(1,max(all_Reffs$'0_0'$time/0.1),100)], labels=all_Reffs$'0_0'$time[seq(1,max(all_Reffs$'0_0'$time/0.1),100)]/10, col="black", cex=2)
lines(all_Sseries$'0_0.5'$avg, all_Reffs$'0_0.5'$avg, col="orange", lty="solid",lwd=2)
text(all_Sseries$'0_0.5'$avg[seq(1,max(all_Sseries$'0_0.5'$time/0.1),100)], all_Reffs$'0_0.5'$avg[seq(1,max(all_Reffs$'0_0.5'$time/0.1),100)], labels=all_Reffs$'0_0.5'$time[seq(1,max(all_Reffs$'0_0.5'$time/0.1),100)]/10, col="orange", cex=2)
#polygon(c(rev(all_Sseries$'0_0.5'$avg), all_Sseries$'0_0.5'$avg), c(rev(all_Reffs$'0_0.5'$low), all_Reffs$'0_0.5'$high),
#        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'0.5_0'$avg, all_Reffs$'0.5_0'$avg, col="orchid", lty="solid",lwd=2)
text(all_Sseries$'0.5_0'$avg[seq(1,max(all_Sseries$'0.5_0'$time/0.1),100)], all_Reffs$'0.5_0'$avg[seq(1,max(all_Reffs$'0.5_0'$time/0.1),100)], labels=all_Reffs$'0.5_0'$time[seq(1,max(all_Reffs$'0.5_0'$time/0.1),100)]/10, col="orchid", cex=2)
#polygon(c(rev(all_Sseries$'0.5_0'$avg), all_Sseries$'0.5_0'$avg), c(rev(all_Reffs$'0.5_0'$low), all_Reffs$'0.5_0'$high),
#        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'0.5_0.5_0'$avg, all_Reffs$'0.5_0.5_0'$avg, col="yellow3", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'0.5_0.5_0'$avg), all_Sseries$'0.5_0.5_0'$avg), c(rev(all_Reffs$'0.5_0.5_0'$low), all_Reffs$'0.5_0.5_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0.5_0.5_0'$avg[seq(1,max(all_Sseries$'0.5_0.5_0'$time/0.1),100)], all_Reffs$'0.5_0.5_0'$avg[seq(1,max(all_Reffs$'0.5_0.5_0'$time/0.1),100)], labels=all_Reffs$'0.5_0.5_0'$time[seq(1,max(all_Reffs$'0.5_0.5_0'$time/0.1),100)]/10, col="yellow3", cex=2)
lines(all_Sseries$'0.5_0.5_-1'$avg, all_Reffs$'0.5_0.5_-1'$avg, col="red",lwd=2)
#polygon(c(rev(all_Sseries$'0.5_0.5_-1'$avg), all_Sseries$'0.5_0.5_-1'$avg), c(rev(all_Reffs$'0.5_0.5_-1'$low), all_Reffs$'0.5_0.5_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0.5_0.5_-1'$avg[seq(1,max(all_Sseries$'0.5_0.5_-1'$time/0.1),100)], all_Reffs$'0.5_0.5_-1'$avg[seq(1,max(all_Reffs$'0.5_0.5_-1'$time/0.1),100)], labels=all_Reffs$'0.5_0.5_-1'$time[seq(1,max(all_Reffs$'0.5_0.5_-1'$time/0.1),100)]/10, col="red", cex=2)
lines(all_Sseries$'0.5_0.5_1'$avg, all_Reffs$'0.5_0.5_1'$avg, type="l", col="blue", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'0.5_0.5_1'$avg), all_Sseries$'0.5_0.5_1'$avg), c(rev(all_Reffs$'0.5_0.5_1'$low), all_Reffs$'0.5_0.5_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0.5_0.5_1'$avg[seq(1,max(all_Sseries$'0.5_0.5_1'$time/0.1),100)], all_Reffs$'0.5_0.5_1'$avg[seq(1,max(all_Reffs$'0.5_0.5_1'$time/0.1),100)], labels=all_Reffs$'0.5_0.5_1'$time[seq(1,max(all_Reffs$'0.5_0.5_1'$time/0.1),100)]/10, col="blue", cex=2)
abline(h=1, col="gray75", lty="dashed")

plot(all_Sseries$'0_0'$avg, all_Reffs$'0_0'$avg, type="l", col="black",lwd=2,
     xlab="Number of susceptible individuals", ylab="", xlim=c(1000,0), ylim=c(0,10), cex.axis=1.3, cex.lab=1.5)
#polygon(c(rev(all_Sseries$'0_0'$avg), all_Sseries$'0_0'$avg), c(rev(all_Reffs$'0_0'$low), all_Reffs$'0_0'$high),
#        col=adjustcolor("black", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0_0'$avg[seq(1,max(all_Sseries$'0_0'$time/0.1),100)], all_Reffs$'0_0'$avg[seq(1,max(all_Reffs$'0_0'$time/0.1),100)], labels=all_Reffs$'0_0'$time[seq(1,max(all_Reffs$'0_0'$time/0.1),100)]/10, col="black", cex=2)
lines(all_Sseries$'0_0.5'$avg, all_Reffs$'0_0.5'$avg, col="orange", lty="solid",lwd=2)
text(all_Sseries$'0_0.5'$avg[seq(1,max(all_Sseries$'0_0.5'$time/0.1),100)], all_Reffs$'0_0.5'$avg[seq(1,max(all_Reffs$'0_0.5'$time/0.1),100)], labels=all_Reffs$'0_0.5'$time[seq(1,max(all_Reffs$'0_0.5'$time/0.1),100)]/10, col="orange", cex=2)
#polygon(c(rev(all_Sseries$'0_0.5'$avg), all_Sseries$'0_0.5'$avg), c(rev(all_Reffs$'0_0.5'$low), all_Reffs$'0_0.5'$high),
#        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'1_0'$avg, all_Reffs$'1_0'$avg, col="orchid", lty="solid",lwd=2)
text(all_Sseries$'1_0'$avg[seq(1,max(all_Sseries$'1_0'$time/0.1),100)], all_Reffs$'1_0'$avg[seq(1,max(all_Reffs$'1_0'$time/0.1),100)], labels=all_Reffs$'1_0'$time[seq(1,max(all_Reffs$'1_0'$time/0.1),100)]/10, col="orchid", cex=2)
#polygon(c(rev(all_Sseries$'1_0'$avg), all_Sseries$'1_0'$avg), c(rev(all_Reffs$'1_0'$low), all_Reffs$'1_0'$high),
#        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'1_0.5_0'$avg, all_Reffs$'1_0.5_0'$avg, col="yellow3", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'1_0.5_0'$avg), all_Sseries$'1_0.5_0'$avg), c(rev(all_Reffs$'1_0.5_0'$low), all_Reffs$'1_0.5_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
text(all_Sseries$'1_0.5_0'$avg[seq(1,max(all_Sseries$'1_0.5_0'$time/0.1),100)], all_Reffs$'1_0.5_0'$avg[seq(1,max(all_Reffs$'1_0.5_0'$time/0.1),100)], all_Reffs$'1_0.5_0'$time[seq(1,max(all_Reffs$'1_0.5_0'$time/0.1),100)]/10, col="yellow3", cex=2)
lines(all_Sseries$'1_0.5_-1'$avg, all_Reffs$'1_0.5_-1'$avg, col="red",lwd=2)
#polygon(c(rev(all_Sseries$'1_0.5_-1'$avg), all_Sseries$'1_0.5_-1'$avg), c(rev(all_Reffs$'1_0.5_-1'$low), all_Reffs$'1_0.5_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
text(all_Sseries$'1_0.5_-1'$avg[seq(1,max(all_Sseries$'1_0.5_-1'$time/0.1),100)], all_Reffs$'1_0.5_-1'$avg[seq(1,max(all_Reffs$'1_0.5_-1'$time/0.1),100)], labels=all_Reffs$'1_0.5_-1'$time[seq(1,max(all_Reffs$'1_0.5_-1'$time/0.1),100)]/10, col="red", cex=2)
lines(all_Sseries$'1_0.5_1'$avg, all_Reffs$'1_0.5_1'$avg, type="l", col="blue", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'1_0.5_1'$avg), all_Sseries$'1_0.5_1'$avg), c(rev(all_Reffs$'1_0.5_1'$low), all_Reffs$'1_0.5_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
text(all_Sseries$'1_0.5_1'$avg[seq(1,max(all_Sseries$'1_0.5_1'$time/0.1),100)], all_Reffs$'1_0.5_1'$avg[seq(1,max(all_Reffs$'1_0.5_1'$time/0.1),100)], labels=all_Reffs$'1_0.5_1'$time[seq(1,max(all_Reffs$'1_0.5_1'$time/0.1),100)]/10, col="blue", cex=2)
abline(h=1, col="gray75", lty="dashed")

plot(all_Sseries$'0_0'$avg, all_Reffs$'0_0'$avg, type="l", col="black",lwd=2,
     xlab="Number of susceptible individuals", ylab="", xlim=c(1000,0), ylim=c(0,10), cex.axis=1.3, cex.lab=1.5)
#polygon(c(rev(all_Sseries$'0_0'$avg), all_Sseries$'0_0'$avg), c(rev(all_Reffs$'0_0'$low), all_Reffs$'0_0'$high),
#        col=adjustcolor("black", alpha.f=0.2) , lty = 0)
text(all_Sseries$'0_0'$avg[seq(1,max(all_Sseries$'0_0'$time/0.1),100)], all_Reffs$'0_0'$avg[seq(1,max(all_Reffs$'0_0'$time/0.1),100)], labels=all_Reffs$'0_0'$time[seq(1,max(all_Reffs$'0_0'$time/0.1),100)]/10, col="black", cex=2)
lines(all_Sseries$'0_0.5'$avg, all_Reffs$'0_0.5'$avg, col="orange", lty="solid",lwd=2)
text(all_Sseries$'0_0.5'$avg[seq(1,max(all_Sseries$'0_0.5'$time/0.1),100)], all_Reffs$'0_0.5'$avg[seq(1,max(all_Reffs$'0_0.5'$time/0.1),100)], labels=all_Reffs$'0_0.5'$time[seq(1,max(all_Reffs$'0_0.5'$time/0.1),100)]/10, col="orange", cex=2)
#polygon(c(rev(all_Sseries$'0_0.5'$avg), all_Sseries$'0_0.5'$avg), c(rev(all_Reffs$'0_0.5'$low), all_Reffs$'0_0.5'$high),
#        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_0'$avg, all_Reffs$'3_0'$avg, col="orchid", lty="solid",lwd=2)
text(all_Sseries$'3_0'$avg[seq(1,max(all_Sseries$'3_0'$time/0.1),100)], all_Reffs$'3_0'$avg[seq(1,max(all_Reffs$'3_0'$time/0.1),100)], labels=all_Reffs$'3_0'$time[seq(1,max(all_Reffs$'3_0'$time/0.1),100)]/10, col="orchid", cex=2)
#polygon(c(rev(all_Sseries$'3_0'$avg), all_Sseries$'3_0'$avg), c(rev(all_Reffs$'3_0'$low), all_Reffs$'3_0'$high),
#        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_0.5_0'$avg, all_Reffs$'3_0.5_0'$avg, col="yellow3", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'3_0.5_0'$avg), all_Sseries$'3_0.5_0'$avg), c(rev(all_Reffs$'3_0.5_0'$low), all_Reffs$'3_0.5_0'$high),
#        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
text(all_Sseries$'3_0.5_0'$avg[seq(1,max(all_Sseries$'3_0.5_0'$time/0.1),100)], all_Reffs$'3_0.5_0'$avg[seq(1,max(all_Reffs$'3_0.5_0'$time/0.1),100)], labels=all_Reffs$'3_0.5_0'$time[seq(1,max(all_Reffs$'3_0.5_0'$time/0.1),100)]/10, col="yellow3", cex=2)
lines(all_Sseries$'3_0.5_-1'$avg, all_Reffs$'3_0.5_-1'$avg, col="red",lwd=2)
#polygon(c(rev(all_Sseries$'3_0.5_-1'$avg), all_Sseries$'3_0.5_-1'$avg), c(rev(all_Reffs$'3_0.5_-1'$low), all_Reffs$'3_0.5_-1'$high),
#        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
text(all_Sseries$'3_0.5_-1'$avg[seq(1,max(all_Sseries$'3_0.5_-1'$time/0.1),100)], all_Reffs$'3_0.5_-1'$avg[seq(1,max(all_Reffs$'3_0.5_-1'$time/0.1),100)], labels=all_Reffs$'3_0.5_-1'$time[seq(1,max(all_Reffs$'3_0.5_-1'$time/0.1),100)]/10, col="red", cex=2)
lines(all_Sseries$'3_0.5_1'$avg, all_Reffs$'3_0.5_1'$avg, type="l", col="blue", lty="solid",lwd=2)
#polygon(c(rev(all_Sseries$'3_0.5_1'$avg), all_Sseries$'3_0.5_1'$avg), c(rev(all_Reffs$'3_0.5_1'$low), all_Reffs$'3_0.5_1'$high),
#        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
text(all_Sseries$'3_0.5_1'$avg[seq(1,max(all_Sseries$'3_0.5_1'$time/0.1),100)], all_Reffs$'3_0.5_1'$avg[seq(1,max(all_Reffs$'3_0.5_1'$time/0.1),100)], labels=all_Reffs$'3_0.5_1'$time[seq(1,max(all_Reffs$'3_0.5_1'$time/0.1),100)]/10, col="blue", cex=2)
abline(h=1, col="gray75", lty="dashed")


par(mai = c(0,0,0,0))

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top", bty="n", inset=0, horiz=T, legend=c("Hom","HiT only","HiS only","Neg cor","No cor","Pos cor"),
       col=c("black","orange","orchid","red","yellow3","blue"),
       lty="solid", lwd=2, cex=1.5)
#legend("top", bty="n", inset=0, horiz=T, legend=c("Hom","Neg cor","No cor","Pos cor"),
#       col=c("black","red","yellow3","blue"),
#       lty="solid", lwd=2, cex=1.5)

# 10 x 9



############################################################################
#### plot time to jth infection for j=50

# first, collect avg time to 50th inf from each parameter combo

# initialize data frame to hold all info
time50thinf <- data.frame(cv_s = double(), 
                          cv_t = double(), 
                          cor = character(), 
                          time = double(),
                          stringsAsFactors = FALSE) 

# extract info from dataframes that have heterogeneity
i <- 1  # initalize counter variable
for(cv_s in c(0.5, 1, 3))
{
  for(cv_t in c(0.5, 1, 3))
  {
    for(corr in c(-1, -0.5, 0, 0.5, 1))
    {
      if(50 %in% all_timestojthinf[[i]]$jinf) # check if have at least 50 indivs infected
        time_j50 <- all_timestojthinf[[i]][which(all_timestojthinf[[i]]$jinf == 50),]$avg  # time for 50th inf
      else # there aren't 50 indivs infected
        time_j50 <- NA
      time50thinf[i,] <- c(cv_s, cv_t, corr, time_j50)  # add data to dataframe
      i <- i + 1
    }
  }
}

# extract info from dataframes with homogeneity
cvs_hom <- c(0,0.5,1,3,0,0,0)
cvt_hom <- c(0,0,0,0,0.5,1,3)

for(j in 1:7)
{
  time_j50 <- all_timestojthinf[[i]][which(all_timestojthinf[[i]]$jinf == 50),]$avg  # time for 50th inf
  time50thinf[i,] <- c(cvs_hom[j], cvt_hom[j], 0, time_j50)  # add data to dataframe
  i <- i + 1
}


# set factors
time50thinf$cv_s <- as.factor(time50thinf$cv_s)
time50thinf$cv_t <- as.factor(time50thinf$cv_t)
time50thinf$cor <- factor(time50thinf$cor, levels = c(-1,-0.5,0,0.5,1))



# plot heat map of avg time to 50th inf
timej_corn1 <- ggplot(time50thinf[which(time50thinf$cor == -1),], aes(x=cv_s, y=cv_t, fill=time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(10,122),  breaks=seq(10,122,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to jth Infection") +
  xlab("") +
  ylab("") +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), plot.margin=unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")
timej_corn1 <- set_panel_size( timej_corn1, height=unit(3*2, "cm"),
                               width=unit(3*2, "cm") )

timej_corn0.5 <- ggplot(time50thinf[which(time50thinf$cor == -0.5),], aes(x=cv_s, y=cv_t, fill=time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(10,122),  breaks=seq(10,122,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to jth Infection") +
  xlab("") +
  ylab("") +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), plot.margin=unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")
timej_corn0.5 <- set_panel_size( timej_corn0.5, height=unit(3*2, "cm"),
                                 width=unit(3*2, "cm") )

timej_cor0 <- ggplot(time50thinf[which(time50thinf$cor == 0 & time50thinf$cv_s != 0 & time50thinf$cv_t != 0),], aes(x=cv_s, y=cv_t, fill=time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(10,122),  breaks=seq(10,122,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to jth Infection") +
  xlab("") +
  ylab("") +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), plot.margin=unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")
timej_cor0 <- set_panel_size( timej_cor0, height=unit(3*2, "cm"),
                              width=unit(3*2, "cm") )

timej_cor0.5 <- ggplot(time50thinf[which(time50thinf$cor == 0.5),], aes(x=cv_s, y=cv_t, fill=time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(10,122),  breaks=seq(10,122,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to jth Infection") +
  xlab("") +
  ylab("") +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), plot.margin=unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")
timej_cor0.5 <- set_panel_size( timej_cor0.5, height=unit(3*2, "cm"),
                                width=unit(3*2, "cm") )

timej_cor1 <- ggplot(time50thinf[which(time50thinf$cor == 1),], aes(x=cv_s, y=cv_t, fill=time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(10,122),  breaks=seq(10,122,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to jth Infection") +
  xlab("") +
  ylab("") +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), plot.margin=unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")
timej_cor1 <- set_panel_size( timej_cor1, height=unit(3*2, "cm"),
                              width=unit(3*2, "cm") )

timej_cs0 <- ggplot(time50thinf[which(time50thinf$cv_s == 0 & time50thinf$cv_t != 0),], aes(x=cor, y=cv_t, fill=time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(10,122),  breaks=seq(10,122,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to jth Infection") +
  xlab("") +
  ylab("") +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), plot.margin=unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")
timej_cs0 <- set_panel_size( timej_cs0, height=unit(3*2, "cm"),
                             width=unit(1*2, "cm") )

timej_cs0ct0 <- ggplot(time50thinf[which(time50thinf$cv_s == 0 & time50thinf$cv_t == 0),], aes(x=cor, y=cv_t, fill=time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(10,122),  breaks=seq(10,122,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to jth Infection") +
  xlab("") +
  ylab("") +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), plot.margin=unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")
timej_cs0ct0 <- set_panel_size( timej_cs0ct0, height=unit(1*2, "cm"),
                                width=unit(1*2, "cm") )


timej_ct0 <- ggplot(time50thinf[which(time50thinf$cv_s != 0 & time50thinf$cv_t == 0),], aes(x=cv_s, y=cv_t, fill=time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(10,122),  breaks=seq(10,122,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to jth Infection") +
  xlab("") +
  ylab("") +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), plot.margin=unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")
timej_ct0 <- set_panel_size( timej_ct0, height=unit(1*2, "cm"),
                             width=unit(3*2, "cm") )


# determine legend for combined plots
timej_leg <- ggplot(time50thinf[which(time50thinf$cor == "0"),], aes(x=cv_t, y=cv_s, fill=time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  #scale_fill_steps(low="white", high="red", limits=c(4,67),  breaks=seq(4,67,7),  # R0=3
  scale_fill_steps(low="white", high="red", limits=c(10,122),  breaks=seq(10,122,8),  # R0=0.8
  #scale_fill_steps(low="white", high="red", limits=c(9,84),  breaks=seq(9,84,5),  # R0=1.1
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to 50th Infection") +
  xlab(expression(paste("Coefficient of variation of transmission (",italic(C)[italic(t)],")"))) +
  ylab(expression(paste("Coefficient of variation of susceptibility ( ",italic(C)[italic(s)],")"))) +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), legend.key.height=unit(2.5,"cm"))
legend <- get_legend(timej_leg)

# plot all combined plots in a grid with common legend
plt <- plot_grid(timej_cs0, timej_corn1, timej_corn0.5, timej_cor0, timej_cor0.5, timej_cor1,
                 NULL, NULL, NULL, NULL, NULL, NULL,
                 timej_cs0ct0, NULL, NULL, timej_ct0, NULL, NULL,
                 nrow=3, ncol=6, rel_heights=c(1,-0.5,1), rel_widths=c(0.4,1,1,1,1,1))
plot_grid(plt, legend,
          nrow=1, rel_widths=c(1,0.2))
# size: 19 x 7.5


