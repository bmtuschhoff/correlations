# read in all files for individual-based SIR model for each parameter combination for the level of heterogeneity in susceptibility (cv_s), heterogeneity in transmission (cv_t), and correlation (cor)
# determine mean, median, 95% CI for SIR dynamics given a major outbreak and use bootstrapping to get variance for probability of a major epidemic
# plot SIR dynamics

library(dplyr)
library(ggplot2)
library(metR)

# number of simulations
sims <- 200

# population size for each simulation
pop_size <- 1000

# set working directory to get data from
setwd("C:/Users/bmt5507/Documents/Cor_HiST_Copula_I1_R0.8")

# function to read in datasets that are tab delimited
read_sav <- function(fileName) {read.table(fileName, header=T, sep="\t")}


# function to get I over epidemic adjusted to longest time series
get_Itimeseries <- function(epi_sim, cv_s, max_length)
{
  Iseries <- rep(0, max_length)
  
  # determine whether there was a major outbreak
  # set threshold for major outbreak based on level of HiS because more HiS = smaller outbreak but not less likely necessarily
  if(cv_s == 0 || cv_s == 0.5)
    threshold <- 200
  if(cv_s == 1)
    threshold <- 100
  if(cv_s == 3)
    threshold <- 50
  
  # check if final epidemic size (FES) is larger than threshold for major epidemic
  FES <- epi_sim$I[nrow(epi_sim)] + epi_sim$R[nrow(epi_sim)]
  if(FES > threshold)  # major epidemic, set number inf for each time step
  {
    i <- 1
    for(time in seq(0, max(epi_sim$time)+0.1, 0.1))
    {
      lessthan_time <- which(epi_sim$time <= time)
      index_time <- lessthan_time[length(lessthan_time)]
      Iseries[i] <- epi_sim[index_time,]$I
      i <- i + 1
    }  
    
    return(Iseries)
  }
  else  # no major epidemic, return all NAs
    return(rep(NA, max_length))
}


# function to get S over epidemic adjusted to longest time series
get_Stimeseries <- function(epi_sim, cv_s, max_length)
{
  Sseries <- rep(epi_sim$S[nrow(epi_sim)], max_length)
  
  # determine whether there was a major outbreak
  # set threshold for major outbreak based on level of HiS because more HiS = smaller outbreak but not less likely necessarily
  if(cv_s == 0 || cv_s == 0.5)
    threshold <- 200
  if(cv_s == 1)
    threshold <- 100
  if(cv_s == 3)
    threshold <- 50
  
  # check if final epidemic size (FES) is larger than threshold for major epidemic
  FES <- epi_sim$I[nrow(epi_sim)] + epi_sim$R[nrow(epi_sim)]
  if(FES > threshold)  # major epidemic, set number inf for each time step
  {
    i <- 1
    for(time in seq(0, max(epi_sim$time)+0.1, 0.1))
    {
      lessthan_time <- which(epi_sim$time <= time)
      index_time <- lessthan_time[length(lessthan_time)]
      Sseries[i] <- epi_sim[index_time,]$S
      i <- i + 1
    }  
    
    return(Sseries)
  }
  else  # no major epidemic, return all NAs
    return(rep(NA, max_length))
}


# levels of heterogeneity in transmission, susceptibility, and correlation
cv_ts <- c(0.5, 1, 3)
cv_ss <- c(0.5, 1, 3)
cors <- c(-1,-0.5,0,0.5,1)

# make empty lists to store all data  
all_Iseries <- vector("list", length=length(cv_ts)*length(cv_ss)*length(cors))
all_Sseries <- vector("list", length=length(cv_ts)*length(cv_ss)*length(cors))


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
      max_length <- length(seq(0,max_time+0.1,0.1))
      
      ## time series for I indivs
      # make data frame containing series of I indivs for each time series of the sims
      temp <- data.frame(sapply(read.all, get_Itimeseries, cv_s=cv_s, max_length=max_length))
      
      # make sure all numeric
      temp <- data.frame(sapply(temp, as.numeric))
      
      # find mean, 95% CIs for each time point
      I_mean <- apply(temp, 1, mean, na.rm=T)
      I_low <- apply(temp, 1, function(x){quantile(x, probs=0.025, na.rm=T)})
      I_high <- apply(temp, 1, function(x){quantile(x, probs=0.975, na.rm=T)})
      
      # make data frame containing time, mean, and 95% CIs
      element_name <- paste(cv_s,cv_t,corr, sep="_")
      Iseries <- data.frame(time=seq(0,max_time+0.1,0.1),
                               avg=I_mean,
                               low=I_low,
                               high=I_high)
      
      # add data frame to list of data frames for each param combo (cv_s, cv_t, cor) and name
      all_Iseries[[combo]] <- Iseries
      names(all_Iseries)[[combo]] <- element_name
      
      
      ## time series for S indivs
      # make data frame containing series of I indivs for each time series of the sims
      temp <- data.frame(sapply(read.all, get_Stimeseries, cv_s=cv_s, max_length=max_length))
      
      # make sure all numeric
      temp <- data.frame(sapply(temp, as.numeric))
      
      # find mean, 95% CIs for each time point
      S_mean <- apply(temp, 1, mean, na.rm=T)
      S_low <- apply(temp, 1, function(x){quantile(x, probs=0.025, na.rm=T)})
      S_high <- apply(temp, 1, function(x){quantile(x, probs=0.975, na.rm=T)})
      
      # make data frame containing time, mean, and 95% CIs
      element_name <- paste(cv_s,cv_t,corr, sep="_")
      Sseries <- data.frame(time=seq(0,max_time+0.1,0.1),
                            avg=S_mean,
                            low=S_low,
                            high=S_high)
      
      # add data frame to list of data frames for each param combo (cv_s, cv_t, cor) and name
      all_Sseries[[combo]] <- Sseries
      names(all_Sseries)[[combo]] <- element_name
      
      combo <- combo + 1
    }
  }
}


## separately get data for homogeneous case (HiS=0, HiT=0) and cases with just HiS or just HiT
# list of file names
filePatterns <- c()
filePatterns[1] <- paste("results_his0hit0cor0g0.1p1000R00.8initI1s*", sep="")
filePatterns[2] <- paste("results_his0.5hit0cor0g0.1p1000R00.8initI1s*", sep="")
filePatterns[3] <- paste("results_his1hit0cor0g0.1p1000R00.8initI1s*", sep="")
filePatterns[4] <- paste("results_his3hit0cor0g0.1p1000R00.8initI1s*", sep="")
filePatterns[5] <- paste("results_his0hit0.5cor0g0.1p1000R00.8initI1s*", sep="")
filePatterns[6] <- paste("results_his0hit1cor0g0.1p1000R00.8initI1s*", sep="")
filePatterns[7] <- paste("results_his0hit3cor0g0.1p1000R00.8initI1s*", sep="")

# levels of heterogeneity
cv_ss <- c(0,0.5,1,3,0,0,0)
cv_ts <- c(0,0,0,0,0.5,1,3)

# get data for each level of HiS, HiT, and correlation
for(i in 1:7)
{
  files <- list.files(pattern=filePatterns[i])
  
  # make a list that contains all files
  read.all <- lapply(files, read_sav)
  
  # determine length of longest time series based on time
  max_time <- max(unlist(lapply(read.all, function(x) {x$time[nrow(x)]})))
  max_length <- length(seq(0,max_time+0.1,0.1))
  
  ## time series for I indivs
  # make data frame containing series of I indivs for each time series of the sims
  temp <- data.frame(sapply(read.all, get_Itimeseries, cv_s=cv_ss[i], max_length=max_length))
  
  # make sure all numeric
  temp <- data.frame(sapply(temp, as.numeric))
  
  # # plot all I lines
  # for(column in 1:100)
  # {
  #   if(!is.na(temp[1,column]))
  #   {
  #     lines(seq(0,max_time+0.1,0.1), temp[,column], col="orange", lwd=2)
  #   }
  # }
  
  # find mean, 95% CIs for each time point
  I_mean <- apply(temp, 1, mean, na.rm=T)
  I_low <- apply(temp, 1, function(x){quantile(x, probs=0.025, na.rm=T)})
  I_high <- apply(temp, 1, function(x){quantile(x, probs=0.975, na.rm=T)})
  
  # make data frame containing time, mean, and 95% CIs
  element_name <- paste(cv_ss[i],cv_ts[i], sep="_")
  Iseries <- data.frame(time=seq(0,max_time+0.1,0.1),
                           avg=I_mean,
                           low=I_low,
                           high=I_high)
  
  # add data frame to list of data frames for each param combo (cv_s, cv_t, cor) and name
  all_Iseries[[combo]] <- Iseries
  names(all_Iseries)[[combo]] <- element_name
  
  
  
  ## time series for S indivs
  # make data frame containing series of I indivs for each time series of the sims
  temp <- data.frame(sapply(read.all, get_Stimeseries, cv_s=cv_ss[i], max_length=max_length))
  
  # make sure all numeric
  temp <- data.frame(sapply(temp, as.numeric))
  
  # find mean, 95% CIs for each time point
  S_mean <- apply(temp, 1, mean, na.rm=T)
  S_low <- apply(temp, 1, function(x){quantile(x, probs=0.025, na.rm=T)})
  S_high <- apply(temp, 1, function(x){quantile(x, probs=0.975, na.rm=T)})
  
  # make data frame containing time, mean, and 95% CIs
  element_name <- paste(cv_ss[i],cv_ts[i], sep="_")
  Sseries <- data.frame(time=seq(0,max_time+0.1,0.1),
                        avg=S_mean,
                        low=S_low,
                        high=S_high)
  
  # add data frame to list of data frames for each param combo (cv_s, cv_t, cor) and name
  all_Sseries[[combo]] <- Sseries
  names(all_Sseries)[[combo]] <- element_name
  
  combo <- combo + 1
}



#########################################################################
########################## plot dynamics ################################
#########################################################################
par(mar=c(5.1,5.1,4.1,2.1)) # change margin sizes

###################### plot I over time #################################
# 8.8" width x 6.7" height
## HiS=0.5, HiT=0.5
# hom, HiS, HiT
plot(all_Iseries$'0_0'$time, all_Iseries$'0_0'$avg, type="l", col="black", lwd=2, 
     xlab="Time", ylab="Number of infected individuals", xlim=c(0,200), ylim=c(0,350), cex.axis=2.2, cex.lab=2.2)
polygon(c(rev(all_Iseries$'0_0'$time), all_Iseries$'0_0'$time), c(rev(all_Iseries$'0_0'$low), all_Iseries$'0_0'$high),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'0.5_0'$time, all_Iseries$'0.5_0'$avg, type="l", col="orchid", lwd=2)
polygon(c(rev(all_Iseries$'0.5_0'$time), all_Iseries$'0.5_0'$time), c(rev(all_Iseries$'0.5_0'$low), all_Iseries$'0.5_0'$high),
        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'0_0.5'$time, all_Iseries$'0_0.5'$avg, type="l", col="orange", lwd=2)
polygon(c(rev(all_Iseries$'0_0.5'$time), all_Iseries$'0_0.5'$time), c(rev(all_Iseries$'0_0.5'$low), all_Iseries$'0_0.5'$high),
        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)
legend("topright", bty="n", legend=c("Hom","HiS only","HiT only"),
       col=c("black","orchid","orange"),
       lty=c("solid","solid","solid"), lwd=2, cex=1.7)


# hom, no cor, neg, pos
plot(all_Iseries$'0_0'$time, all_Iseries$'0_0'$avg, type="l", col="black", lwd=2, 
     xlab="Time", ylab="Number of infected individuals", xlim=c(0,200), ylim=c(0,350), cex.axis=2.2, cex.lab=2.2)
polygon(c(rev(all_Iseries$'0_0'$time), all_Iseries$'0_0'$time), c(rev(all_Iseries$'0_0'$low), all_Iseries$'0_0'$high),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'0.5_0.5_0'$time, all_Iseries$'0.5_0.5_0'$avg, type="l", col="yellow3", lwd=2)
polygon(c(rev(all_Iseries$'0.5_0.5_0'$time), all_Iseries$'0.5_0.5_0'$time), c(rev(all_Iseries$'0.5_0.5_0'$low), all_Iseries$'0.5_0.5_0'$high),
        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'0.5_0.5_1'$time, all_Iseries$'0.5_0.5_1'$avg, type="l", col="blue", lwd=2)
polygon(c(rev(all_Iseries$'0.5_0.5_1'$time), all_Iseries$'0.5_0.5_1'$time), c(rev(all_Iseries$'0.5_0.5_1'$low), all_Iseries$'0.5_0.5_1'$high),
        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'0.5_0.5_-1'$time, all_Iseries$'0.5_0.5_-1'$avg, type="l", col="red", lwd=2)
polygon(c(rev(all_Iseries$'0.5_0.5_-1'$time), all_Iseries$'0.5_0.5_-1'$time), c(rev(all_Iseries$'0.5_0.5_-1'$low), all_Iseries$'0.5_0.5_-1'$high),
        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
legend("topright", bty="n", legend=c("Hom","Neg cor","No cor","Pos cor"),
       col=c("black","red","yellow3","blue"),
       lty=c("solid","solid","solid","solid"), lwd=2, cex=1.7)


## HiS=3, HiT=3
# hom, HiS, HiT
plot(all_Iseries$'0_0'$time, all_Iseries$'0_0'$avg, type="l", col="black", lwd=2, 
     xlab="Time", ylab="Number of infected individuals", xlim=c(0,200), ylim=c(0,350), cex.axis=2.2, cex.lab=2.2)
polygon(c(rev(all_Iseries$'0_0'$time), all_Iseries$'0_0'$time), c(rev(all_Iseries$'0_0'$low), all_Iseries$'0_0'$high),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'3_0'$time, all_Iseries$'3_0'$avg, type="l", col="orchid", lwd=2)
polygon(c(rev(all_Iseries$'3_0'$time), all_Iseries$'3_0'$time), c(rev(all_Iseries$'3_0'$low), all_Iseries$'3_0'$high),
        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'0_3'$time, all_Iseries$'0_3'$avg, type="l", col="orange", lwd=2)
polygon(c(rev(all_Iseries$'0_3'$time), all_Iseries$'0_3'$time), c(rev(all_Iseries$'0_3'$low), all_Iseries$'0_3'$high),
        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)
legend("topright", bty="n", legend=c("Hom","HiS only","HiT only"),
       col=c("black","orchid","orange"),
       lty=c("solid","solid","solid"), lwd=2, cex=1.7)

# hom, no cor, neg, pos
plot(all_Iseries$'0_0'$time, all_Iseries$'0_0'$avg, type="l", col="black", lwd=2, 
     xlab="Time", ylab="Number of infected individuals", xlim=c(0,200), ylim=c(0,350), cex.axis=2.2, cex.lab=2.2)
polygon(c(rev(all_Iseries$'0_0'$time), all_Iseries$'0_0'$time), c(rev(all_Iseries$'0_0'$low), all_Iseries$'0_0'$high),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'3_3_0'$time, all_Iseries$'3_3_0'$avg, type="l", col="yellow3", lwd=2)
polygon(c(rev(all_Iseries$'3_3_0'$time), all_Iseries$'3_3_0'$time), c(rev(all_Iseries$'3_3_0'$low), all_Iseries$'3_3_0'$high),
        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'3_3_1'$time, all_Iseries$'3_3_1'$avg, type="l", col="blue", lwd=2)
polygon(c(rev(all_Iseries$'3_3_1'$time), all_Iseries$'3_3_1'$time), c(rev(all_Iseries$'3_3_1'$low), all_Iseries$'3_3_1'$high),
        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'3_3_-1'$time, all_Iseries$'3_3_-1'$avg, type="l", col="red", lwd=2)
polygon(c(rev(all_Iseries$'3_3_-1'$time), all_Iseries$'3_3_-1'$time), c(rev(all_Iseries$'3_3_-1'$low), all_Iseries$'3_3_-1'$high),
        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
legend("topright", bty="n", legend=c("Hom","Neg cor","No cor","Pos cor"),
       col=c("black","red","yellow3","blue"),
       lty=c("solid","solid","solid","solid"), lwd=2, cex=1.7)


########### plot I for R0=0.8, 1.1, combined plots for just HiS=HiT=3 ###############
plot(all_Iseries$'0_0'$time, all_Iseries$'0_0'$avg, type="l", col="black", lwd=2, 
     xlab="Time", ylab="Number of infected individuals", xlim=c(0,400), ylim=c(0,350), cex.axis=2.2, cex.lab=2.2)
polygon(c(rev(all_Iseries$'0_0'$time), all_Iseries$'0_0'$time), c(rev(all_Iseries$'0_0'$low), all_Iseries$'0_0'$high),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'3_0'$time, all_Iseries$'3_0'$avg, type="l", col="orchid", lwd=2)
polygon(c(rev(all_Iseries$'3_0'$time), all_Iseries$'3_0'$time), c(rev(all_Iseries$'3_0'$low), all_Iseries$'3_0'$high),
        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'0_3'$time, all_Iseries$'0_3'$avg, type="l", col="orange", lwd=2)
polygon(c(rev(all_Iseries$'0_3'$time), all_Iseries$'0_3'$time), c(rev(all_Iseries$'0_3'$low), all_Iseries$'0_3'$high),
        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'3_3_0'$time, all_Iseries$'3_3_0'$avg, type="l", col="yellow3", lwd=2)
polygon(c(rev(all_Iseries$'3_3_0'$time), all_Iseries$'3_3_0'$time), c(rev(all_Iseries$'3_3_0'$low), all_Iseries$'3_3_0'$high),
        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'3_3_1'$time, all_Iseries$'3_3_1'$avg, type="l", col="blue", lwd=2)
polygon(c(rev(all_Iseries$'3_3_1'$time), all_Iseries$'3_3_1'$time), c(rev(all_Iseries$'3_3_1'$low), all_Iseries$'3_3_1'$high),
        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Iseries$'3_3_-1'$time, all_Iseries$'3_3_-1'$avg, type="l", col="red", lwd=2)
polygon(c(rev(all_Iseries$'3_3_-1'$time), all_Iseries$'3_3_-1'$time), c(rev(all_Iseries$'3_3_-1'$low), all_Iseries$'3_3_-1'$high),
        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
legend("topright", bty="n", legend=c("Hom", "HiS only", "HiT only","Neg cor","No cor","Pos cor"),
       col=c("black","orchid","orange","red","yellow3","blue"),
       lty=c("solid","solid","solid","solid","solid","solid"), lwd=2, cex=1.7)


###################### plot S over time #################################
# 8.8 x 7.4
## HiS=0.5, HiT=0.5
# hom, HiS, HiT
plot(all_Sseries$'0_0'$time, pop_size-all_Sseries$'0_0'$avg, type="l", col="black", lwd=2, 
     xlab="Time", ylab="Cumulative number infected", xlim=c(0,200), ylim=c(0,pop_size), cex.axis=2.2, cex.lab=2.2)
polygon(c(rev(all_Sseries$'0_0'$time), all_Sseries$'0_0'$time), c(rev(pop_size-all_Sseries$'0_0'$low), pop_size-all_Sseries$'0_0'$high),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'0.5_0'$time, pop_size-all_Sseries$'0.5_0'$avg, type="l", col="orchid", lwd=2)
polygon(c(rev(all_Sseries$'0.5_0'$time), all_Sseries$'0.5_0'$time), c(rev(pop_size-all_Sseries$'0.5_0'$low), pop_size-all_Sseries$'0.5_0'$high),
        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'0_0.5'$time, pop_size-all_Sseries$'0_0.5'$avg, type="l", col="orange", lwd=2)
polygon(c(rev(all_Sseries$'0_0.5'$time), all_Sseries$'0_0.5'$time), c(rev(pop_size-all_Sseries$'0_0.5'$low), pop_size-all_Sseries$'0_0.5'$high),
        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)

# hom, no cor, neg, pos
plot(all_Sseries$'0_0'$time, pop_size-all_Sseries$'0_0'$avg, type="l", col="black", lwd=2, 
     xlab="Time", ylab="Cumulative number infected", xlim=c(0,200), ylim=c(0,pop_size), cex.axis=2.2, cex.lab=2.2)
polygon(c(rev(all_Sseries$'0_0'$time), all_Sseries$'0_0'$time), c(rev(pop_size-all_Sseries$'0_0'$low), pop_size-all_Sseries$'0_0'$high),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'0.5_0.5_0'$time, pop_size-all_Sseries$'0.5_0.5_0'$avg, type="l", col="yellow3", lwd=2)
polygon(c(rev(all_Sseries$'0.5_0.5_0'$time), all_Sseries$'0.5_0.5_0'$time), c(rev(pop_size-all_Sseries$'0.5_0.5_0'$low), pop_size-all_Sseries$'0.5_0.5_0'$high),
        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'0.5_0.5_1'$time, pop_size-all_Sseries$'0.5_0.5_1'$avg, type="l", col="blue", lwd=2)
polygon(c(rev(all_Sseries$'0.5_0.5_1'$time), all_Sseries$'0.5_0.5_1'$time), c(rev(pop_size-all_Sseries$'0.5_0.5_1'$low), pop_size-all_Sseries$'0.5_0.5_1'$high),
        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'0.5_0.5_-1'$time, pop_size-all_Sseries$'0.5_0.5_-1'$avg, type="l", col="red", lwd=2)
polygon(c(rev(all_Sseries$'0.5_0.5_-1'$time), all_Sseries$'0.5_0.5_-1'$time), c(rev(pop_size-all_Sseries$'0.5_0.5_-1'$low), pop_size-all_Sseries$'0.5_0.5_-1'$high),
        col=adjustcolor("red", alpha.f=0.2) , lty = 0)

## HiS=3, HiT=3
# hom, HiS, HiT
plot(all_Sseries$'0_0'$time, pop_size-all_Sseries$'0_0'$avg, type="l", col="black", lwd=2, 
     xlab="Time", ylab="Cumulative number infected", xlim=c(0,200), ylim=c(0,pop_size), cex.axis=2.2, cex.lab=2.2)
polygon(c(rev(all_Sseries$'0_0'$time), all_Sseries$'0_0'$time), c(rev(pop_size-all_Sseries$'0_0'$low), pop_size-all_Sseries$'0_0'$high),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_0'$time, pop_size-all_Sseries$'3_0'$avg, type="l", col="orchid", lwd=2)
polygon(c(rev(all_Sseries$'3_0'$time), all_Sseries$'3_0'$time), c(rev(pop_size-all_Sseries$'3_0'$low), pop_size-all_Sseries$'3_0'$high),
        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'0_3'$time, pop_size-all_Sseries$'0_3'$avg, type="l", col="orange", lwd=2)
polygon(c(rev(all_Sseries$'0_3'$time), all_Sseries$'0_3'$time), c(rev(pop_size-all_Sseries$'0_3'$low), pop_size-all_Sseries$'0_3'$high),
        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)

# hom, no cor, neg, pos
plot(all_Sseries$'0_0'$time, pop_size-all_Sseries$'0_0'$avg, type="l", col="black", lwd=2, 
     xlab="Time", ylab="Cumulative number infected", xlim=c(0,200), ylim=c(0,pop_size), cex.axis=2.2, cex.lab=2.2)
polygon(c(rev(all_Sseries$'0_0'$time), all_Sseries$'0_0'$time), c(rev(pop_size-all_Sseries$'0_0'$low), pop_size-all_Sseries$'0_0'$high),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_3_0'$time, pop_size-all_Sseries$'3_3_0'$avg, type="l", col="yellow3", lwd=2)
polygon(c(rev(all_Sseries$'3_3_0'$time), all_Sseries$'3_3_0'$time), c(rev(pop_size-all_Sseries$'3_3_0'$low), pop_size-all_Sseries$'3_3_0'$high),
        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_3_1'$time, pop_size-all_Sseries$'3_3_1'$avg, type="l", col="blue", lwd=2)
polygon(c(rev(all_Sseries$'3_3_1'$time), all_Sseries$'3_3_1'$time), c(rev(pop_size-all_Sseries$'3_3_1'$low), pop_size-all_Sseries$'3_3_1'$high),
        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_3_-1'$time, pop_size-all_Sseries$'3_3_-1'$avg, type="l", col="red", lwd=2)
polygon(c(rev(all_Sseries$'3_3_-1'$time), all_Sseries$'3_3_-1'$time), c(rev(pop_size-all_Sseries$'3_3_-1'$low), pop_size-all_Sseries$'3_3_-1'$high),
        col=adjustcolor("red", alpha.f=0.2) , lty = 0)


########### plot S for R0=0.8, 1.1, combined plots for just HiS=HiT=3 ###############
plot(all_Sseries$'0_0'$time, pop_size-all_Sseries$'0_0'$avg, type="l", col="black", lwd=2, 
     xlab="Time", ylab="Cumulative number infected", xlim=c(0,400), ylim=c(0,pop_size), cex.axis=2.2, cex.lab=2.2)
polygon(c(rev(all_Sseries$'0_0'$time), all_Sseries$'0_0'$time), c(rev(pop_size-all_Sseries$'0_0'$low), pop_size-all_Sseries$'0_0'$high),
        col=adjustcolor("gray50", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_0'$time, pop_size-all_Sseries$'3_0'$avg, type="l", col="orchid", lwd=2)
polygon(c(rev(all_Sseries$'3_0'$time), all_Sseries$'3_0'$time), c(rev(pop_size-all_Sseries$'3_0'$low), pop_size-all_Sseries$'3_0'$high),
        col=adjustcolor("orchid", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'0_3'$time, pop_size-all_Sseries$'0_3'$avg, type="l", col="orange", lwd=2)
polygon(c(rev(all_Sseries$'0_3'$time), all_Sseries$'0_3'$time), c(rev(pop_size-all_Sseries$'0_3'$low), pop_size-all_Sseries$'0_3'$high),
        col=adjustcolor("orange", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_3_0'$time, pop_size-all_Sseries$'3_3_0'$avg, type="l", col="yellow3", lwd=2)
polygon(c(rev(all_Sseries$'3_3_0'$time), all_Sseries$'3_3_0'$time), c(rev(pop_size-all_Sseries$'3_3_0'$low), pop_size-all_Sseries$'3_3_0'$high),
        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_3_1'$time, pop_size-all_Sseries$'3_3_1'$avg, type="l", col="blue", lwd=2)
polygon(c(rev(all_Sseries$'3_3_1'$time), all_Sseries$'3_3_1'$time), c(rev(pop_size-all_Sseries$'3_3_1'$low), pop_size-all_Sseries$'3_3_1'$high),
        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Sseries$'3_3_-1'$time, pop_size-all_Sseries$'3_3_-1'$avg, type="l", col="red", lwd=2)
polygon(c(rev(all_Sseries$'3_3_-1'$time), all_Sseries$'3_3_-1'$time), c(rev(pop_size-all_Sseries$'3_3_-1'$low), pop_size-all_Sseries$'3_3_-1'$high),
        col=adjustcolor("red", alpha.f=0.2) , lty = 0)



############################################################################

############### probability of a major epidemic ############################

############################################################################
## bootstrap samples to get mean and var for prob of major epidemics
pop_size <- 1000
boot_size <- 100

# take prob of major outbreaks dataframe from HeatMap_CompProbPeakSizeTimeFES_CorHiST.R
probs_majorouts <- data_R3_probout
probs_majorouts_hom <- data_R3_hom_probout

# make empty dataframe to store this data 
bootstrap_probmajout <- data.frame(cv_s = double(), 
                                   cv_t = double(), 
                                   cor = character(), 
                                   resamp_mean = double(),
                                   resamp_sd = double(),
                                   stringsAsFactors = FALSE)

# bootstrap for each parameter combo
for(row in 1:nrow(probs_majorouts))
{
  prob <- probs_majorouts$major_out[row]
  orig_sample <- c(rep(1, prob*pop_size), rep(0, (1-prob)*pop_size))
  
  resamp <- rep(NA, 1000) # vector of all resampled probs
  
  for(i in 1:1000)
    # resample with replacement 1000 times and sum to get sampled probability
    resamp[i] <- sum(sample(orig_sample, boot_size, replace=T)) / boot_size
  
  # calculate mean and variance of resampled population
  resamp_mean <- mean(resamp)
  resamp_sd <- sd(resamp)
  
  # save data
  bootstrap_probmajout[row,] <- c(as.numeric(paste(probs_majorouts$cv_s[row])), as.numeric(paste(probs_majorouts$cv_t[row])), as.numeric(probs_majorouts$cor[row]), resamp_mean, resamp_sd)
}



# make empty dataframe to store this data 
bootstrap_probmajout_hom <- data.frame(cv_s = double(), 
                                       cv_t = double(), 
                                       cor = character(), 
                                       resamp_mean = double(),
                                       resamp_sd = double(),
                                       stringsAsFactors = FALSE)

# bootstrap for each parameter combo
for(row in 1:nrow(probs_majorouts_hom))
{
  prob <- probs_majorouts_hom$major_out[row]
  orig_sample <- c(rep(1, prob*pop_size), rep(0, (1-prob)*pop_size))
  
  resamp <- rep(NA, 1000) # vector of all resampled probs
  
  for(i in 1:1000)
    # resample with replacement 1000 times and sum to get sampled probability
    resamp[i] <- sum(sample(orig_sample, boot_size, replace=T)) / boot_size
  
  # calculate mean and variance of resampled population
  resamp_mean <- mean(resamp)
  resamp_sd <- sd(resamp)
  
  # save data
  bootstrap_probmajout_hom[row,] <- c(as.numeric(paste(probs_majorouts_hom$cv_s[row])), as.numeric(paste(probs_majorouts_hom$cv_t[row])), as.numeric(probs_majorouts_hom$cor[row]), resamp_mean, resamp_sd)
}


# combine data frames into one for C=0.5,3
bootstrap_probmajout0.5 <- rbind(bootstrap_probmajout[which(bootstrap_probmajout$cv_s == 0.5 & bootstrap_probmajout$cv_t == 0.5),],
                                 bootstrap_probmajout_hom[which(bootstrap_probmajout_hom$cv_s == 0.5 | bootstrap_probmajout_hom$cv_t == 0.5),],
                                 bootstrap_probmajout_hom[which(bootstrap_probmajout_hom$cv_s == 0 & bootstrap_probmajout_hom$cv_t == 0),])

bootstrap_probmajout3 <- rbind(bootstrap_probmajout[which(bootstrap_probmajout$cv_s == 3 & bootstrap_probmajout$cv_t == 3),],
                               bootstrap_probmajout_hom[which(bootstrap_probmajout_hom$cv_s == 3 | bootstrap_probmajout_hom$cv_t == 3),],
                               bootstrap_probmajout_hom[which(bootstrap_probmajout_hom$cv_s == 0 & bootstrap_probmajout_hom$cv_t == 0),])

# remove intermediate correlations
bootstrap_probmajout0.5 <- bootstrap_probmajout0.5[-c(2,4),]
bootstrap_probmajout3 <- bootstrap_probmajout3[-c(2,4),]

# add labels
bootstrap_probmajout0.5$name <- c("Neg \ncor", "No \ncor", "Pos \ncor", "HiT \nonly", "HiS \nonly", "Hom")
bootstrap_probmajout3$name <- c("Neg \ncor", "No \ncor", "Pos \ncor", "HiT \nonly", "HiS \nonly", "Hom")

# put labels in correct order
bootstrap_probmajout0.5$name <- factor(bootstrap_probmajout0.5$name, levels=c("Hom", "HiS \nonly", "HiT \nonly", "Neg \ncor", "No \ncor", "Pos \ncor"))
bootstrap_probmajout3$name <- factor(bootstrap_probmajout3$name, levels=c("Hom", "HiS \nonly", "HiT \nonly", "Neg \ncor", "No \ncor", "Pos \ncor"))

# plot prob of major outbreak
# size: 7.3 x 5.9
ggplot(data=bootstrap_probmajout0.5, aes(x=name, y=resamp_mean, col=name)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=resamp_mean-2*resamp_sd, ymax=resamp_mean+2*resamp_sd), size=2) +
  scale_y_continuous("Probability of major epidemic", limits=c(0,0.8)) +
  scale_x_discrete("") +
  scale_color_manual(values=c("black", "orchid", "orange", "red", "yellow3", "blue")) +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=35), axis.ticks=element_line(), legend.position="none")

ggplot(data=bootstrap_probmajout3, aes(x=name, y=resamp_mean, col=name)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=pmax(0,resamp_mean-2*resamp_sd), ymax=resamp_mean+2*resamp_sd), size=2) +
  scale_y_continuous("Probability of major epidemic", limits=c(0,0.8)) +
  scale_x_discrete("") +
  scale_color_manual(values=c("black", "orchid", "orange", "red", "yellow3", "blue")) +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=35), axis.ticks=element_line(), legend.position="none")

