# read in all files for individual-based SIR model for each parameter combination for the level of heterogeneity in susceptibility (cv_s), heterogeneity in transmission (cv_t), and correlation (cor)
# determine mean, median, 95% CI for SIR dynamics given a major outbreak and use bootstrapping to get variance for probability of a major epidemic
# plot SIR dynamics

library(dplyr)
library(ggplot2)
library(metR)
library(lubridate)

# number of simulations
sims <- 500

# population size for each simulation
pop_size <- 70180

# set working directory to get data from
setwd("C:/Users/bmt5507/Documents/Cor_HiST_mpox")

# function to read in datasets that are tab delimited
read_sav <- function(fileName) {read.table(fileName, header=T, sep="\t")}


# function to get I over epidemic adjusted to longest time series
get_Itimeseries <- function(epi_sim, cor_level, max_length)
{
  Iseries <- rep(0, max_length)
  
  # set threshold for major outbreak based on if epidemic lasts long enough to match the data
  # if correlation isn't positive, won't get a long epidemic so use all data regardless
  if(cor_level == 1 | R0 == 5.2)
    threshold <- 365
  else
    threshold <- 0
  
  # check if epidemic lasts longer than threshold time
  if(epi_sim$time[nrow(epi_sim)] >= threshold)  # major epidemic, set number inf for each time step
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

# function to get I observations (based on case detection rate) over epidemic adjusted to longest time series
get_Iobstimeseries <- function(epi_sim, cor_level, R0, max_length)
{
  Iseries <- rep(0, max_length)
  
  # set threshold for major outbreak based on if epidemic lasts long enough to match the data
  # if correlation isn't positive, won't get a long epidemic so use all data regardless
  if(cor_level == 1 | R0 == 5.2)
    threshold <- 365
  else
    threshold <- 0  
  
  # check if epidemic lasts longer than threshold time
  if(epi_sim$time[nrow(epi_sim)] >= threshold)  # major epidemic, set number inf for each time step
  {
    i <- 1
    for(time in seq(0, max(epi_sim$time)+0.1, 0.1))
    {
      lessthan_time <- which(epi_sim$time <= time)
      index_time <- lessthan_time[length(lessthan_time)]
      Iseries[i] <- rbinom(1, epi_sim[index_time,]$I, 0.05) # apply case detection rate to determine number of cases observed - only detect 5% of cases
      i <- i + 1
    }
    
    return(Iseries)
  }
  else  # no major epidemic, return all NAs
    return(rep(NA, max_length))
}


# function to get E over epidemic adjusted to longest time series
get_Etimeseries <- function(epi_sim, cv_s, max_length)
{
  Eseries <- rep(0, max_length)
  
  # set threshold for major outbreak based on if epidemic lasts long enough to match the data
  # if correlation isn't positive, won't get a long epidemic so use all data regardless
  if(cor_level == 1 | R0 == 5.2)
    threshold <- 365
  else
    threshold <- 0
  
  # check if epidemic lasts longer than threshold time
  if(epi_sim$time[nrow(epi_sim)] >= threshold)  # major epidemic, set number inf for each time step
  {
    i <- 1
    for(time in seq(0, max(epi_sim$time)+0.1, 0.1))
    {
      lessthan_time <- which(epi_sim$time <= time)
      index_time <- lessthan_time[length(lessthan_time)]
      Eseries[i] <- epi_sim[index_time,]$E
      i <- i + 1
    }  
    
    return(Eseries)
  }
  else  # no major epidemic, return all NAs
    return(rep(NA, max_length))
}



# function to get S over epidemic adjusted to longest time series
get_Stimeseries <- function(epi_sim, cv_s, max_length)
{
  Sseries <- rep(epi_sim$S[nrow(epi_sim)], max_length)
  
  # set threshold for major outbreak based on if epidemic lasts long enough to match the data
  # if correlation isn't positive, won't get a long epidemic so use all data regardless
  if(cor_level == 1 | R0 == 5.2)
    threshold <- 365
  else
    threshold <- 0
  
  # check if epidemic lasts longer than threshold time
  if(epi_sim$time[nrow(epi_sim)] >= threshold)  # major epidemic, set number inf for each time step
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
cv_ts <- c(3)
cv_ss <- c(3)
cors <- c(-1, 0, 1)
R0 <- 0.52

# make empty lists to store all data  
all_Iseries <- vector("list", length=length(cv_ts)*length(cv_ss)*length(cors))
all_Iobsseries <- vector("list", length=length(cv_ts)*length(cv_ss)*length(cors))
all_Eseries <- vector("list", length=length(cv_ts)*length(cv_ss)*length(cors))
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
      filePattern <- paste("results_his",cv_s,"hit",cv_t,"cor",corr,"g0.071e0.13b70d0.000997p70180R00.52initE19s*", sep="")
      files <- list.files(pattern=filePattern)
      
      # make a list that contains all files
      read.all <- lapply(files, read_sav)
      
      # determine length of longest time series based on time
      max_time <- max(unlist(lapply(read.all, function(x) {x$time[nrow(x)]})))
      max_length <- length(seq(0,max_time+0.1,0.1))
      
      ## time series for I indivs
      # make data frame containing series of I indivs for each time series of the sims
      temp <- data.frame(sapply(read.all, get_Itimeseries, cor_level=corr, max_length=max_length))

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
      
      
      ## time series for I obs indivs
      # make data frame containing series of I indivs for each time series of the sims
      temp <- data.frame(sapply(read.all, get_Iobstimeseries, cor_level=corr, R0=R0, max_length=max_length))
      
      # make sure all numeric
      temp <- data.frame(sapply(temp, as.numeric))
      
      # find mean, 95% CIs for each time point
      I_mean <- apply(temp, 1, mean, na.rm=T)
      I_low <- apply(temp, 1, function(x){quantile(x, probs=0.025, na.rm=T)})
      I_high <- apply(temp, 1, function(x){quantile(x, probs=0.975, na.rm=T)})
      
      # make data frame containing time, mean, and 95% CIs
      element_name <- paste(cv_s,cv_t,corr,R0, sep="_")
      Iobsseries <- data.frame(time=seq(0,max_time+0.1,0.1),
                               avg=I_mean,
                               low=I_low,
                               high=I_high)
      
      # add data frame to list of data frames for each param combo (cv_s, cv_t, cor) and name
      all_Iobsseries[[combo]] <- Iobsseries
      names(all_Iobsseries)[[combo]] <- element_name
      
      
      ## time series for E indivs
      # make data frame containing series of I indivs for each time series of the sims
      temp <- data.frame(sapply(read.all, get_Etimeseries, cv_s=cv_s, max_length=max_length))

      # make sure all numeric
      temp <- data.frame(sapply(temp, as.numeric))

      # find mean, 95% CIs for each time point
      E_mean <- apply(temp, 1, mean, na.rm=T)
      E_low <- apply(temp, 1, function(x){quantile(x, probs=0.025, na.rm=T)})
      E_high <- apply(temp, 1, function(x){quantile(x, probs=0.975, na.rm=T)})

      # make data frame containing time, mean, and 95% CIs
      element_name <- paste(cv_s,cv_t,corr, sep="_")
      Eseries <- data.frame(time=seq(0,max_time+0.1,0.1),
                           avg=E_mean,
                           low=E_low,
                           high=E_high)

      # add data frame to list of data frames for each param combo (cv_s, cv_t, cor) and name
      all_Eseries[[combo]] <- Eseries
      names(all_Eseries)[[combo]] <- element_name


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



# levels of heterogeneity in transmission, susceptibility, and correlation
cv_ts <- c(3)
cv_ss <- c(3)
cors <- c(0)
R0 <- 5.2


## get data for each level of HiS, HiT, and correlation for a specific R0
for(cv_s in cv_ss)
{
  for(cv_t in cv_ts)
  {
    for(corr in cors)
    {
      # list of file names
      filePattern <- paste("results_his",cv_s,"hit",cv_t,"cor",corr,"g0.071e0.13b70d0.000997p70180R05.2initE19s*", sep="")
      files <- list.files(pattern=filePattern)
      
      # make a list that contains all files
      read.all <- lapply(files, read_sav)
      
      # determine length of longest time series based on time
      max_time <- max(unlist(lapply(read.all, function(x) {x$time[nrow(x)]})))
      max_length <- length(seq(0,max_time+0.1,0.1))
      
      ## time series for I obs indivs
      # make data frame containing series of I indivs for each time series of the sims
      temp <- data.frame(sapply(read.all, get_Iobstimeseries, cor_level=corr, R0=R0, max_length=max_length))
      
      # make sure all numeric
      temp <- data.frame(sapply(temp, as.numeric))
      
      # find mean, 95% CIs for each time point
      I_mean <- apply(temp, 1, mean, na.rm=T)
      I_low <- apply(temp, 1, function(x){quantile(x, probs=0.025, na.rm=T)})
      I_high <- apply(temp, 1, function(x){quantile(x, probs=0.975, na.rm=T)})
      
      # make data frame containing time, mean, and 95% CIs
      element_name <- paste(cv_s,cv_t,corr,R0, sep="_")
      Iobsseries <- data.frame(time=seq(0,max_time+0.1,0.1),
                               avg=I_mean,
                               low=I_low,
                               high=I_high)
      
      # add data frame to list of data frames for each param combo (cv_s, cv_t, cor) and name
      all_Iobsseries[[combo]] <- Iobsseries
      names(all_Iobsseries)[[combo]] <- element_name
      
      combo <- combo + 1
    }
  }
}




# levels of heterogeneity in transmission, susceptibility, and correlation
cv_ts <- c(3)
cv_ss <- c(3)
cors <- c(-1)
R0 <- 5.2


## get data for each level of HiS, HiT, and correlation for a specific R0
for(cv_s in cv_ss)
{
  for(cv_t in cv_ts)
  {
    for(corr in cors)
    {
      # list of file names
      filePattern <- paste("Iobs_results_his",cv_s,"hit",cv_t,"cor",corr,"g0.071e0.13b70d0.000997p70180R05.2initE19s*", sep="")
      files <- list.files(pattern=filePattern)
      
      # make a list that contains all files
      read.all <- lapply(files, read_sav)
      
      ## time series for I obs indivs
      # make data frame containing series of I indivs for each time series of the sims
      temp <- data.frame(read.all)
      
      # make sure all numeric
      temp <- data.frame(sapply(temp, as.numeric))
      
      # find mean, 95% CIs for each time point
      I_mean <- apply(temp, 1, mean, na.rm=T)
      I_low <- apply(temp, 1, function(x){quantile(x, probs=0.025, na.rm=T)})
      I_high <- apply(temp, 1, function(x){quantile(x, probs=0.975, na.rm=T)})
      
      # make data frame containing time, mean, and 95% CIs
      element_name <- paste(cv_s,cv_t,corr,R0, sep="_")
      Iobsseries <- data.frame(time=seq(0,max_time+0.1,0.1),
                               avg=I_mean,
                               low=I_low,
                               high=I_high)
      
      # add data frame to list of data frames for each param combo (cv_s, cv_t, cor) and name
      all_Iobsseries[[combo]] <- Iobsseries
      names(all_Iobsseries)[[combo]] <- element_name
      
      combo <- combo + 1
    }
  }
}







############ plot case counts and model simulations for NYC #######################
# set working directory to get data from
setwd("C:/Users/bmt5507/Documents")

# read in data
case_counts_nyc22 <- read.csv("mpox_nyc_22.csv", header=T)
case_counts_nyc23plus <- read.csv("mpox_nyc_23plus.csv", header=T)

# convert date into readable date
case_counts_nyc22 <- case_counts_nyc22 %>% mutate(diagnosis_date = mdy(diagnosis_date)) 
case_counts_nyc22 <- case_counts_nyc22 %>% mutate(int_time = row_number())

colnames(case_counts_nyc23plus) <- c("diagnosis_date", "count")
case_counts_nyc23plus <- case_counts_nyc23plus %>% mutate(diagnosis_date = mdy(diagnosis_date)) 
case_counts_nyc23plus <- case_counts_nyc23plus %>% mutate(int_time = row_number()+nrow(case_counts_nyc22))

# combine all case count data across years
case_counts_nyc_all <- rbind(case_counts_nyc22[,c("diagnosis_date","count","int_time")], case_counts_nyc23plus)



### plot model simulations and NYC case count data for R0 defined by Eq 6
par(mfrow=c(1,1), mar=c(9.5,5.1,1.1,0.5), mgp=c(3,1,0))
plot(case_counts_nyc_all$int_time+16, case_counts_nyc_all$count, xlim=c(0,1042), ylim=c(0,140), xaxt="n", xlab="", ylab="Case count", cex.axis=1.5, cex.lab=2)
lines(all_Iobsseries$'3_3_1_0.52'$time[which(all_Iobsseries$'3_3_1_0.52'$time <= 1042)], all_Iobsseries$'3_3_1_0.52'$avg[which(all_Iobsseries$'3_3_1_0.52'$time <= 1042)], col="blue", lwd=1)
polygon(c(rev(all_Iobsseries$'3_3_1_0.52'$time[which(all_Iobsseries$'3_3_1_0.52'$time <= 1042)]), all_Iobsseries$'3_3_1_0.52'$time[which(all_Iobsseries$'3_3_1_0.52'$time <= 1042)]), c(rev(all_Iobsseries$'3_3_1_0.52'$low[which(all_Iobsseries$'3_3_1_0.52'$time <= 1042)]), all_Iobsseries$'3_3_1_0.52'$high[which(all_Iobsseries$'3_3_1_0.52'$time <= 1042)]),
        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
lines(all_Iobsseries$'3_3_0_0.52'$time[which(all_Iobsseries$'3_3_0_0.52'$time <= 1042)], all_Iobsseries$'3_3_0_0.52'$avg[which(all_Iobsseries$'3_3_0_0.52'$time <= 1042)], col="yellow3", lwd=1, lty="solid")
polygon(c(rev(all_Iobsseries$'3_3_0_0.52'$time[which(all_Iobsseries$'3_3_0_0.52'$time <= 1042)]), all_Iobsseries$'3_3_0_0.52'$time[which(all_Iobsseries$'3_3_0_0.52'$time <= 1042)]), c(rev(all_Iobsseries$'3_3_0_0.52'$low[which(all_Iobsseries$'3_3_0_0.52'$time <= 1042)]), all_Iobsseries$'3_3_0_0.52'$high[which(all_Iobsseries$'3_3_0_0.52'$time <= 1042)]),
        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Iobsseries$'3_3_-1_0.52'$time[which(all_Iobsseries$'3_3_-1_0.52'$time <= 1042)], all_Iobsseries$'3_3_-1_0.52'$avg[which(all_Iobsseries$'3_3_-1_0.52'$time <= 1042)], col="red", lwd=1, lty="solid")
polygon(c(rev(all_Iobsseries$'3_3_-1_0.52'$time[which(all_Iobsseries$'3_3_-1_0.52'$time <= 1042)]), all_Iobsseries$'3_3_-1_0.52'$time[which(all_Iobsseries$'3_3_-1_0.52'$time <= 1042)]), c(rev(all_Iobsseries$'3_3_-1_0.52'$low[which(all_Iobsseries$'3_3_-1_0.52'$time <= 1042)]), all_Iobsseries$'3_3_-1_0.52'$high[which(all_Iobsseries$'3_3_-1_0.52'$time <= 1042)]),
        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
axis(1, at=c(0,cumsum(as.numeric(diff(month_dates)))), labels=format(month_dates, "%b %Y"), las = 2, cex.axis=1.5)
legend("topright", bty="n", legend=c("Model, pos cor", "Model, no cor", "Model, neg cor", "Data"),
       col=c("blue","yellow3","red","black"),
       lty=c("solid","solid","solid",NA), pch=c(NA,NA,NA,1), lwd=2, cex=1.7)
title(xlab="Date", line=8.3, cex.lab=2.2)

# save pdf: 7.39 x 7


### plot model simulations and NYC case count data for R0 defined by Eq 6 or 7
# set up x axis
start_date <- as.Date("2022-05-02")
month_dates <- seq(start_date, by="3 months", length.out=12)
# set up plot configuration
par(mfrow=c(1,3), mar=c(9.5,5.1,1.1,0.5), mgp=c(3,1,0))

# plot model simulated with positive correlation
plot(case_counts_nyc_all$int_time+16, case_counts_nyc_all$count, xlim=c(0,1042), ylim=c(0,140), xaxt="n", xlab="", ylab="Case count", cex.axis=1.5, cex.lab=2)
lines(all_Iobsseries$'3_3_1_0.52'$time[which(all_Iobsseries$'3_3_1_0.52'$time <= 1042)], all_Iobsseries$'3_3_1_0.52'$avg[which(all_Iobsseries$'3_3_1_0.52'$time <= 1042)], col="blue", lwd=1)
polygon(c(rev(all_Iobsseries$'3_3_1_0.52'$time[which(all_Iobsseries$'3_3_1_0.52'$time <= 1042)]), all_Iobsseries$'3_3_1_0.52'$time[which(all_Iobsseries$'3_3_1_0.52'$time <= 1042)]), c(rev(all_Iobsseries$'3_3_1_0.52'$low[which(all_Iobsseries$'3_3_1_0.52'$time <= 1042)]), all_Iobsseries$'3_3_1_0.52'$high[which(all_Iobsseries$'3_3_1_0.52'$time <= 1042)]),
        col=adjustcolor("blue", alpha.f=0.2) , lty = 0)
axis(1, at=c(0,cumsum(as.numeric(diff(month_dates)))), labels=format(month_dates, "%b %Y"), las = 2, cex.axis=1.5)
legend("topright", bty="n", legend=c("Model, pos cor", "Data"),
       col=c("blue","black"),
       lty=c("solid",NA), pch=c(NA,1), lwd=2, cex=1.7)
title(xlab="Date", line=8.3, cex.lab=2.2)
text(x=0, y=138, label="a", cex=2)

# plot model simulated with no correlation
plot(case_counts_nyc_all$int_time+16, case_counts_nyc_all$count, xlim=c(0,1042), ylim=c(0,140), xaxt="n", xlab="", ylab="", cex.axis=1.5, cex.lab=2)
lines(all_Iobsseries$'3_3_0_0.52'$time[which(all_Iobsseries$'3_3_0_0.52'$time <= 1042)], all_Iobsseries$'3_3_0_0.52'$avg[which(all_Iobsseries$'3_3_0_0.52'$time <= 1042)], col="yellow3", lwd=1, lty="solid")
polygon(c(rev(all_Iobsseries$'3_3_0_0.52'$time[which(all_Iobsseries$'3_3_0_0.52'$time <= 1042)]), all_Iobsseries$'3_3_0_0.52'$time[which(all_Iobsseries$'3_3_0_0.52'$time <= 1042)]), c(rev(all_Iobsseries$'3_3_0_0.52'$low[which(all_Iobsseries$'3_3_0_0.52'$time <= 1042)]), all_Iobsseries$'3_3_0_0.52'$high[which(all_Iobsseries$'3_3_0_0.52'$time <= 1042)]),
        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
lines(all_Iobsseries$'3_3_0_5.2'$time[which(all_Iobsseries$'3_3_0_5.2'$time <= 1042)], all_Iobsseries$'3_3_0_5.2'$avg[which(all_Iobsseries$'3_3_0_5.2'$time <= 1042)], col="yellow3", lwd=1)
polygon(c(rev(all_Iobsseries$'3_3_0_5.2'$time[which(all_Iobsseries$'3_3_0_5.2'$time <= 1042)]), all_Iobsseries$'3_3_0_5.2'$time[which(all_Iobsseries$'3_3_0_5.2'$time <= 1042)]), c(rev(all_Iobsseries$'3_3_0_5.2'$low[which(all_Iobsseries$'3_3_0_5.2'$time <= 1042)]), all_Iobsseries$'3_3_0_5.2'$high[which(all_Iobsseries$'3_3_0_5.2'$time <= 1042)]),
        col=adjustcolor("yellow3", alpha.f=0.2) , lty = 0)
axis(1, at=c(0,cumsum(as.numeric(diff(month_dates)))), labels=format(month_dates, "%b %Y"), las = 2, cex.axis=1.5)
legend("topright", bty="n", legend=c("Model, no cor, 1", "Model, no cor, 2", "Data"),
       col=c("yellow3","yellow3","black"),
       lty=c("dashed","solid",NA), pch=c(NA,NA,1), lwd=2, cex=1.7)
title(xlab="Date", line=8.3, cex.lab=2.2)
text(x=0, y=138, label="b", cex=2)

# plot model simulated with negative correlation
plot(case_counts_nyc_all$int_time+16, case_counts_nyc_all$count, xlim=c(0,1042), ylim=c(0,140), xaxt="n", xlab="", ylab="", cex.axis=1.5, cex.lab=2)
lines(all_Iobsseries$'3_3_-1_0.52'$time[which(all_Iobsseries$'3_3_-1_0.52'$time <= 1042)], all_Iobsseries$'3_3_-1_0.52'$avg[which(all_Iobsseries$'3_3_-1_0.52'$time <= 1042)], col="red", lwd=1, lty="solid")
polygon(c(rev(all_Iobsseries$'3_3_-1_0.52'$time[which(all_Iobsseries$'3_3_-1_0.52'$time <= 1042)]), all_Iobsseries$'3_3_-1_0.52'$time[which(all_Iobsseries$'3_3_-1_0.52'$time <= 1042)]), c(rev(all_Iobsseries$'3_3_-1_0.52'$low[which(all_Iobsseries$'3_3_-1_0.52'$time <= 1042)]), all_Iobsseries$'3_3_-1_0.52'$high[which(all_Iobsseries$'3_3_-1_0.52'$time <= 1042)]),
        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
lines(all_Iobsseries$'3_3_-1_5.2'$time[which(all_Iobsseries$'3_3_-1_5.2'$time <= 1042)], all_Iobsseries$'3_3_-1_5.2'$avg[which(all_Iobsseries$'3_3_-1_5.2'$time <= 1042)], col="red", lwd=1)
polygon(c(rev(all_Iobsseries$'3_3_-1_5.2'$time[which(all_Iobsseries$'3_3_-1_5.2'$time <= 1042)]), all_Iobsseries$'3_3_-1_5.2'$time[which(all_Iobsseries$'3_3_-1_5.2'$time <= 1042)]), c(rev(all_Iobsseries$'3_3_-1_5.2'$low[which(all_Iobsseries$'3_3_-1_5.2'$time <= 1042)]), all_Iobsseries$'3_3_-1_5.2'$high[which(all_Iobsseries$'3_3_-1_5.2'$time <= 1042)]),
        col=adjustcolor("red", alpha.f=0.2) , lty = 0)
axis(1, at=c(0,cumsum(as.numeric(diff(month_dates)))), labels=format(month_dates, "%b %Y"), las = 2, cex.axis=1.5)
legend("topright", bty="n", legend=c("Model, neg cor, 1", "Model, neg cor, 2", "Data"),
       col=c("red","red","black"),
       lty=c("dashed","solid",NA), pch=c(NA,NA,1), lwd=2, cex=1.7)
title(xlab="Date", line=8.3, cex.lab=2.2)
text(x=0, y=138, label="c", cex=2)

# save pdf: 15 x 5.75




