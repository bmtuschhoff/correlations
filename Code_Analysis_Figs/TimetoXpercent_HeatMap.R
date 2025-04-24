# read in all files for individual-based SIR model for each parameter combination for the level of heterogeneity in susceptibility (cv_s), heterogeneity in transmission (cv_t), and correlation (cor)
# determine mean, median, 95% CI for time to X percent infected of the final epidemic size given a major outbreak
# plot results

library(dplyr)
library(ggplot2)
library(metR)
library(ggpubr)
library(egg)
library(cowplot)

# number of simulations
sims <- 500

# population size for each simulation
pop_size <- 1000

# set working directory to get data from
setwd("C:/Users/bmt5507/Documents/Cor_HiST_Copula_I1_R0.8")

# function to read in datasets that are tab delimited
read_sav <- function(fileName) {read.table(fileName, header=T, sep="\t")}


# function to get time to X% of final epidemic size
get_timetoX <- function(epi_sim, pop_size, cv_s, cv_t, corr, sim)
{
  # find final epidemic size (FES)
  FES <- epi_sim$I[nrow(epi_sim)] + epi_sim$R[nrow(epi_sim)]
  
  # check if cumulative number of inf indivs is at least 20% of final epidemic size
  for(row in 1:nrow(epi_sim))
  {
    if((pop_size - epi_sim$S[row]) >= 0.2*FES)
    {
      timeto20 <- epi_sim$time[row]
      break
    }
  }
  
  # check if cumulative number of inf indivs is at least 50% of final epidemic size
  for(row in 1:nrow(epi_sim))
  {
    if((pop_size - epi_sim$S[row]) >= 0.5*FES)
    {
      timeto50 <- epi_sim$time[row]
      break
    }
  }
  
  # check if cumulative number of inf indivs is at least 75% of final epidemic size
  for(row in 1:nrow(epi_sim))
  {
    if((pop_size - epi_sim$S[row]) >= 0.75*FES)
    {
      timeto75 <- epi_sim$time[row]
      break
    }
  }
  
  # determine whether there was a major outbreak
  # set threshold for major outbreak based on level of HiS because more HiS = smaller outbreak but not less likely necessarily
  if(cv_s == 0 || cv_s == 0.5)
    threshold <- 200
  if(cv_s == 1)
    threshold <- 100
  if(cv_s == 3)
    threshold <- 50
  
  # check if final epidemic size (FES) is larger than threshold for major epidemic
  if(FES > threshold)  # major epidemic
    major_out <- 1
  else  # no major epidemic
    major_out <- 0
  
  return(c(cv_s, cv_t, corr, sim, timeto20, timeto50, timeto75, major_out))
}


# levels of heterogeneity in transmission, susceptibility, and correlation
cv_ts <- c(0.5, 1, 3)
cv_ss <- c(0.5, 1, 3)
cors <- c(-1,-0.5,0,0.5,1)


# make empty dataframe to store all data  
timetoX_R3 <- data.frame(cv_s = double(), 
                      cv_t = double(), 
                      cor = character(), 
                      sim = integer(),
                      timeto20 = double(),
                      timeto50 = double(),
                      timeto75 = double(),
                      major_out = integer(),
                      stringsAsFactors = FALSE) 

## get data for each level of HiS, HiT, and correlation for a specific R0
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
      
      # make data frame containing peak size, peak time, FES of each sim for one param combo
      temp <- data.frame(t(as.matrix(sapply(read.all, get_timetoX, pop_size=pop_size, cv_s=cv_s, cv_t=cv_t, corr=corr, sim=0))))
      
      # make all numeric columns numeric and add column names
      temp[,-3] <- data.frame(sapply(temp[,-3], as.numeric))
      colnames(temp) <- c("cv_s", "cv_t", "cor", "sim", "timeto20", "timeto50", "timeto75", "major_out")
      
      # add correct simulation numbers
      temp$sim <- seq(1,500,1)
      
      # remove all sims with NA values
      temp <- na.omit(temp)
      
      # add temp to overall dataframe
      timetoX_R3 <- rbind(timetoX_R3, temp)
    }
  }
}


########################################################################################
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

# levels of heterogeneity in transmission and susceptibility
cv_ss <- c(0,0.5,1,3,0,0,0)
cv_ts <- c(0,0,0,0,0.5,1,3)

# make empty dataframe to store this data separately
timetoX_R3_hom <- data.frame(cv_s = double(), 
                         cv_t = double(), 
                         cor = character(), 
                         sim = integer(),
                         timeto20 = double(),
                         timeto50 = double(),
                         timeto75 = double(),
                         major_out = integer(),
                         stringsAsFactors = FALSE) 

# get data for each level of HiS, HiT, and correlation for a specific R0
for(i in 1:7)
{
  files <- list.files(pattern=filePatterns[i])
  
  # make a list that contains all files
  read.all <- lapply(files, read_sav)
  
  # make data frame containing peak size, peak time, FES of each sim for one param combo
  temp <- data.frame(t(as.matrix(sapply(read.all, get_timetoX, pop_size=pop_size, cv_s=cv_ss[i], cv_t=cv_ts[i], corr="0", sim=0))))
  
  # make all numeric columns numeric and add column names
  temp[,-3] <- data.frame(sapply(temp[,-3], as.numeric))
  colnames(temp) <- c("cv_s", "cv_t", "cor", "sim", "timeto20", "timeto50", "timeto75", "major_out")
  
  # add correct simulation numbers
  temp$sim <- seq(1,500,1)
  
  # remove all sims with NA values
  temp <- na.omit(temp)
  
  # add temp to overall dataframe
  timetoX_R3_hom <- rbind(timetoX_R3_hom, temp)
}



# summarize data
timetoX_R3_mean <- timetoX_R3 %>% filter(major_out != 0) %>% 
  group_by(across(all_of(c("cv_s", "cv_t", "cor")))) %>%
  summarise_at(c("timeto20", "timeto50", "timeto75"), mean, na.rm = TRUE)

timetoX_R3_hom_mean <- timetoX_R3_hom %>% filter(major_out != 0) %>% 
  group_by(across(all_of(c("cv_s", "cv_t", "cor")))) %>%
  summarise_at(c("timeto20", "timeto50", "timeto75"), mean, na.rm = TRUE)

# add back in data points that don't have major epidemics as "NA"
for(cv_s in c(0.5,1,3))
{
  for(cv_t in c(0.5,1,3))
  {
    for(corr in cors)
    {
      # check if this parameter set is in the summary dataframe, if not (returns TRUE) add with NA for each epi measure
      if(identical(which(timetoX_R3_mean$cv_s == cv_s & timetoX_R3_mean$cv_t == cv_t & timetoX_R3_mean$cor == corr), integer(0)))
      {
        new_row <- data.frame(cv_s=cv_s, cv_t=cv_t, cor=corr, timeto20=NA, timeto50=NA, timeto75=NA)
        timetoX_R3_mean <- rbind(timetoX_R3_mean, new_row)
      } 
    }
  }
}

for(i in 1:7)
{
  # check if this parameter set is in the summary dataframe, if not (returns TRUE) add with NA for each epi measure
  if(identical(which(timetoX_R3_hom_mean$cv_s == cv_ss[i] & timetoX_R3_hom_mean$cv_t == cv_ts[i]), integer(0)))
  {
    new_row <- data.frame(cv_s=cv_ss[i], cv_t=cv_ts[i], cor="0", timeto20=NA, timeto50=NA, timeto75=NA)
    timetoX_R3_hom_mean <- rbind(timetoX_R3_hom_mean, new_row)
  } 
}


# create factors for plotting
timetoX_R3_mean$cv_s <- as.factor(timetoX_R3_mean$cv_s)
timetoX_R3_mean$cv_t <- as.factor(timetoX_R3_mean$cv_t)
timetoX_R3_mean$cor <- factor(timetoX_R3_mean$cor, levels = c(-1,-0.5,0,0.5,1))

timetoX_R3_hom_mean$cv_s <- as.factor(timetoX_R3_hom_mean$cv_s)
timetoX_R3_hom_mean$cv_t <- as.factor(timetoX_R3_hom_mean$cv_t)
timetoX_R3_hom_mean$cor <- factor(timetoX_R3_hom_mean$cor, levels = c(-1,-0.5,0,0.5,1))


######################################################################################
################################# plot results #######################################
######################################################################################

#### plot time to X% of FES
timetoX_corn1 <- ggplot(timetoX_R3_mean[which(timetoX_R3_mean$cor == -1),], aes(x=cv_s, y=cv_t, fill=timeto50)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(13,101),  breaks=seq(13,101,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to 50%") +
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
timetoX_corn1 <- set_panel_size( timetoX_corn1, height=unit(3*2, "cm"),
                                 width=unit(3*2, "cm") )

timetoX_corn0.5 <- ggplot(timetoX_R3_mean[which(timetoX_R3_mean$cor == -0.5),], aes(x=cv_s, y=cv_t, fill=timeto50)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(13,101),  breaks=seq(13,101,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to 50%") +
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
timetoX_corn0.5 <- set_panel_size( timetoX_corn0.5, height=unit(3*2, "cm"),
                                   width=unit(3*2, "cm") )

timetoX_cor0 <- ggplot(timetoX_R3_mean[which(timetoX_R3_mean$cor == 0 & timetoX_R3_mean$cv_s != 0 & timetoX_R3_mean$cv_t != 0),], aes(x=cv_s, y=cv_t, fill=timeto50)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(13,101),  breaks=seq(13,101,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to 50%") +
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
timetoX_cor0 <- set_panel_size( timetoX_cor0, height=unit(3*2, "cm"),
                                width=unit(3*2, "cm") )

timetoX_cor0.5 <- ggplot(timetoX_R3_mean[which(timetoX_R3_mean$cor == 0.5),], aes(x=cv_s, y=cv_t, fill=timeto50)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(13,101),  breaks=seq(13,101,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to 50%") +
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
timetoX_cor0.5 <- set_panel_size( timetoX_cor0.5, height=unit(3*2, "cm"),
                                  width=unit(3*2, "cm") )

timetoX_cor1 <- ggplot(timetoX_R3_mean[which(timetoX_R3_mean$cor == 1),], aes(x=cv_s, y=cv_t, fill=timeto50)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(13,101),  breaks=seq(13,101,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to 50%") +
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
timetoX_cor1 <- set_panel_size( timetoX_cor1, height=unit(3*2, "cm"),
                                width=unit(3*2, "cm") )

timetoX_cs0 <- ggplot(timetoX_R3_hom_mean[which(timetoX_R3_hom_mean$cv_s == 0 & timetoX_R3_hom_mean$cv_t != 0),], aes(x=cor, y=cv_t, fill=timeto50)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(13,101),  breaks=seq(13,101,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to 50%") +
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
timetoX_cs0 <- set_panel_size( timetoX_cs0, height=unit(3*2, "cm"),
                               width=unit(1*2, "cm") )

timetoX_cs0ct0 <- ggplot(timetoX_R3_hom_mean[which(timetoX_R3_hom_mean$cv_s == 0 & timetoX_R3_hom_mean$cv_t == 0),], aes(x=cor, y=cv_t, fill=timeto50)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(13,101),  breaks=seq(13,101,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to 50%") +
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
timetoX_cs0ct0 <- set_panel_size( timetoX_cs0ct0, height=unit(1*2, "cm"),
                                  width=unit(1*2, "cm") )


timetoX_ct0 <- ggplot(timetoX_R3_hom_mean[which(timetoX_R3_hom_mean$cv_s != 0 & timetoX_R3_hom_mean$cv_t == 0),], aes(x=cv_s, y=cv_t, fill=timeto50)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(13,101),  breaks=seq(13,101,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to 50%") +
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
timetoX_ct0 <- set_panel_size( timetoX_ct0, height=unit(1*2, "cm"),
                               width=unit(3*2, "cm") )


# determine legend for combined plots
timetoX_leg <- ggplot(timetoX_R3_mean[which(timetoX_R3_mean$cor == 0),], aes(x=cv_t, y=cv_s, fill=timeto50)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  #scale_fill_steps(low="white", high="red", limits=c(6,70),  breaks=seq(6,70,4),  # R0=3
  scale_fill_steps(low="white", high="red", limits=c(13,101),  breaks=seq(13,101,8),  # R0=0.8
  #scale_fill_steps(low="white", high="red", limits=c(11,137),  breaks=seq(11,137,14),  # R0=1.1
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Time to 50%") +
  xlab(expression(paste("Coefficient of variation of transmission (",italic(C)[italic(t)],")"))) +
  ylab(expression(paste("Coefficient of variation of susceptibility ( ",italic(C)[italic(s)],")"))) +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), legend.key.height=unit(2.5,"cm"))
legend <- get_legend(timetoX_leg)

# plot all combined plots in a grid with common legend
plt <- plot_grid(timetoX_cs0, timetoX_corn1, timetoX_corn0.5, timetoX_cor0, timetoX_cor0.5, timetoX_cor1,
                 NULL, NULL, NULL, NULL, NULL, NULL,
                 timetoX_cs0ct0, NULL, NULL, timetoX_ct0, NULL, NULL,
                 nrow=3, ncol=6, rel_heights=c(1,-0.5,1), rel_widths=c(0.4,1,1,1,1,1))

plot_grid(plt, legend,
          nrow=1, rel_widths=c(1,0.2))
# size: 19 x 7.5

