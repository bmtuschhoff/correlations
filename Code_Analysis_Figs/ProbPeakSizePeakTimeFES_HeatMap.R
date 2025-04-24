# read in all files for individual-based SIR model for each parameter combination for the level of heterogeneity in susceptibility (cv_s), heterogeneity in transmission (cv_t), and correlation (cor)
# determine mean, median, 95% CI for probability of a major outbreak and peak size, peak time, and final epidemic size given a major outbreak
# plot results

library(dplyr)
library(ggplot2)
library(ggpubr)
library(egg)
library(cowplot)

# number of simulations
sims <- 500

# set working directory to get data from
setwd("C:/Users/bmt5507/Documents/Cor_HiST_Copula_I1")

# function to read in datasets that are tab delimited
read_sav <- function(fileName) {read.table(fileName, header=T, sep="\t")}


# function to get peak size, peak time, and final epidemic size (FES) from each sim for a param combo
get_data <- function(epi_sim, cv_s, cv_t, cor, sim)
{
  peak_size <- max(epi_sim$I)
  peak_time <- epi_sim$time[which.max(epi_sim$I)]
  FES <- epi_sim$I[nrow(epi_sim)] + epi_sim$R[nrow(epi_sim)]
  
  # determine whether there was a major outbreak
  # set threshold for major outbreak based on level of HiS because more HiS = smaller outbreak
  if(cv_s == 0 || cv_s == 0.5)
    threshold <- 200
  if(cv_s == 1)
    threshold <- 100
  if(cv_s == 3)
    threshold <- 50
  
  if(FES > threshold)
    major_out <- 1
  else
    major_out <- 0
  
  return(c(cv_s, cv_t, cor, sim, peak_size, peak_time, FES, major_out))
}


# levels of heterogeneity in transmission, susceptibility, and correlation
cv_ts <- c(0.5, 1, 3)
cv_ss <- c(0.5, 1, 3)
cors <- c(-1,-0.5,0,0.5,1)

# make empty dataframe to store all data  
data_R3 <- data.frame(cv_s = double(), 
                 cv_t = double(), 
                 cor = character(), 
                 sim = integer(),
                 peak_size = integer(),
                 peak_time = double(),
                 FES = integer(),
                 major_out = double(),
                 stringsAsFactors = FALSE) 

## get data for each level of HiS, HiT, and correlation for a specific R0
for(cv_s in cv_ss)
{
  for(cv_t in cv_ts)
  {
    for(corr in cors)
    {
      # list of file names
      filePattern <- paste("results_his",cv_s,"hit",cv_t,"cor",corr,"g0.1p1000R03initI1s*", sep="")
      files <- list.files(pattern=filePattern)
      
      # make a list that contains all files
      read.all <- lapply(files, read_sav)
      
      # make data frame containing peak size, peak time, FES of each sim for one param combo
      temp <- data.frame(t(as.matrix(sapply(read.all, get_data, cv_s=cv_s, cv_t=cv_t, cor=corr, sim=0))))
      
      # make all numeric columns numeric and add column names
      temp[,-3] <- data.frame(sapply(temp[,-3], as.numeric))
      colnames(temp) <- c("cv_s", "cv_t", "cor", "sim", "peak_size", "peak_time", "FES", "major_out")
      
      # add correct simulation numbers
      temp$sim <- seq(1,500,1)
      
      # remove all sims with NA values
      temp <- na.omit(temp)
      
      # add temp to overall dataframe
      data_R3 <- rbind(data_R3, temp)
    }
  }
}

## separately get data for homogeneous case (HiS=0, HiT=0) and cases with just HiS or just HiT
# list of file names
filePatterns <- c()
filePatterns[1] <- paste("results_his0hit0cor0g0.1p1000R03initI1s*", sep="")
filePatterns[2] <- paste("results_his0.5hit0cor0g0.1p1000R03initI1s*", sep="")
filePatterns[3] <- paste("results_his1hit0cor0g0.1p1000R03initI1s*", sep="")
filePatterns[4] <- paste("results_his3hit0cor0g0.1p1000R03initI1s*", sep="")
filePatterns[5] <- paste("results_his0hit0.5cor0g0.1p1000R03initI1s*", sep="")
filePatterns[6] <- paste("results_his0hit1cor0g0.1p1000R03initI1s*", sep="")
filePatterns[7] <- paste("results_his0hit3cor0g0.1p1000R03initI1s*", sep="")

# levels of heterogeneity
cv_ss <- c(0,0.5,1,3,0,0,0)
cv_ts <- c(0,0,0,0,0.5,1,3)

# make empty dataframe to store this data separately
data_R3_hom <- data.frame(cv_s = double(), 
                      cv_t = double(), 
                      cor = character(), 
                      sim = integer(),
                      peak_size = integer(),
                      peak_time = double(),
                      FES = integer(),
                      major_out = double(),
                      stringsAsFactors = FALSE) 

# get data for each level of HiS, HiT, and correlation
for(i in 1:7)
{
  files <- list.files(pattern=filePatterns[i])
  
  # make a list that contains all files
  read.all <- lapply(files, read_sav)
  
  # make data frame containing peak size, peak time, FES of each sim for one param combo
  temp <- data.frame(t(as.matrix(sapply(read.all, get_data, cv_s=cv_ss[i], cv_t=cv_ts[i], cor="0", sim=0))))
  
  # make all numeric columns numeric and add column names
  temp[,-3] <- data.frame(sapply(temp[,-3], as.numeric))
  colnames(temp) <- c("cv_s", "cv_t", "cor", "sim", "peak_size", "peak_time", "FES", "major_out")
  
  # add correct simulation numbers
  temp$sim <- seq(1,500,1)
  
  # remove all sims with NA values
  temp <- na.omit(temp)
  
  # add temp to overall dataframe
  data_R3_hom <- rbind(data_R3_hom, temp)
}


# check FES for determining major outbreaks
hist(data_R3[which(data_R3$cor == 0),]$FES, ylim=c(0,1000), xlab="Final epidemic size", ylab="Number of epidemics", main="")
abline(v=200, lty="dashed", col="gray70")
hist(data_R3_hom[which(data_R3_hom$cv_s == 1),]$FES)

# summarize data for peak size, peak time, and final epidemic size given a major epidemic and probability of a major epidemic
data_R3_summary <- data_R3 %>% filter(major_out != 0) %>% 
  group_by(across(all_of(c("cv_s", "cv_t", "cor")))) %>%
  summarise_at(c("peak_size", "peak_time", "FES"), mean, na.rm = TRUE)
  
data_R3_hom_summary <- data_R3_hom %>% filter(major_out != 0) %>% 
  group_by(across(all_of(c("cv_s","cv_t", "cor")))) %>%
  summarise_at(c("peak_size", "peak_time", "FES"), mean, na.rm = TRUE)

data_R3_probout <- data_R3 %>% 
  group_by(across(all_of(c("cv_s", "cv_t", "cor")))) %>%
  summarise_at(c("major_out"), mean, na.rm = TRUE)

data_R3_hom_probout <- data_R3_hom %>% 
  group_by(across(all_of(c("cv_s", "cv_t", "cor")))) %>%
  summarise_at(c("major_out"), mean, na.rm = TRUE)


# add back in data points that don't have major epidemics as "NA" for plotting
for(cv_s in c(0.5,1,3))
{
  for(cv_t in c(0.5,1,3))
  {
    for(corr in cors)
    {
      # check if this parameter set is in the summary dataframe, if not (returns TRUE) add with NA for each epi measure
      if(identical(which(data_R3_summary$cv_s == cv_s & data_R3_summary$cv_t == cv_t & data_R3_summary$cor == corr), integer(0)))
      {
        new_row <- data.frame(cv_s=cv_s, cv_t=cv_t, cor=corr, peak_size=NA, peak_time=NA, FES=NA)
        data_R3_summary <- rbind(data_R3_summary, new_row)
      } 
    }
  }
}

for(i in 1:7)
{
  # check if this parameter set is in the summary dataframe, if not (returns TRUE) add with NA for each epi measure
  if(identical(which(data_R3_hom_summary$cv_s == cv_ss[i] & data_R3_hom_summary$cv_t == cv_ts[i]), integer(0)))
  {
    new_row <- data.frame(cv_s=cv_ss[i], cv_t=cv_ts[i], cor="0", peak_size=NA, peak_time=NA, FES=NA)
    data_R3_hom_summary <- rbind(data_R3_hom_summary, new_row)
  } 
}


# create factors to order results for plotting
data_R3_summary$cv_s <- as.factor(data_R3_summary$cv_s)
data_R3_summary$cv_t <- as.factor(data_R3_summary$cv_t)
data_R3_summary$cor <- factor(data_R3_summary$cor, levels = c(-1,-0.5,0,0.5,1))

data_R3_hom_summary$cv_s <- as.factor(data_R3_hom_summary$cv_s)
data_R3_hom_summary$cv_t <- as.factor(data_R3_hom_summary$cv_t)
data_R3_hom_summary$cor <- factor(data_R3_hom_summary$cor, levels = c(-1,-0.5,0,0.5,1))

data_R3_probout$cv_s <- as.factor(data_R3_probout$cv_s)
data_R3_probout$cv_t <- as.factor(data_R3_probout$cv_t)
data_R3_probout$cor <- factor(data_R3_probout$cor, levels = c(-1,-0.5,0,0.5,1))

data_R3_hom_probout$cv_s <- as.factor(data_R3_hom_probout$cv_s)
data_R3_hom_probout$cv_t <- as.factor(data_R3_hom_probout$cv_t)
data_R3_hom_probout$cor <- factor(data_R3_hom_probout$cor, levels = c(-1,-0.5,0,0.5,1))




############################################################################
######################### plot results #####################################
############################################################################
#### plot peak size
peaksize_corn1 <- ggplot(data_R3_summary[which(data_R3_summary$cor == -1),], aes(x=cv_s, y=cv_t, fill=peak_size)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(12,116),  breaks=seq(12,116,8), # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Size") +
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
peaksize_corn1 <- set_panel_size( peaksize_corn1, height=unit(3*2, "cm"),
                                  width=unit(3*2, "cm") )

peaksize_corn0.5 <- ggplot(data_R3_summary[which(data_R3_summary$cor == -0.5),], aes(x=cv_s, y=cv_t, fill=peak_size)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(12,116),  breaks=seq(12,116,8), # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Size") +
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
peaksize_corn0.5 <- set_panel_size( peaksize_corn0.5, height=unit(3*2, "cm"),
                                    width=unit(3*2, "cm") )

peaksize_cor0 <- ggplot(data_R3_summary[which(data_R3_summary$cor == 0),], aes(x=cv_s, y=cv_t, fill=peak_size)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(12,116),  breaks=seq(12,116,8), # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Size") +
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
peaksize_cor0 <- set_panel_size( peaksize_cor0, height=unit(3*2, "cm"),
                                 width=unit(3*2, "cm") )

peaksize_cor0.5 <- ggplot(data_R3_summary[which(data_R3_summary$cor == 0.5),], aes(x=cv_s, y=cv_t, fill=peak_size)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(12,116),  breaks=seq(12,116,8), # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Size") +
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
peaksize_cor0.5 <- set_panel_size( peaksize_cor0.5, height=unit(3*2, "cm"),
                                   width=unit(3*2, "cm") )

peaksize_cor1 <- ggplot(data_R3_summary[which(data_R3_summary$cor == 1),], aes(x=cv_s, y=cv_t, fill=peak_size)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(12,116),  breaks=seq(12,116,8), # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Size") +
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
peaksize_cor1 <- set_panel_size( peaksize_cor1, height=unit(3*2, "cm"),
                                 width=unit(3*2, "cm") )

peaksize_cs0 <- ggplot(data_R3_hom_summary[which(data_R3_hom_summary$cv_s == 0 & data_R3_hom_summary$cv_t != 0),], aes(x=cv_s, y=cv_t, fill=peak_size)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(12,116),  breaks=seq(12,116,8), # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Size") +
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
peaksize_cs0 <- set_panel_size( peaksize_cs0, height=unit(3*2, "cm"),
                                width=unit(1*2, "cm") )

peaksize_cs0ct0 <- ggplot(data_R3_hom_summary[which(data_R3_hom_summary$cv_s == 0 & data_R3_hom_summary$cv_t == 0),], aes(x=cv_s, y=cv_t, fill=peak_size)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(12,116),  breaks=seq(12,116,8), # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Size") +
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
peaksize_cs0ct0 <- set_panel_size( peaksize_cs0ct0, height=unit(1*2, "cm"),
                                   width=unit(1*2, "cm") )


peaksize_ct0 <- ggplot(data_R3_hom_summary[which(data_R3_hom_summary$cv_s != 0 & data_R3_hom_summary$cv_t == 0),], aes(x=cv_s, y=cv_t, fill=peak_size)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(12,116),  breaks=seq(12,116,8), # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Size") +
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
peaksize_ct0 <- set_panel_size( peaksize_ct0, height=unit(1*2, "cm"),
                                width=unit(3*2, "cm") )


# determine legend for combined plots
peaksize_leg <- ggplot(data_R3_summary[which(data_R3_summary$cor == "0"),], aes(x=cv_t, y=cv_s, fill=peak_size)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  #scale_fill_steps(low="white", high="red", limits=c(15,399),  breaks=seq(15,399,24), # R0=3
  scale_fill_steps(low="white", high="red", limits=c(12,116),  breaks=seq(12,116,8), # R0=0.8
  #scale_fill_steps(low="white", high="red", limits=c(12,166),  breaks=seq(12,166,14), # R0=1.1
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Size") +
  xlab(expression(paste("Coefficient of variation of transmission (",italic(C)[italic(t)],")"))) +
  ylab(expression(paste("Coefficient of variation of susceptibility ( ",italic(C)[italic(s)],")"))) +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), legend.key.height=unit(2.5,"cm"))
legend <- get_legend(peaksize_leg)

# plot all combined plots in a grid with common legend
plt <- plot_grid(peaksize_cs0, peaksize_corn1, peaksize_corn0.5, peaksize_cor0, peaksize_cor0.5, peaksize_cor1,
                 NULL, NULL, NULL, NULL, NULL, NULL,
                 peaksize_cs0ct0, NULL, NULL, peaksize_ct0, NULL, NULL,
                 nrow=3, ncol=6, rel_heights=c(1,-0.5,1), rel_widths=c(0.4,1,1,1,1,1))
plot_grid(plt, legend,
          nrow=1, rel_widths=c(1,0.2))
# size: 19 x 7.5



###############################################################################
#### plot peak time
peaktime_corn1 <- ggplot(data_R3_summary[which(data_R3_summary$cor == -1),], aes(x=cv_s, y=cv_t, fill=peak_time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(16,96),  breaks=seq(16,96,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Time") +
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
peaktime_corn1 <- set_panel_size( peaktime_corn1, height=unit(3*2, "cm"),
                                  width=unit(3*2, "cm") )

peaktime_corn0.5 <- ggplot(data_R3_summary[which(data_R3_summary$cor == -0.5),], aes(x=cv_s, y=cv_t, fill=peak_time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(16,96),  breaks=seq(16,96,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Time") +
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
peaktime_corn0.5 <- set_panel_size( peaktime_corn0.5, height=unit(3*2, "cm"),
                                    width=unit(3*2, "cm") )

peaktime_cor0 <- ggplot(data_R3_summary[which(data_R3_summary$cor == 0),], aes(x=cv_s, y=cv_t, fill=peak_time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(16,96),  breaks=seq(16,96,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Time") +
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
peaktime_cor0 <- set_panel_size( peaktime_cor0, height=unit(3*2, "cm"),
                                 width=unit(3*2, "cm") )

peaktime_cor0.5 <- ggplot(data_R3_summary[which(data_R3_summary$cor == 0.5),], aes(x=cv_s, y=cv_t, fill=peak_time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(16,96),  breaks=seq(16,96,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Time") +
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
peaktime_cor0.5 <- set_panel_size( peaktime_cor0.5, height=unit(3*2, "cm"),
                                   width=unit(3*2, "cm") )

peaktime_cor1 <- ggplot(data_R3_summary[which(data_R3_summary$cor == 1),], aes(x=cv_s, y=cv_t, fill=peak_time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(16,96),  breaks=seq(16,96,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Time") +
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
peaktime_cor1 <- set_panel_size( peaktime_cor1, height=unit(3*2, "cm"),
                                 width=unit(3*2, "cm") )

peaktime_cs0 <- ggplot(data_R3_hom_summary[which(data_R3_hom_summary$cv_s == 0 & data_R3_hom_summary$cv_t != 0),], aes(x=cv_s, y=cv_t, fill=peak_time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(16,96),  breaks=seq(16,96,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Time") +
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
peaktime_cs0 <- set_panel_size( peaktime_cs0, height=unit(3*2, "cm"),
                                width=unit(1*2, "cm") )

peaktime_cs0ct0 <- ggplot(data_R3_hom_summary[which(data_R3_hom_summary$cv_s == 0 & data_R3_hom_summary$cv_t == 0),], aes(x=cv_s, y=cv_t, fill=peak_time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(16,96),  breaks=seq(16,96,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Time") +
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
peaktime_cs0ct0 <- set_panel_size( peaktime_cs0ct0, height=unit(1*2, "cm"),
                                   width=unit(1*2, "cm") )


peaktime_ct0 <- ggplot(data_R3_hom_summary[which(data_R3_hom_summary$cv_s != 0 & data_R3_hom_summary$cv_t == 0),], aes(x=cv_s, y=cv_t, fill=peak_time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(16,96),  breaks=seq(16,96,8),  # R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Time") +
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
peaktime_ct0 <- set_panel_size( peaktime_ct0, height=unit(1*2, "cm"),
                                width=unit(3*2, "cm") )


# determine legend for combined plots
peaktime_leg <- ggplot(data_R3_summary[which(data_R3_summary$cor == "0"),], aes(x=cv_t, y=cv_s, fill=peak_time)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  #scale_fill_steps(low="white", high="red", limits=c(8,74),  breaks=seq(8,74,6),  # R0=3
  scale_fill_steps(low="white", high="red", limits=c(16,96),  breaks=seq(16,96,8),  # R0=0.8
  #scale_fill_steps(low="white", high="red", limits=c(15,132),  breaks=seq(15,132,13),  # R0=1.1
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Peak Time") +
  xlab(expression(paste("Coefficient of variation of transmission (",italic(C)[italic(t)],")"))) +
  ylab(expression(paste("Coefficient of variation of susceptibility ( ",italic(C)[italic(s)],")"))) +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), legend.key.height=unit(2.5,"cm"))
legend <- get_legend(peaktime_leg)

# plot all combined plots in a grid with common legend
plt <- plot_grid(peaktime_cs0, peaktime_corn1, peaktime_corn0.5, peaktime_cor0, peaktime_cor0.5, peaktime_cor1,
                 NULL, NULL, NULL, NULL, NULL, NULL,
                 peaktime_cs0ct0, NULL, NULL, peaktime_ct0, NULL, NULL,
                 nrow=3, ncol=6, rel_heights=c(1,-0.5,1), rel_widths=c(0.4,1,1,1,1,1))
plot_grid(plt, legend,
          nrow=1, rel_widths=c(1,0.2))
# size: 19 x 7.5






###################################################################################
#### plot final epidemic size
FES_corn1 <- ggplot(data_R3_summary[which(data_R3_summary$cor == -1),], aes(x=cv_s, y=cv_t, fill=FES)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(53,405),  breaks=seq(53,405,32),  #R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Final Epidemic Size") +
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
FES_corn1 <- set_panel_size( FES_corn1, height=unit(3*2, "cm"),
                             width=unit(3*2, "cm") )

FES_corn0.5 <- ggplot(data_R3_summary[which(data_R3_summary$cor == -0.5),], aes(x=cv_s, y=cv_t, fill=FES)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(53,405),  breaks=seq(53,405,32),  #R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Final Epidemic Size") +
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
FES_corn0.5 <- set_panel_size( FES_corn0.5, height=unit(3*2, "cm"),
                               width=unit(3*2, "cm") )

FES_cor0 <- ggplot(data_R3_summary[which(data_R3_summary$cor == 0),], aes(x=cv_s, y=cv_t, fill=FES)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(53,405),  breaks=seq(53,405,32),  #R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Final Epidemic Size") +
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
FES_cor0 <- set_panel_size( FES_cor0, height=unit(3*2, "cm"),
                            width=unit(3*2, "cm") )

FES_cor0.5 <- ggplot(data_R3_summary[which(data_R3_summary$cor == 0.5),], aes(x=cv_s, y=cv_t, fill=FES)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(53,405),  breaks=seq(53,405,32),  #R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Final Epidemic Size") +
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
FES_cor0.5 <- set_panel_size( FES_cor0.5, height=unit(3*2, "cm"),
                              width=unit(3*2, "cm") )

FES_cor1 <- ggplot(data_R3_summary[which(data_R3_summary$cor == 1),], aes(x=cv_s, y=cv_t, fill=FES)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(53,405),  breaks=seq(53,405,32),  #R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Final Epidemic Size") +
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
FES_cor1 <- set_panel_size( FES_cor1, height=unit(3*2, "cm"),
                            width=unit(3*2, "cm") )

FES_cs0 <- ggplot(data_R3_hom_summary[which(data_R3_hom_summary$cv_s == 0 & data_R3_hom_summary$cv_t != 0),], aes(x=cv_s, y=cv_t, fill=FES)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(53,405),  breaks=seq(53,405,32),  #R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Final Epidemic Size") +
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
FES_cs0 <- set_panel_size( FES_cs0, height=unit(3*2, "cm"),
                           width=unit(1*2, "cm") )

FES_cs0ct0 <- ggplot(data_R3_hom_summary[which(data_R3_hom_summary$cv_s == 0 & data_R3_hom_summary$cv_t == 0),], aes(x=cv_s, y=cv_t, fill=FES)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(53,405),  breaks=seq(53,405,32),  #R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Final Epidemic Size") +
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
FES_cs0ct0 <- set_panel_size( FES_cs0ct0, height=unit(1*2, "cm"),
                              width=unit(1*2, "cm") )


FES_ct0 <- ggplot(data_R3_hom_summary[which(data_R3_hom_summary$cv_s != 0 & data_R3_hom_summary$cv_t == 0),], aes(x=cv_s, y=cv_t, fill=FES)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(53,405),  breaks=seq(53,405,32),  #R0=0.8
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Final Epidemic Size") +
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
FES_ct0 <- set_panel_size( FES_ct0, height=unit(1*2, "cm"),
                           width=unit(3*2, "cm") )


# determine legend for combined plots
FES_leg <- ggplot(data_R3_summary[which(data_R3_summary$cor == "0"),], aes(x=cv_t, y=cv_s, fill=FES)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  #scale_fill_steps(low="white", high="red", limits=c(67,941),  breaks=seq(67,941,46),  #R0=3
  scale_fill_steps(low="white", high="red", limits=c(53,405),  breaks=seq(53,405,32),  #R0=0.8
  #scale_fill_steps(low="white", high="red", limits=c(68,552),  breaks=seq(68,552,44),  #R0=1.1
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="gray") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Final Epidemic Size") +
  xlab(expression(paste("Coefficient of variation of transmission (",italic(C)[italic(t)],")"))) +
  ylab(expression(paste("Coefficient of variation of susceptibility ( ",italic(C)[italic(s)],")"))) +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), legend.key.height=unit(2.5,"cm"))
legend <- get_legend(FES_leg)


# plot all combined plots in a grid with common legend
plt <- plot_grid(FES_cs0, FES_corn1, FES_corn0.5, FES_cor0, FES_cor0.5, FES_cor1,
                 NULL, NULL, NULL, NULL, NULL, NULL,
                 FES_cs0ct0, NULL, NULL, FES_ct0, NULL, NULL,
                 nrow=3, ncol=6, rel_heights=c(1,-0.5,1), rel_widths=c(0.4,1,1,1,1,1))
plot_grid(plt, legend,
          nrow=1, rel_widths=c(1,0.2))
# size: 19 x 7.5





###################################################################################
#### plot probability of a major outbreak
# for R0=0.8 set 0 values as NA to distinguish
data_R3_probout$major_out <- ifelse(data_R3_probout$major_out==0,NA,data_R3_probout$major_out)
data_R3_hom_probout$major_out <- ifelse(data_R3_hom_probout$major_out==0,NA,data_R3_hom_probout$major_out)

probout_corn1 <- ggplot(data_R3_probout[which(data_R3_probout$cor == -1),], aes(x=cv_s, y=cv_t, fill=major_out)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(0,0.694),  breaks=seq(0.006,0.694,0.043),  # R0=3
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="white") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Prob Major Epidemic") +
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
probout_corn1 <- set_panel_size( probout_corn1, height=unit(3*2, "cm"),
                                 width=unit(3*2, "cm") )

probout_corn0.5 <- ggplot(data_R3_probout[which(data_R3_probout$cor == -0.5),], aes(x=cv_s, y=cv_t, fill=major_out)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(0,0.694),  breaks=seq(0.006,0.694,0.043),  # R0=3
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="white") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Prob Major Epidemic") +
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
probout_corn0.5 <- set_panel_size( probout_corn0.5, height=unit(3*2, "cm"),
                                   width=unit(3*2, "cm") )

probout_cor0 <- ggplot(data_R3_probout[which(data_R3_probout$cor == 0),], aes(x=cv_s, y=cv_t, fill=major_out)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(0,0.694),  breaks=seq(0.006,0.694,0.043),  # R0=3
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="white") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Prob Major Epidemic") +
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
probout_cor0 <- set_panel_size( probout_cor0, height=unit(3*2, "cm"),
                                width=unit(3*2, "cm") )

probout_cor0.5 <- ggplot(data_R3_probout[which(data_R3_probout$cor == 0.5),], aes(x=cv_s, y=cv_t, fill=major_out)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(0,0.694),  breaks=seq(0.006,0.694,0.043),  # R0=3
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="white") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Prob Major Epidemic") +
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
probout_cor0.5 <- set_panel_size( probout_cor0.5, height=unit(3*2, "cm"),
                                  width=unit(3*2, "cm") )

probout_cor1 <- ggplot(data_R3_probout[which(data_R3_probout$cor == 1),], aes(x=cv_s, y=cv_t, fill=major_out)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(0,0.694),  breaks=seq(0.006,0.694,0.043),  # R0=3
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="white") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Prob Major Epidemic") +
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
probout_cor1 <- set_panel_size( probout_cor1, height=unit(3*2, "cm"),
                                width=unit(3*2, "cm") )

probout_cs0 <- ggplot(data_R3_hom_probout[which(data_R3_hom_probout$cv_s == 0 & data_R3_hom_probout$cv_t != 0),], aes(x=cv_s, y=cv_t, fill=major_out)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(0,0.694),  breaks=seq(0.006,0.694,0.043),  # R0=3
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="white") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Prob Major Epidemic") +
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
probout_cs0 <- set_panel_size( probout_cs0, height=unit(3*2, "cm"),
                               width=unit(1*2, "cm") )

probout_cs0ct0 <- ggplot(data_R3_hom_probout[which(data_R3_hom_probout$cv_s == 0 & data_R3_hom_probout$cv_t == 0),], aes(x=cv_s, y=cv_t, fill=major_out)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(0,0.694),  breaks=seq(0.006,0.694,0.043),  # R0=3
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="white") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Prob Major Epidemic") +
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
probout_cs0ct0 <- set_panel_size( probout_cs0ct0, height=unit(1*2, "cm"),
                                  width=unit(1*2, "cm") )


probout_ct0 <- ggplot(data_R3_hom_probout[which(data_R3_hom_probout$cv_s != 0 & data_R3_hom_probout$cv_t == 0),], aes(x=cv_s, y=cv_t, fill=major_out)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(0,0.694),  breaks=seq(0.006,0.694,0.043),  # R0=3
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="white") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Prob Major Epidemic") +
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
probout_ct0 <- set_panel_size( probout_ct0, height=unit(1*2, "cm"),
                               width=unit(3*2, "cm") )


# determine legend for combined plots
probout_leg <- ggplot(data_R3_probout[which(data_R3_probout$cor == "0"),], aes(x=cv_t, y=cv_s, fill=major_out)) +
  geom_tile(color="gray40", lwd=1, linetype=1) +
  #scale_fill_gradient(low="white", high="red", breaks=seq(15,365,50)) +
  scale_fill_steps(low="white", high="red", limits=c(0,0.694),  breaks=seq(0.006,0.694,0.043),  # R0=3
  #scale_fill_steps(low="white", high="red", limits=c(0.0,0.265),  breaks=seq(0.001,0.265,0.024),  # R0=0.8
  #scale_fill_steps(low="white", high="red", limits=c(0.0,0.353),  breaks=seq(0.001,0.353,0.032),  # R0=1.1
                   guide=guide_colorbar(frame.colour="gray40", ticks.colour="gray40"), na.value="white") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="Prob Major Epidemic") +
  xlab(expression(paste("Coefficient of variation of transmission (",italic(C)[italic(t)],")"))) +
  ylab(expression(paste("Coefficient of variation of susceptibility ( ",italic(C)[italic(s)],")"))) +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = "gray"), panel.border = element_rect(color="black", fill=NA),
        text=element_text(size=20), legend.text=element_text(size=18), legend.title = element_text(size=20),
        axis.text.x=element_text(hjust=0.5), axis.ticks=element_line(),
        legend.key=element_rect(), legend.key.height=unit(2.5,"cm"))
legend <- get_legend(probout_leg)


# plot all combined plots in a grid with common legend
plt <- plot_grid(probout_cs0, probout_corn1, probout_corn0.5, probout_cor0, probout_cor0.5, probout_cor1,
                 NULL, NULL, NULL, NULL, NULL, NULL,
                 probout_cs0ct0, NULL, NULL, probout_ct0, NULL, NULL,
                 nrow=3, ncol=6, rel_heights=c(1,-0.5,1), rel_widths=c(0.4,1,1,1,1,1))
plot_grid(plt, legend,
          nrow=1, rel_widths=c(1,0.2))
# size: 19 x 7.5
