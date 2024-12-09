# set up risk and force of infection distributions for stochastic individual-based SIR model
# set each as gamma distributions with normal copula to set level of correlation

install.packages("copula", lib="../R_libs", repos="http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
library("copula", lib="../R_libs")

## read in population size, CVs for dists, level of correlation, seed
arg <- commandArgs(trailingOnly=T)

# population size
pop_size <- as.numeric(arg[1])

# coefficients of variation for susceptibility (cv_s) and transmission (cv_t)
cv_s <- as.numeric(arg[2])
cv_t <- as.numeric(arg[3])

# level of correlation
cor_level <- as.numeric(arg[4])

# seed
seed <- as.numeric(arg[5])


## set shape parameters for gamma dists from CVs, scale is 1/shape
k <- 1 / (cv_s)^2
m <- 1 / (cv_t)^2


## set up and draw from distribution
set.seed(seed)

# set up copula with desired correlation level
my_copula <- normalCopula(cor_level)

# set up gamma dists with desired parameters and copula
my_dist <- mvdc(my_copula, margins=rep("gamma",2),
                paramMargins=list(list(shape=k,scale=1/k),list(shape=m,scale=1/m)))

# randomly draw risk and force of infection values from distribution
risk_FoIs <- rMvdc(pop_size, my_dist)


## read out risk and force of infection values to file to read into SIR model
filename <- paste("risksFoIs_his",cv_s,"_hit",cv_t,"_cor",cor_level,"_p",pop_size,"_seed",seed,".txt",sep="")
write.table(risk_FoIs, file=filename, sep="\t", row.names=F, col.names=F)
