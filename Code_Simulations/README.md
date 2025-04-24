The code in this directory was used to simulate data for the theoretical model. Code in the Mpox directory was used to simulate data for mpox dynamics. Below are brief descriptions of the files in this directory.

GenRisksFoIs_CorHiST_copula.R generates a risk of being infected and force of infection for each individual in the population with correlations between these parameters set using a normal copula with marginal gamma distributions.
The output is a file containing the risk and force of infection for each individual.

IBM_SIR_CorHiST.cpp simulates a stochastic, individual-based SIR model via the Gillespie algorithm including both heterogeneity in transmission and susceptibility that flexibly allows for positive or negative correlations between transmissibility and susceptibility.
The file GenRisksFoIs_CorHiST_copula.R must first be run to generate the risks and forces of infection for each individual.
The output is a file with the number of S, I, and R individuals at each time point in the epidemic.

RunSims.sb is the script used to set all parameters and run the simulations over different parameter combinations. This script was set up to submit jobs to a cluster that uses Slurm for job scheduling and management.
