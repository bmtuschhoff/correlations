The code in this directory was used to simulate data for mpox dynamics using an SEIR model. Below are brief descriptions of the files.

GenRisksFoIs_CorHiST_copula.R generates a risk of being infected and force of infection for each individual in the population with correlations between these parameters set using a normal copula with marginal gamma distributions.
The output is a file containing the risk and force of infection for each individual.

GenRisksFoIs_CorHiST_copulab.R generates a risk of being infected and force of infection for many individuals that could enter the population with correlations between these parameters set using a normal copula with marginal gamma distributions.
The output is a file containing the risk and force of infection for each individual.

IBM_SEIR_CorHiST.cpp simulates a stochastic, individual-based SEIR model via the Gillespie algorithm including both heterogeneity in transmission and susceptibility that flexibly allows for positive or negative correlations between transmissibility and susceptibility.
This model also allows for individuals to enter and leave the population.
The files GenRisksFoIs_CorHiST_copula.R and GenRisksFoIs_CorHiST_copulab.R must first be run to generate the risks and forces of infection for each individual.
The output is a file with the number of S, E, I, and R individuals at each time point in the epidemic.

RunSims_mpox.sb is the script used to set all parameters and run the simulations over different parameter combinations. This script was set up to submit jobs to a cluster that uses Slurm for job scheduling and management.