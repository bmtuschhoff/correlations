// individual based SIR Gillespie model including different forms of heterogeneity
// heterogeneity can be in contact rate or given contact and it can be HiT or HiS
// also set up different correlations between HiT and HiS (positive, negative, uncorrelated)
// want to know if dynamics look different for different types of heterogeneity and what can be distinguished

#include <iostream>
#include <cmath>		// sqrt and log
#include <math.h>
#include <gsl/gsl_rng.h>	// to set seeds for random dists
#include <gsl/gsl_randist.h>	// random dists
#include <fstream>		// file input/output
#include <random>		// set up distributions of rates to choose indiv infected/recovered
#include <cstring>		// to compare strings for inputting parameters from command line
#include <algorithm> 		// contains sort() function
#include <sstream>		// to read in risks and FoIs
#include <string>		// to read in risks and FoIs

using namespace std;
using std::istringstream;
using std::string;


bool parseOptions(int argc, char *argv[]);
double find_tot_rate(double *infect_rate, double *recover_rate, char current_pop_status[], vector<double> risk_S, vector<double> c_rate_S, vector<double> lamb_I, vector<double> c_rate_I, double beta_avg, int numS, int numI);
void infect(char current_pop_status[], vector<double> &risk_S, vector<double> &lamb_S, vector<double> &c_rate_S, vector<int> &index_S, vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, double beta_avg, int *numS, int *numI, mt19937 &rng2);
void recover(char current_pop_status[], vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, int *numI, mt19937 &rng2);
double c_rate(int S_index, int I_index);


// default parameter settings, can change with flags
string output = "results.txt";	// output file name
string input = "input.txt";	// input file name
char his = 'F';			// set whether to include HiS given contact as false
char hit = 'F';			// set whether to include HiT given contact as false
char hic = 'F';			// set whether to include het in contact as false
double cor = 0.0;		// set correlation between HiS and HiT -- could be anything in [-1,1]
int pop_size = 4;		// population size
int init_inf = 20;		// initial number of infected indivs
double R0 = 3;			// R0
double gam = 0.1;		// recovery rate
double cv_s = 1.3;		// coefficient of variation of risk (HiS)
double cv_t = 1.0 / sqrt(0.5);		// coefficient of variation (HiT)
double cv_c = 1.0;		// coefficient of variation for contact rate
int seed = 3;			// seed for random number generators


int main(int argc, char *argv[])
{
	parseOptions(argc, argv);	// set up parameters from command line input and flags

	double k, m, psi;	// shape params for HiS, HiT, HiC distributions
	double beta_avg;	// avg transmission rate
	int time = 1000;	// total length of time to run the simulation
	int numS, numI;		// number of S and I indivs in the population
	//int init_inf = 20;	// initial number of infected indivs
	vector<double> risk_S(pop_size), lamb_S(pop_size), c_rate_S(pop_size);	// vectors of risks, forces of inf, and contact rates for each S indiv
	vector<double> risk_I(pop_size), lamb_I(pop_size), c_rate_I(pop_size);	// vectors of risks, forces of inf, and contact rates for each I indiv
	vector<int> index_S(pop_size), index_I(pop_size);	// vectors of indices for each S, I indiv
	int indiv, i, j;	// variables for iteration
	int index;		// index of node in population
	char current_pop_status[pop_size];	// vector holding current status of each indiv (S,I,R)
	int init_inf_indivs[init_inf];		// vector holding the indices of the indivs initially infected
	int pop_indices[pop_size];		// vector of indices of all indivs in the pop
	double t = 0.0;		// track time throughout epidemic
	double infect_rate, recover_rate, tot_rate;	// rates of infection, recovery, and both across all indivs
	double event;	// uniform RV to determine which event occurs next
	double Reff;	// reproductive number R effective
	int inf_num = 0;	// number of infections that have occurred
	string line;	// each line read in from input file



	// declare file input stream and open
	ifstream ifs;
	ifs.open(input.c_str());

	// read all risks and FoIs from input file
	// get one line at a time and read the two numbers on each line
	// already know there will be pop_size number of lines with risk and FoI on each line
	for(indiv = 0; indiv < pop_size; indiv++)
	{
		getline(ifs, line);	// get a line from input file and save in "line"
		istringstream line_stream(line);	// declare line stream to get numbers from that line
		line_stream >> risk_S[indiv] >> lamb_S[indiv];	// read risk and FoI from line into appropriate variables
	}

	// close input file
	ifs.close();


	// calculate other params from those read-in
	psi = 1.0 / pow(cv_c, 2);

	// calculate avg transmission rate from R0, gam, pop size
	beta_avg = (gam * R0) / (pop_size - init_inf);

	// set up output file stream
	// check if file already exists
	bool file_exists = false;
	if(FILE *file = fopen(output.c_str(), "r"))
	{
		fclose(file);
		file_exists = true;
	}

	// declare file output stream
	ofstream ofs;
	ofs.open(output.c_str(), ofstream::app);	// set up to append new data to end of file

	// write header to file if it doesn't already exist
	if(file_exists == false)
	{
		ofs << "time" << "\tS" << "\tI" << "\tR" << "\tReff" << "\tinf_num" << endl;
	}
	


	// set up random number generator for using with GSL functions - GSL's Taus generator
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
	// initialize generator with seed
	gsl_rng_set(rng, seed);
	// set up second random number generator for using with discrete_distribution
	mt19937 rng2(seed);


	// draw contact rate for each indiv from gamma dist depending on whether including that het
	for(indiv = 0; indiv < pop_size; indiv++)
	{
		if(hic == 'T')	// there is HiC
			c_rate_S[indiv] = gsl_ran_gamma(rng, psi, 1.0/psi);
		else	// no HiC
			c_rate_S[indiv] = 1.0;
	}



	// initialize vector with current status of each indiv
	// randomly select init_inf indiv(s) to be infected at the start
	for(i = 0; i < pop_size; i++)	// initialize vector of indices for each indiv in pop, for randomly choosing indivs
	{				// also initialize all indivs in pop as S
		pop_indices[i] = i;
		index_S[i] = i;
		current_pop_status[i] = 'S';
	}	

	gsl_ran_choose(rng, init_inf_indivs, init_inf, pop_indices, pop_size, sizeof(int));	// randomly select the indices of init_inf indivs out of pop_size to start off as inf

	for(i = 0; i < init_inf; i++)	// set up indivs to be infected at start
		current_pop_status[init_inf_indivs[i]] = 'I';

	// move each initially infected indiv from S to I vectors
	for(i = 0; i < init_inf; i++)
	{
		// move each param from S to first available row in I (all rows available currently)
		risk_I[i] = risk_S[init_inf_indivs[i]];
		lamb_I[i] = lamb_S[init_inf_indivs[i]];
		c_rate_I[i] = c_rate_S[init_inf_indivs[i]];
		index_I[i] = index_S[init_inf_indivs[i]];

		// replace data for each param in S with data from last row of S
		risk_S[init_inf_indivs[i]] = risk_S[pop_size-1-i];
		lamb_S[init_inf_indivs[i]] = lamb_S[pop_size-1-i];
		c_rate_S[init_inf_indivs[i]] = c_rate_S[pop_size-1-i];
		index_S[init_inf_indivs[i]] = index_S[pop_size-1-i];
	}



	// set num S and num I indivs at start of epidemic
	numS = pop_size - init_inf;
	numI = init_inf;
	
	
	// write out data at time 0: time and current compartment status of each indiv in population
	ofs << t << "\t" << numS << "\t" << numI << "\t" << pop_size-numS-numI << "\t" << R0 << "\t" << inf_num << endl;

	// start Gillespie algorithm - draw time, determine which event happened, carry out event, repeat
	while(t <= time)
	{

//cout << "numS: " << numS << " numI: " << numI << endl;

		// determine rate of each event and combined total rate -- events are infection or recovery
		tot_rate = find_tot_rate(&infect_rate, &recover_rate, current_pop_status, risk_S, c_rate_S, lamb_I, c_rate_I, beta_avg, numS, numI);

		if(tot_rate == -1)	// no inf indivs, so epidemic is over - close output file and end simulation
		{
			ofs.close();
			return(0);
		}

		// determine time of next event
		t += gsl_ran_exponential(rng, 1 / tot_rate);	// rate lambda for exp dist is tot_rate, this fntn takes the mean 1/lambda

		// determine next event and carry out
		// draw uniform RV and compare to probabilities of events
		event = gsl_ran_flat(rng, 0, 1);
		
		
		if(event <= infect_rate / tot_rate)	// event is infection
		{
			infect(current_pop_status, risk_S, lamb_S, c_rate_S, index_S, risk_I, lamb_I, c_rate_I, index_I, beta_avg, &numS, &numI, rng2);
			inf_num++;
		}
		else	// event is recovery
		{
			recover(current_pop_status, risk_I, lamb_I, c_rate_I, index_I, &numI, rng2);
		}

		// calculate Reff at this time point
		tot_rate = find_tot_rate(&infect_rate, &recover_rate, current_pop_status, risk_S, c_rate_S, lamb_I, c_rate_I, beta_avg, numS, numI);	// update infect_rate
		Reff = (infect_rate / (numI)) * (1.0 / gam);

		// write out time and current compartment status of each indiv in population
		ofs << t << "\t" << numS << "\t" << numI << "\t" << pop_size-numS-numI << "\t" << Reff << "\t" << inf_num << endl;
	}

	ofs.close();	// close ouput file
	return(0);
}




// find total rates of next event
double find_tot_rate(double *infect_rate, double *recover_rate, char current_pop_status[], vector<double> risk_S, vector<double> c_rate_S, vector<double> lamb_I, vector<double> c_rate_I, double beta_avg, int numS, int numI)
{
	int i, j;	// variables for iteration
	double prob_inf;	// probability of infection given contact
	
	*infect_rate = 0.0;
	*recover_rate = 0.0;
	for(i = 0; i < numI; i++)	// calculate infection and recovery rates across all infected indivs
	{

		*recover_rate += gam;	// each inf indiv could recover with rate gamma

		for(j = 0; j < numS; j++)	// determine infection rate across all S indivs that could have contacted this indiv
		{
			*infect_rate += beta_avg * c_rate(j, i) * risk_S[j] * lamb_I[i];
		}
	}

	if(numI == 0)	// no one is infected, epidemic is over, flag total rate
	{
		return(-1);
	}

	return(*infect_rate + *recover_rate);	// total rate of next event
}




// know event is that an indiv is inf, need to decide which one and implement
void infect(char current_pop_status[], vector<double> &risk_S, vector<double> &lamb_S, vector<double> &c_rate_S, vector<int> &index_S, vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, double beta_avg, int *numS, int *numI, mt19937 &rng2)
{
	int i, j;	// variables for iterations
	double prob_inf;	// probability of infection given contact
	int I_indiv;	// index of indiv that will be infected, index in terms of indivs currently S/S vector
	vector<double> infect_rate(*numS, 0);	// vector of infection rate for each contact
	
	for(i = 0; i < *numI; i++)
	{
		for(j = 0; j < *numS; j++)	// determine infection rate across all S indivs that could have contacted this indiv
		{
			infect_rate[j] += beta_avg * c_rate(j, i) * risk_S[j] * lamb_I[i];
		}
	}


	// have vector that contains infection rate for each contact
	// need to pick which indiv is inf
	
	// choose which contact is infected
	discrete_distribution<int> pick_inf(infect_rate.begin(), infect_rate.end());	// set up dist of infection rates to choose who is inf
	I_indiv = pick_inf(rng2);	// pick index of indiv who is infected

	current_pop_status[index_S[I_indiv]] = 'I';	// change status of chosen indiv to be inf

	// move each param from S to first available row in I - since there are numI I indivs, the first row would be at index numI
	risk_I[*numI] = risk_S[I_indiv];
	lamb_I[*numI] = lamb_S[I_indiv];
	c_rate_I[*numI] = c_rate_S[I_indiv];
	index_I[*numI] = index_S[I_indiv];

	// replace data for each param in S with data from last row of numS S indivs
	risk_S[I_indiv] = risk_S[*numS-1];
	lamb_S[I_indiv] = lamb_S[*numS-1];
	c_rate_S[I_indiv] = c_rate_S[*numS-1];
	index_S[I_indiv] = index_S[*numS-1];

	// keep track of number of S and I indivs
	--*numS;
	++*numI;
}




// know event is that an indiv recovers, need to decide which one and implement
void recover(char current_pop_status[], vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, int *numI, mt19937 &rng2)
{
	int R_indiv;	// index of indiv that will recover, index in terms of indivs currently I/I vector

	// all indivs have same recovery rate gam so just have to pick number from 0 to numI-1
	uniform_int_distribution<int> pick_recover(0, *numI-1);
	R_indiv = pick_recover(rng2);	// pick index of indiv who recovers

	current_pop_status[index_I[R_indiv]] = 'R';	// change status of chosen indiv to be recovered

	// replace data for each param in I with data from last row of numI I indivs since indiv is no longer I
	risk_I[R_indiv] = risk_I[*numI-1];
	lamb_I[R_indiv] = lamb_I[*numI-1];
	c_rate_I[R_indiv] = c_rate_I[*numI-1];
	index_I[R_indiv] = index_I[*numI-1];

	// keep track of number of I indivs
	--*numI;
}



double c_rate(int S_index, int I_index)
{
	return(1.0);
}



// read in parameters and parse
bool parseOptions(int argc, char *argv[])
{
	for(int a = 1; a < argc; a++)
	{
		if(strcmp(argv[a], "-o") == 0 && a + 1 < argc)	// set output file name
			output = argv[++a];
		if(strcmp(argv[a], "-in") == 0 && a + 1 < argc)	// set input file name
			input = argv[++a];
		if(strcmp(argv[a], "-his") == 0)	// include heterogeneity in susceptibility given contact
			his = 'T';
		if(strcmp(argv[a], "-hit") == 0)	// include heterogeneity in transmission given contact
			hit = 'T';
		if(strcmp(argv[a], "-hic") == 0)	// include heterogeneity in contact rate
			hic = 'T';
		if(strcmp(argv[a], "-cor") == 0 && a + 1 < argc)	// set correlation between HiS and HiT
			cor = atof(argv[++a]);
		if(strcmp(argv[a], "-p") == 0 && a + 1 < argc)	// set population size
			pop_size = atoi(argv[++a]);
		if(strcmp(argv[a], "-i") == 0 && a + 1 < argc)	// set initial number of infected individuals
			init_inf = atoi(argv[++a]);
		if(strcmp(argv[a], "-r") == 0 && a + 1 < argc)	// set R0
			R0 = atof(argv[++a]);
		if(strcmp(argv[a], "-g") == 0 && a + 1 < argc)	// set recovery rate gamma
			gam = atof(argv[++a]);
		if(strcmp(argv[a], "-cs") == 0 && a + 1 < argc)	// set cv for HiS
			cv_s = atof(argv[++a]);
		if(strcmp(argv[a], "-ct") == 0 && a + 1 < argc)	// set cv for HiT
			cv_t = atof(argv[++a]);
		if(strcmp(argv[a], "-s") == 0 && a + 1 < argc)	// set seed for simulation
			seed = atoi(argv[++a]);
	}

	return(false);
}

