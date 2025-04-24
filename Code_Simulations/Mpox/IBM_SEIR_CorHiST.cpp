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
#include <cstdio>		// for getting risks and FoIs for new indivs
#include <memory>		// for getting risks and FoIs for new indivs
#include <chrono>		// track run time for events


using namespace std;
using std::istringstream;
using std::string;


bool parseOptions(int argc, char *argv[]);
double find_tot_rate(double *birth_rate, double *infect_rate, double *EtoI_rate, double *recover_rate, double *death_rate, char current_pop_status[], vector<double> risk_S, vector<double> c_rate_S, vector<double> lamb_I, vector<double> c_rate_I, double beta_avg, int numS, int numE, int numI, int numR);
void birth(vector<double> &risk_birth, vector<double> &lamb_birth, int *num_birth, vector<double> &risk_S, vector<double> &lamb_S, vector<double> &c_rate_S, vector<int> &index_S, vector<double> &lamb_I, int *numS, int numE, int numI, int numR, double beta_avg, double *infect_rate);
void infect(char current_pop_status[], vector<double> &risk_S, vector<double> &lamb_S, vector<double> &c_rate_S, vector<int> &index_S, vector<double> &risk_E, vector<double> &lamb_E, vector<double> &c_rate_E, vector<int> &index_E, vector<double> &lamb_I, double beta_avg, int *numS, int *numE, int numI, double *infect_rate, mt19937 &rng2);
void EtoI(char current_pop_status[], vector<double> &risk_E, vector<double> &lamb_E, vector<double> &c_rate_E, vector<int> &index_E, vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, vector<double> &risk_S, int numS, int *numE, int *numI, double beta_avg, double *infect_rate, mt19937 &rng2);
void recover(char current_pop_status[], vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, vector<double> &risk_S, int numS, int *numI, int *numR, double beta_avg, double *infect_rate, mt19937 &rng2);
void death(vector<double> &risk_S, vector<double> &lamb_S, vector<double> &c_rate_S, vector<int> &index_S, vector<double> &risk_E, vector<double> &lamb_E, vector<double> &c_rate_E, vector<int> &index_E, vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, int *numS, int *numE, int *numI, int *numR, double beta_avg, double *infect_rate, mt19937 &rng2);
double c_rate(int S_index, int I_index);
string exec_Rscript(const char *runR);


// default parameter settings, can change with flags
string output = "results.txt";	// output file name
string input = "input.txt";	// input file name
string inputb = "inputb.txt";	// input file name for births
char his = 'F';			// set whether to include HiS given contact as false
char hit = 'F';			// set whether to include HiT given contact as false
char hic = 'F';			// set whether to include het in contact as false
double cor = 0.0;		// set correlation between HiS and HiT -- could be anything in [-1,1]
int pop_size = 40;		// population size
int init_inf = 1;		// initial number of infected indivs
double R0 = 3;			// R0
double gam = 0.1;		// recovery rate
double delt = 0.1;		// rate at which individuals move from E to I
double b = 0.1;			// birth rate
double d = 0.1;			// death rate
double cv_s = 1.3;		// coefficient of variation of risk (HiS)
double cv_t = 1.0 / sqrt(0.5);		// coefficient of variation (HiT)
double cv_c = 1.0;		// coefficient of variation for contact rate
int seed = 3;			// seed for random number generators


int main(int argc, char *argv[])
{
	parseOptions(argc, argv);	// set up parameters from command line input and flags

	double k, m, psi;	// shape params for HiS, HiT, HiC distributions
	double beta_avg;	// avg transmission rate
	int time = 365*3;	// total length of time to run the simulation
	int numS, numE, numI, numR;		// number of S, E, and I indivs in the population
	int num_birth = 0;	// number of individuals born
	//int init_inf = 20;	// initial number of infected indivs
	vector<double> risk_birth(b*365*3), lamb_birth(b*365*3);	// vectors of risks and forces of inf for all indivs to be born
	vector<double> risk_S(pop_size), lamb_S(pop_size), c_rate_S(pop_size);	// vectors of risks, forces of inf, and contact rates for each S indiv
	vector<double> risk_E(pop_size), lamb_E(pop_size), c_rate_E(pop_size);	// vectors of risks, forces of inf, and contact rates for each E indiv
	vector<double> risk_I(pop_size), lamb_I(pop_size), c_rate_I(pop_size);	// vectors of risks, forces of inf, and contact rates for each I indiv
	vector<int> index_S(pop_size), index_E(pop_size), index_I(pop_size);	// vectors of indices for each S, E, I indiv
	int indiv, i, j;	// variables for iteration
	int index;		// index of node in population
	char current_pop_status[pop_size];	// vector holding current status of each indiv (S,I,R)
	int init_inf_indivs[init_inf];		// vector holding the indices of the indivs initially infected
	int pop_indices[pop_size];		// vector of indices of all indivs in the pop
	double t = 0.0;		// track time throughout epidemic
	double birth_rate, infect_rate, EtoI_rate, recover_rate, death_rate, tot_rate;	// rates of infection, progressing from exposed to infectious, recovery, and all across all indivs
	double event;	// uniform RV to determine which event occurs next
	double Reff;	// reproductive number R effective
	int inf_num = 0;	// number of infections that have occurred
	string line;	// each line read in from input file
	int num_events = 0;	// track number of events that have occurred



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


	// open file input stream for births
	ifs.open(inputb.c_str());

	// read all risks and FoIs from input file to be used for all births
	// get one line at a time and read the two numbers on each line
	// already know there will be b*365*3 number of lines with risk and FoI on each line
	for(indiv = 0; indiv < b*365*3; indiv++)
	{
		getline(ifs, line);	// get a line from input file and save in "line"
		istringstream line_stream(line);	// declare line stream to get numbers from that line
		line_stream >> risk_birth[indiv] >> lamb_birth[indiv];	// read risk and FoI from line into appropriate variables
	}

	// close input file
	ifs.close();



	// calculate other params from those read-in
	psi = 1.0 / pow(cv_c, 2);

	// calculate avg transmission rate from R0, gam, pop size
	beta_avg = ((delt + d) * (gam + d) * R0) / (delt * pop_size);

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
		ofs << "time" << "\tS" << "\tE" << "\tI" << "\tR" << "\tReff" << "\tinf_num" << endl;
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
		current_pop_status[init_inf_indivs[i]] = 'E';

	// move each initially infected indiv from S to E vectors
	for(i = 0; i < init_inf; i++)
	{
		// move each param from S to first available row in I (all rows available currently)
		risk_E[i] = risk_S[init_inf_indivs[i]];
		lamb_E[i] = lamb_S[init_inf_indivs[i]];
		c_rate_E[i] = c_rate_S[init_inf_indivs[i]];
		index_E[i] = index_S[init_inf_indivs[i]];

		// replace data for each param in S with data from last row of S
		risk_S[init_inf_indivs[i]] = risk_S[pop_size-1-i];
		lamb_S[init_inf_indivs[i]] = lamb_S[pop_size-1-i];
		c_rate_S[init_inf_indivs[i]] = c_rate_S[pop_size-1-i];
		index_S[init_inf_indivs[i]] = index_S[pop_size-1-i];
	}



	// set num S, E, I, and R indivs at start of epidemic
	numS = pop_size - init_inf;
	numE = init_inf;
	numI = 0;
	numR = 0;

	// calculate infection rate at start of epidemic for baseline
	infect_rate = 0.0;
	for(i = 0; i < numI; i++)	// calculate infection rates across all infected indivs
	{

		for(j = 0; j < numS; j++)	// determine infection rate across all S indivs that could have contacted this indiv
		{
			infect_rate += beta_avg * c_rate(j, i) * risk_S[j] * lamb_I[i];
		}
	}

	
	// write out data at time 0: time and current compartment status of each indiv in population
	ofs << t << "\t" << numS << "\t" << numE << "\t" << numI << "\t" << numR << "\t" << R0 << "\t" << inf_num << endl;


	// start Gillespie algorithm - draw time, determine which event happened, carry out event, repeat
	while(t <= time)
	{

		// determine rate of each event and combined total rate -- events are infection or recovery
		tot_rate = find_tot_rate(&birth_rate, &infect_rate, &EtoI_rate, &recover_rate, &death_rate, current_pop_status, risk_S, c_rate_S, lamb_I, c_rate_I, beta_avg, numS, numE, numI, numR);

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

	

		if(event <= birth_rate / tot_rate)	// event is birth
		{
			birth(risk_birth, lamb_birth, &num_birth, risk_S, lamb_S, c_rate_S, index_S, lamb_I, &numS, numE, numI, numR, beta_avg, &infect_rate);
		}
		else
		{
			if(event <= (birth_rate + infect_rate) / tot_rate)	// event is infection
			{
				infect(current_pop_status, risk_S, lamb_S, c_rate_S, index_S, risk_E, lamb_E, c_rate_E, index_E, lamb_I, beta_avg, &numS, &numE, numI, &infect_rate, rng2);
				inf_num++;
			//	num_events++;
			}
			else
			{
				if(event <= (birth_rate + infect_rate + EtoI_rate) / tot_rate)	// event is individual moving from exposed to infectious class
				{
					EtoI(current_pop_status, risk_E, lamb_E, c_rate_E, index_E, risk_I, lamb_I, c_rate_I, index_I, risk_S, numS, &numE, &numI, beta_avg, &infect_rate, rng2);
				}
				else	// event is recovery
				{
					if(event <= (birth_rate + infect_rate + EtoI_rate + recover_rate) / tot_rate)	// event is individual moving from exposed to infectious class
					{
						recover(current_pop_status, risk_I, lamb_I, c_rate_I, index_I, risk_S, numS, &numI, &numR, beta_avg, &infect_rate, rng2);
					}
					else	// event is death
					{
						death(risk_S, lamb_S, c_rate_S, index_S, risk_E, lamb_E, c_rate_E, index_E, risk_I, lamb_I, c_rate_I, index_I, &numS, &numE, &numI, &numR, beta_avg, &infect_rate, rng2);
					}
				}
			}
		}

		// calculate Reff at this time point
		//tot_rate = find_tot_rate(&birth_rate, &infect_rate, &EtoI_rate, &recover_rate, current_pop_status, risk_S, c_rate_S, lamb_I, c_rate_I, beta_avg, numS, numE, numI, numR);	// update infect_rate
		Reff = (infect_rate / numI) * (1.0 / gam);

		// write out time and current compartment status of each indiv in population
		ofs << t << "\t" << numS << "\t" << numE << "\t" << numI << "\t" << numR << "\t" << Reff << "\t" << inf_num << endl;
	}

	ofs.close();	// close ouput file
	return(0);
}




// find total rates of next event
double find_tot_rate(double *birth_rate, double *infect_rate, double *EtoI_rate, double *recover_rate, double *death_rate, char current_pop_status[], vector<double> risk_S, vector<double> c_rate_S, vector<double> lamb_I, vector<double> c_rate_I, double beta_avg, int numS, int numE, int numI, int numR)
{
	int i, j;	// variables for iteration
	double prob_inf;	// probability of infection given contact
	
	*birth_rate = b;	// b is a "birth" rate not dependent on current pop size
	*EtoI_rate = numE * delt;	// each exposed individual moves to infected with rate delta
	*recover_rate = numI * gam;	// each inf indiv could recover with rate gamma
	*death_rate = d * (numS + numE + numI + numR);	// d is "death" rate per indiv

	if(numE + numI == 0)	// no one is infected, epidemic is over, flag total rate
	{
		return(-1);
	}

	return(*birth_rate + *infect_rate + *EtoI_rate + *recover_rate + *death_rate);	// total rate of next event
}




// know event is that an indiv enters the population, always enter as S
void birth(vector<double> &risk_birth, vector<double> &lamb_birth, int *num_birth, vector<double> &risk_S, vector<double> &lamb_S, vector<double> &c_rate_S, vector<int> &index_S, vector<double> &lamb_I, int *numS, int numE, int numI, int numR, double beta_avg, double *infect_rate)
{
	int i;	// param for iteration

	// assign new indiv's params to vectors in the S population
	// if numS < size of vector for S population, can put new indiv at index numS
	if(*numS < risk_S.size())
	{
		risk_S[*numS] = risk_birth[*num_birth];
		lamb_S[*numS] = lamb_birth[*num_birth];
		c_rate_S[*numS] = 1.0;
		index_S[*numS] = *numS + numE + numI + numR;
	}
	else	// if numS == size of vector for S population, need to append new indiv to the end of vector
	{
		risk_S.push_back(risk_birth[*num_birth]);
		lamb_S.push_back(lamb_birth[*num_birth]);
		c_rate_S.push_back(1.0);
		index_S.push_back(*numS + numE + numI + numR);
	}

	// update infect_rate: S indiv is added so need to multiply their risk with all current I indivs
	for(i = 0; i < numI; i++)
		*infect_rate += beta_avg * c_rate(*numS + numE + numI + numR, i) * risk_birth[*num_birth] * lamb_I[i];
	
	++*numS;
	++*num_birth;
}




// know event is that an indiv is inf, need to decide which one and implement
void infect(char current_pop_status[], vector<double> &risk_S, vector<double> &lamb_S, vector<double> &c_rate_S, vector<int> &index_S, vector<double> &risk_E, vector<double> &lamb_E, vector<double> &c_rate_E, vector<int> &index_E, vector<double> &lamb_I, double beta_avg, int *numS, int *numE, int numI, double *infect_rate, mt19937 &rng2)
{
	int i, j;	// variables for iterations
	double prob_inf;	// probability of infection given contact
	int I_indiv;	// index of indiv that will be infected, index in terms of indivs currently S/S vector
	vector<double> susc(*numS, 0);	// vector of susceptibility for each contact
	
	for(j = 0; j < *numS; j++)	// determine susceptibility for each S indiv
	{
		susc[j] = risk_S[j];
	}


	// have vector that contains infection rate for each contact
	// need to pick which indiv is inf
	
	// choose which contact is infected
	discrete_distribution<int> pick_inf(susc.begin(), susc.end());	// set up dist of infection rates to choose who is inf
	I_indiv = pick_inf(rng2);	// pick index of indiv who is infected

	// move each param from S to first available row in E - there are numE E indivs
	// if numE < size of vector for E population, can put new indiv at index numE
	if(*numE < risk_E.size())
	{
		risk_E[*numE] = risk_S[I_indiv];
		lamb_E[*numE] = lamb_S[I_indiv];
		c_rate_E[*numE] = c_rate_S[I_indiv];
		index_E[*numE] = index_S[I_indiv];
	}
	else	// if numE == size of vector for E population, need to append new indiv to the end of vector
	{
		risk_E.push_back(risk_S[I_indiv]);
		lamb_E.push_back(lamb_S[I_indiv]);
		c_rate_E.push_back(c_rate_S[I_indiv]);
		index_E.push_back(index_S[I_indiv]);
	}

	// update infect_rate: S indiv is removed so need to multiply their risk with all current I indivs and subtract
	for(i = 0; i < numI; i++)
		*infect_rate -= beta_avg * c_rate(index_S[I_indiv], i) * risk_S[I_indiv] * lamb_I[i];

	// replace data for each param in S with data from last row of numS S indivs
	risk_S[I_indiv] = risk_S[*numS-1];
	lamb_S[I_indiv] = lamb_S[*numS-1];
	c_rate_S[I_indiv] = c_rate_S[*numS-1];
	index_S[I_indiv] = index_S[*numS-1];


	// keep track of number of S and E indivs
	--*numS;
	++*numE;
}




// know event is that an indiv moves to infectious class from exposed, need to decide which one and implement
void EtoI(char current_pop_status[], vector<double> &risk_E, vector<double> &lamb_E, vector<double> &c_rate_E, vector<int> &index_E, vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, vector<double> &risk_S, int numS, int *numE, int *numI, double beta_avg, double *infect_rate, mt19937 &rng2)
{
	int I_indiv;	// index of indiv that will become infectious, index in terms of indivs currently E/E vector
	int i;	// param for iteration

	// all indivs have same rate delt so just have to pick number from 0 to numE-1
	uniform_int_distribution<int> pick_EtoI(0, *numE-1);
	I_indiv = pick_EtoI(rng2);	// pick index of indiv who becomes infectious

	// move each param from E to first available row in I - there are numI I indivs
	// if numI < size of vector for I population, can put new indiv at index numI
	if(*numI < risk_I.size())
	{
		risk_I[*numI] = risk_E[I_indiv];
		lamb_I[*numI] = lamb_E[I_indiv];
		c_rate_I[*numI] = c_rate_E[I_indiv];
		index_I[*numI] = index_E[I_indiv];
	}
	else	// if numI == size of vector for I population, need to append new indiv to the end of vector
	{
		risk_I.push_back(risk_E[I_indiv]);
		lamb_I.push_back(lamb_E[I_indiv]);
		c_rate_I.push_back(c_rate_E[I_indiv]);
		index_I.push_back(index_E[I_indiv]);
	}

	// update infect_rate: E indiv is added to I so need to multiply their FoI with all current S indivs and add
	for(i = 0; i < numS; i++)
		*infect_rate += beta_avg * c_rate(i, index_E[I_indiv]) * risk_S[i] * lamb_E[I_indiv];


	// replace data for each param in E with data from last row of numE E indivs since indiv is no longer E
	risk_E[I_indiv] = risk_E[*numE-1];
	lamb_E[I_indiv] = lamb_E[*numE-1];
	c_rate_E[I_indiv] = c_rate_E[*numE-1];
	index_E[I_indiv] = index_E[*numE-1];

	// keep track of number of E and I indivs
	--*numE;
	++*numI;
}




// know event is that an indiv recovers, need to decide which one and implement
void recover(char current_pop_status[], vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, vector<double> &risk_S, int numS, int *numI, int *numR, double beta_avg, double *infect_rate, mt19937 &rng2)
{
	int R_indiv;	// index of indiv that will recover, index in terms of indivs currently I/I vector
	int i;	// param for iteration

	// all indivs have same recovery rate gam so just have to pick number from 0 to numI-1
	uniform_int_distribution<int> pick_recover(0, *numI-1);
	R_indiv = pick_recover(rng2);	// pick index of indiv who recovers

	// update infect_rate: I indiv removed so need to multiply their FoI with all current S indivs and subtract
	for(i = 0; i < numS; i++)
		*infect_rate -= beta_avg * c_rate(i, index_I[R_indiv]) * risk_S[i] * lamb_I[R_indiv];


	// replace data for each param in I with data from last row of numI I indivs since indiv is no longer I
	risk_I[R_indiv] = risk_I[*numI-1];
	lamb_I[R_indiv] = lamb_I[*numI-1];
	c_rate_I[R_indiv] = c_rate_I[*numI-1];
	index_I[R_indiv] = index_I[*numI-1];

	// keep track of number of I indivs
	--*numI;
	++*numR;
}




// know event is that an indiv leaves pop, need to decide which one and implement
void death(vector<double> &risk_S, vector<double> &lamb_S, vector<double> &c_rate_S, vector<int> &index_S, vector<double> &risk_E, vector<double> &lamb_E, vector<double> &c_rate_E, vector<int> &index_E, vector<double> &risk_I, vector<double> &lamb_I, vector<double> &c_rate_I, vector<int> &index_I, int *numS, int *numE, int *numI, int *numR, double beta_avg, double *infect_rate, mt19937 &rng2)
{
	int D_indiv;	// index of indiv that will be removed from pop
	int i;	// param for iteration

	// all indivs have same death rate so just have to pick number from 0 to total pop size-1
	uniform_int_distribution<int> pick_dead(0, *numS+*numE+*numI+*numR-1);
	D_indiv = pick_dead(rng2);	// pick index of indiv who leaves pop

	if(D_indiv < *numS)	// removed indiv is type S
	{
		// update infect_rate: S indiv is removed so need to multiply their risk with all current I indivs and subtract
		for(i = 0; i < *numI; i++)
			*infect_rate -= beta_avg * c_rate(index_S[D_indiv], i) * risk_S[D_indiv] * lamb_I[i];

		// replace data for each param in S with data from last row of numS S indivs
		risk_S[D_indiv] = risk_S[*numS-1];
		lamb_S[D_indiv] = lamb_S[*numS-1];
		c_rate_S[D_indiv] = c_rate_S[*numS-1];
		index_S[D_indiv] = index_S[*numS-1];

		// keep track of number of S indivs
		--*numS;
	}
	else
	{
		if(D_indiv < *numS + *numE)	// removed indiv is type E
		{
			D_indiv -= *numS;	// subtract off number of S indivs to get index in terms of vectors for E indivs

			// replace data for each param in E with data from last row of numE E indivs
			risk_E[D_indiv] = risk_E[*numE-1];
			lamb_E[D_indiv] = lamb_E[*numE-1];
			c_rate_E[D_indiv] = c_rate_E[*numE-1];
			index_E[D_indiv] = index_E[*numE-1];

			// keep track of number of E indivs
			--*numE;
		}
		else
		{
			if(D_indiv < *numS + *numE + *numI)	// removed indiv is type I
			{
				D_indiv -= *numS + *numE;	// subtract off number of S and E indivs to get index in terms of vectors for I indivs

				// update infect_rate: I indiv removed so need to multiply their FoI with all current S indivs and subtract
				for(i = 0; i < *numS; i++)
					*infect_rate -= beta_avg * c_rate(i, index_I[D_indiv]) * risk_S[i] * lamb_I[D_indiv];

				// replace data for each param in I with data from last row of numI I indivs since indiv is no longer I
				risk_I[D_indiv] = risk_I[*numI-1];
				lamb_I[D_indiv] = lamb_I[*numI-1];
				c_rate_I[D_indiv] = c_rate_I[*numI-1];
				index_I[D_indiv] = index_I[*numI-1];

				// keep track of number of I indivs
				--*numI;
			}
			else	// removed indiv is type R
			{
				// keep track of number of R indivs
				--*numR;
			}
		}
	}

}




double c_rate(int S_index, int I_index)
{
	return(1.0);
}




// function to execute R script in a separate file and capture the output of risk and FoI for a new indiv
string exec_Rscript(const char *runR)
{
	char buffer[128];	// buffer to store output in temporarily
	string new_indiv = "";	// output new indiv of R script (risk and FoI)
	FILE *fp = popen(runR, "r");	// run R script

	while(fgets(buffer, sizeof(buffer), fp) != nullptr)	// read output from the R script by line, store in buffer, and append to new_indiv string
		new_indiv += buffer;

	fclose(fp);	// close file stream

	return(new_indiv);
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
		if(strcmp(argv[a], "-inb") == 0 && a + 1 < argc)	// set input file name for births
			inputb = argv[++a];
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
		if(strcmp(argv[a], "-e") == 0 && a + 1 < argc)	// set rate at which individuals move from exposed to infectious
			delt = atof(argv[++a]);
		if(strcmp(argv[a], "-b") == 0 && a + 1 < argc)	// set "birth" rate'
			b = atof(argv[++a]);
		if(strcmp(argv[a], "-d") == 0 && a + 1 < argc)	// set "death" rate'
			d = atof(argv[++a]);
		if(strcmp(argv[a], "-cs") == 0 && a + 1 < argc)	// set cv for HiS
			cv_s = atof(argv[++a]);
		if(strcmp(argv[a], "-ct") == 0 && a + 1 < argc)	// set cv for HiT
			cv_t = atof(argv[++a]);
		if(strcmp(argv[a], "-s") == 0 && a + 1 < argc)	// set seed for simulation
			seed = atoi(argv[++a]);
	}

	return(false);
}

