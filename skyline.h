/* skyline.h
Skyline class definition
This object stores and manipulates skyline data, a time series of parameter estimates measured across
a number of MCMC runs
*/

#ifndef SKYLINE_H
#define SKYLINE_H

#include <map>
using std::map;
#include <set>
using std::set;
#include <vector>
using std::vector;

class Skyline {

public:
	Skyline();												// constructor
															// don't know vector sizes ahead of time
	
	void appendIteration(vector<double>, vector<double>);	// takes a vector of time points and a vector
															// of values and appends them to the skyline object
	void printStatistics();									// goes through Skyline object and prints:
															// timepoint   95%lower   median   95%upper
	void printMeans();										// goes through Skyline object and print means
	
					
private:
	vector< vector<double> > skylineindices;				// holds time points as [iteration] by [step]
	vector< vector<double> > skylinevalues;					// holds values as [iteration] by [step]
	int maxiter;											// iteration with largest number of steps
	int maxsteps;											// largest number of steps over all iterations
	int iter;												// last iteration in skyline object
	
	vector<double> stats(int);								// returns the <95%lower,median,95%upper> based on row number
	double mean(int);										// returns the mean based on row number
	
};

#endif