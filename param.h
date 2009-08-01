/* param.h
This contains all the global variables used by pact, can be modified later if necessary
*/

#ifndef PARAM_H
#define PARAM_H

#include <cmath>

#define INF pow(double(10),double(100)) // 10^100 (~infinity)

// defaults
// defaulting times to infinity
bool push_times_back = false;
double push_times_back_start = INF;
double push_times_back_stop = INF;

bool prune_to_trunk = false;

bool prune_to_label = false;
int prune_to_label_label = 0;

bool trim_ends = false;
double trim_ends_start = INF;
double trim_ends_stop = INF;

bool section = false;
double section_start = INF;
double section_window = INF;
double section_step = INF;

bool time_slice = false;
double time_slice_time = INF;

bool print_hp_tree = true;

bool summary_tmrca = false;		
bool summary_length = false;			
bool summary_proportions = false;	
bool summary_coal_rates = true;		
bool summary_mig_rates = true;		
bool summary_diversity = false;		
bool summary_fst = false;				
bool summary_tajima_d = false;	

double skyline_start = INF;
double skyline_stop = INF;
double skyline_step = INF;

bool skyline_tmrca = false;	
bool skyline_proportions = false;
bool skyline_coal_rates = false;
bool skyline_mig_rates = false;
bool skyline_diversity = false;
bool skyline_fst = false;
bool skyline_tajima_d = false;


class Parameters {

public:
	Parameters();						// constructor, imports parameters from in.param if available

	void print();						// prints parameter listing
	
private:
	void importLine(string);			// reads a string and attempts to extract parameters from it
	
	
};

#endif