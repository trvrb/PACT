/* param.h
This contains all the global variables used by pact, can be modified later if necessary
*/

#ifndef PARAM_H
#define PARAM_H

#include <cmath>

#define INF pow(double(10),double(100)) // 10^100 (~infinity)

class Parameters {

public:
	Parameters();						// constructor, imports parameters from in.param if available

	void print();						// prints parameter listing

	// PARAMETERS
	bool push_times_back;
	double push_times_back_start;
	double push_times_back_stop;
	
	bool prune_to_trunk;
	
	bool prune_to_label;
	int prune_to_label_label;
	
	bool trim_ends;
	double trim_ends_start;
	double trim_ends_stop;
	
	bool section_tree;
	double section_tree_start;
	double section_tree_window;
	double section_tree_step;
	
	bool time_slice;
	double time_slice_time;
	
	bool print_hp_tree;
	
	bool summary_tmrca;		
	bool summary_length;			
	bool summary_proportions;	
	bool summary_coal_rates;		
	bool summary_mig_rates;		
	bool summary_diversity;		
	bool summary_fst;				
	bool summary_tajima_d;	
	
	double skyline_start;
	double skyline_stop;
	double skyline_step;
	
	bool skyline_tmrca;	
	bool skyline_proportions;
	bool skyline_coal_rates;
	bool skyline_mig_rates;
	bool skyline_diversity;
	bool skyline_fst;
	bool skyline_tajima_d;
	
private:
	void importLine(string);			// reads a string and attempts to extract parameters from it
	
	
};

#endif