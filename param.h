/* param.h
Copyright 2009-2010 Trevor Bedford <bedfordt@umich.edu>
This contains all the global variables used by pact, can be modified later if necessary
*/

/*
This file is part of PACT.

PACT is free software: you can redistribute it and/or modify it under the terms of the GNU General 
Public License as published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PACT is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the 
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
Public License for more details.

You should have received a copy of the GNU General Public License along with PACT.  If not, see 
<http://www.gnu.org/licenses/>.
*/

#ifndef PARAM_H
#define PARAM_H

#include <cmath>

#include <vector>
using std::vector;

class Parameters {

public:
	Parameters();						// constructor, imports parameters from in.param if available

	void print();						// prints parameter listing

	bool general();						// are any of the general parameters true?
	bool manip();						// are any tree manipulation parameters true?
	bool summary();						// are any summary parameters true?
	bool tips();
	bool skyline();						// are any skyline parameters true?

	// PARAMETERS
	
	bool burnin;
	vector<double> burnin_values;			// count
	
	bool push_times_back;
	vector<double> push_times_back_values;	// start, stop
	
	bool reduce_tips;
	vector<double> reduce_tips_values;		// count
	
	bool renew_trunk;
	vector<double> renew_trunk_values;		// time
	
	bool prune_to_trunk;
	
	bool prune_to_label;
	vector<string> prune_to_label_values;	// label
	
	bool prune_to_time;
	vector<double> prune_to_time_values;	// start, stop	
	
	bool collapse_labels;
	
	bool trim_ends;
	vector<double> trim_ends_values;		// start, stop
	
	bool section_tree;
	vector<double> section_tree_values;		// start, window, step
	
	bool time_slice;
	vector<double> time_slice_values;		// time
	
	bool print_rule_tree;
	
	bool summary_tmrca;		
	bool summary_length;			
	bool summary_proportions;	
	bool summary_coal_rates;		
	bool summary_mig_rates;		
	bool summary_diversity;		
	bool summary_fst;				
	bool summary_tajima_d;	
	
	bool tips_time_to_trunk;
		
	vector<double> skyline_values;			// start, stop, step
	
	bool skyline_tmrca;	
	bool skyline_length;
	bool skyline_proportions;
	bool skyline_coal_rates;
	bool skyline_mig_rates;
	bool skyline_diversity;
	bool skyline_fst;
	bool skyline_tajima_d;
	bool skyline_timetofix;	
	
private:
	void importLine(string);			// reads a string and attempts to extract parameters from it
	
	
};

#endif