/* param.h
Copyright 2009-2012 Trevor Bedford <t.bedford@ed.ac.uk>
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
	bool printtree();					// are any of the tree printing parameters true?
	bool summary();						// are any summary parameters true?
	bool tips();						// are any of the tip statistics true?
	bool skyline();						// are any skyline parameters true?
	bool pairs();						// are any of the pair statistics true?

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
	
	bool prune_to_tips;
	vector<string> prune_to_tips_values;	// list of tips
	
	bool remove_tips;
	vector<string> remove_tips_values;		// list of tips	
	
	bool prune_to_time;
	vector<double> prune_to_time_values;	// start, stop	
	
	bool collapse_labels;
	
	bool pad_migration_events;
	
	bool trim_ends;
	vector<double> trim_ends_values;		// start, stop
	
	bool section_tree;
	vector<double> section_tree_values;		// start, window, step
	
	bool time_slice;
	vector<double> time_slice_values;		// time
	
	bool rotate;
	vector<double> rotate_values;			// degrees	
	
	bool add_tail;
	vector<double> add_tail_values;			// time		
	
	bool print_tree;
	bool print_circular_tree;	
	bool print_all_trees;	
	
	bool summary_tmrca;		
	bool summary_length;			
	bool summary_root_proportions;	
	bool summary_proportions;	
	bool summary_coal_rates;		
	bool summary_mig_rates;	
	bool summary_sub_rates;		
	bool summary_diversity;		
	bool summary_fst;				
	bool summary_tajima_d;	
	bool summary_diffusion_coefficient;
	bool summary_drift_rate;
	bool summary_persistence;
	
	bool tips_time_to_trunk;
	bool x_loc_history;
	vector<double> x_loc_history_values;	// start, stop, step
	bool y_loc_history;
	vector<double> y_loc_history_values;	// start, stop, step	
	bool coord_history;
	vector<double> coord_history_values;	// start, stop, step	
		
	vector<double> skyline_values;			// start, stop, step
	
	bool skyline_tmrca;	
	bool skyline_length;
	bool skyline_proportions;
	bool skyline_coal_rates;
	bool skyline_mig_rates;
	bool skyline_pro_history_from_tips;
	bool skyline_diversity;
	bool skyline_fst;
	bool skyline_tajima_d;
	bool skyline_timetofix;
	bool skyline_xmean;
	bool skyline_ymean;
	bool skyline_xdrift;
	bool skyline_ratemean;
	bool skyline_xtrunkdiff;
	bool skyline_locsample;	
	bool skyline_locgrid;		
	bool skyline_drift_rate_from_tips;
	
	bool ordering;
	vector<string> ordering_values;
	
	bool pairs_diversity;
	vector<double> pairs_diversity_values;	// time diff	
	
private:
	void importLine(string);			// reads a string and attempts to extract parameters from it
	
	
};

#endif