/* param.cpp
Copyright 2009-2012 Trevor Bedford <t.bedford@ed.ac.uk>
Member function definitions for Parameter class
Parameter values are all global
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

#include <iostream>
#include <fstream>
using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::ios;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <cstdlib>
using std::atof;

#include <stdexcept>
using std::runtime_error;
using std::out_of_range;

#include "param.h"

Parameters::Parameters() {
	
	// default parameter values
	// leaving value vectors empty purposely
	burnin = false;
	push_times_back = false;
	reduce_tips = false;
	renew_trunk = false;
	prune_to_trunk = false;
	prune_to_time = false;
	prune_to_label = false;
	collapse_labels = false;
	trim_ends = false;
	section_tree = false;	
	time_slice = false;
	rotate = false;
	add_tail = false;
	
	print_rule_tree = false;
	print_all_trees = false;	
	
	summary_tmrca = false;		
	summary_length = false;			
	summary_proportions = false;	
	summary_coal_rates = false;		
	summary_mig_rates = false;		
	summary_diversity = false;		
	summary_fst = false;				
	summary_tajima_d = false;	
	summary_diffusion_coefficient = false;	
	
	tips_time_to_trunk = false;
	x_loc_history = false;	
	y_loc_history = false;		
	coord_history = false;		
	
	skyline_tmrca = false;	
	skyline_length = false;
	skyline_proportions = false;
	skyline_coal_rates = false;
	skyline_mig_rates = false;
	skyline_diversity = false;
	skyline_fst = false;
	skyline_tajima_d = false;
	skyline_timetofix = false;
	skyline_xmean = false;
	skyline_ymean = false;
	skyline_xdrift = false;	
	skyline_ratemean = false;
	skyline_xtrunkdiff = false;	
	skyline_locsample = false;
	skyline_locgrid = false;	
	skyline_drift_rate_from_tips = false;
	
	ordering = false;
	
	pairs_diversity = false;
	
	// check to see file exists
	string paramString;
	ifstream paramFile ("in.param");
	if (paramFile.is_open()) {
	
		cout << "Reading parameters from in.param" << endl;
		cout << endl;
	
		while (! paramFile.eof() ) {
			getline (paramFile,paramString);
			importLine(paramString);
		}
	}
	else { 
		throw runtime_error("parameter file in.param not found");
//		cout << "Unable to find in.param" << endl; 
//		cout << "Running with default parameters" << endl;
//		cout << endl;
		
//		print_rule_tree = true;		
//		summary_coal_rates = true;		
//		summary_mig_rates = true;
		
	}
	
}

/* Reads a string and attempts to extract parameters from it */
void Parameters::importLine(string line) {
	
	// READING LINE STRING
	string pstring = "";				// fill with a-z or _
	string vstring = "";				// fill with 0-9 or . or - or A-Z
	vector<double> values;				// convert vstring to double and push here
	vector<string> svalues;

	for (string::iterator is = line.begin(); is != line.end(); is++) {
	
		if (*is == '#')		// ignore rest of line after comment
			break; 
		else if (*is >= 'a' && *is <= 'z' && vstring.size() == 0)
    		pstring += *is;
		else if ( (*is >= '0' && *is <= '9') || (*is >= 'a' && *is <= 'z') || (*is >= 'A' && *is <= 'Z') || *is == '.' || *is == '-')
    		vstring += *is;
    	else if (vstring.size() > 0) {
    		values.push_back( atof(vstring.c_str()) );
    		svalues.push_back(vstring);
    		vstring = "";
    	}
		
	}
	
	if (vstring.size() > 0) {
		values.push_back( atof(vstring.c_str()) );
	}	
	if (vstring.size() > 0) {
		svalues.push_back(vstring);
	}		
		
		
	// SETTING PARAMETERS	

	if (pstring == "burnin") { 
		if (values.size() == 1) {
			burnin = true; 
			burnin_values = values;			
		}
	}		
	
	if (pstring == "pushtimesback") { 
		if (values.size() == 1 || values.size() == 2) {
			push_times_back = true; 
			push_times_back_values = values;			
		}	
	}
	
	if (pstring == "reducetips") { 
		if (values.size() == 1) {
			reduce_tips = true; 
			reduce_tips_values = values;			
		}
	}	

	if (pstring == "renewtrunk") { 
		if (values.size() == 1) {
			renew_trunk = true;
			renew_trunk_values = values;			
		} 
	}
	
	if (pstring == "prunetotrunk") { 
		prune_to_trunk = true;
	}
	
	if (pstring == "prunetolabel") { 
		if (values.size() == 1) {
			prune_to_label = true; 			
			prune_to_label_values = svalues;			
		}
	}	

	if (pstring == "prunetotime") { 
		if (values.size() == 2) {
			prune_to_time = true; 
			prune_to_time_values = values;
		}
	}		
	
	if (pstring == "collapselabels") { 
		collapse_labels = true;
	}	
	
	if (pstring == "trimends") { 
		if (values.size() == 2) {
			trim_ends = true; 
			trim_ends_values = values;
		}
	}	
	
	if (pstring == "sectiontree") { 
		if (values.size() == 3) {
			section_tree = true; 
			section_tree_values = values;
		}
	}	

	if (pstring == "timeslice") { 
		if (values.size() == 1) {
			time_slice = true; 
			time_slice_values = values;	
		}
	}	
	
	if (pstring == "rotate") { 
		if (values.size() == 1) {
			rotate = true; 
			rotate_values = values;	
		}
	}	
	
	if (pstring == "addtail") { 
		if (values.size() == 1) {
			add_tail = true; 
			add_tail_values = values;	
		}
	}		
	
	if (pstring == "ordering") { 
		ordering = true; 				
		ordering_values = svalues;	
	}		
	
	if (pstring == "printruletree") { 
		print_rule_tree = true; 
	}	
	
	if (pstring == "printalltrees") { 
		print_all_trees = true; 
	}		

	if (pstring == "summarytmrca") { summary_tmrca = true; }
	if (pstring == "summarylength") { summary_length = true; }
	if (pstring == "summaryproportions") { summary_proportions = true; }
	if (pstring == "summarycoalrates") { summary_coal_rates = true; }
	if (pstring == "summarymigrates") { summary_mig_rates = true; }
	if (pstring == "summarysubrates") { summary_sub_rates = true; }	
	if (pstring == "summarydiversity") { summary_diversity = true; }
	if (pstring == "summaryfst") { summary_fst = true; }
	if (pstring == "summarytajimad") { summary_tajima_d = true; }
	if (pstring == "summarydiffusioncoefficient") { summary_diffusion_coefficient = true; }
	if (pstring == "summarydriftrate") { summary_drift_rate = true; }	
	
	if (pstring == "tipstimetotrunk") { tips_time_to_trunk = true; }	
	
	if (pstring == "tipsxlochistory") { 
		if (values.size() == 3) {
			x_loc_history = true; 
			x_loc_history_values = values;	
		}
		
	}	

	if (pstring == "tipsylochistory") { 
		if (values.size() == 3) {
			y_loc_history = true; 
			y_loc_history_values = values;	
		}
		
	}		

	if (pstring == "coordhistory") { 
		if (values.size() == 3) {
			x_loc_history = true; 
			x_loc_history_values = values;	
		}
		
	}		
	
	if (pstring == "skylinesettings") { 
		if (values.size() == 3) { 
			skyline_values = values;			
		}			
	}

	if (pstring == "skylinetmrca") { skyline_tmrca = true; }
	if (pstring == "skylinelength") { skyline_length = true; }
	if (pstring == "skylineproportions") { skyline_proportions = true; }
	if (pstring == "skylinecoalrates") { skyline_coal_rates = true; }
	if (pstring == "skylinemigrates") { skyline_mig_rates = true; }
	if (pstring == "skylinediversity") { skyline_diversity = true; }
	if (pstring == "skylinefst") { skyline_fst = true; }
	if (pstring == "skylinetajimad") { skyline_tajima_d = true; }
	if (pstring == "skylinetimetofix") { skyline_timetofix = true; }
	if (pstring == "skylinexmean") { skyline_xmean = true; }	
	if (pstring == "skylineymean") { skyline_ymean = true; }
	if (pstring == "skylinexdrift") { skyline_xdrift = true; }		
	if (pstring == "skylineratemean") { skyline_ratemean = true; }
	if (pstring == "skylinextrunkdiff") { skyline_xtrunkdiff = true; }		
	if (pstring == "skylinelocsample") { skyline_locsample = true; }	
	if (pstring == "skylinelocgrid") { skyline_locgrid = true; }	
	if (pstring == "skylinedriftratefromtips") { skyline_drift_rate_from_tips = true; }		
	
	if (pstring == "pairsdiversity") { 
		if (values.size() == 1) {
			pairs_diversity = true; 
			pairs_diversity_values = values;	
		}
	}	
		
}

/* prints parameters */
void Parameters::print() {

	// GENERAL
	if ( general() ) {
		
		cout << "General:" << endl;
		
		if (burnin) {
			cout << "burnin " << burnin_values[0] << endl;
		}	
	
		cout << endl;
	
	}

	// TREE MANIPULATION
	if ( manip() ) {
	
		cout << "Tree manipulation:" << endl;
		
		if (push_times_back) {
			cout << "push times back ";
			for (int i = 0; i < push_times_back_values.size(); i++) {
				cout << push_times_back_values[i] << " ";
			}
			cout << endl;
		}
		
		if (reduce_tips) {
			cout << "reduce tips " << reduce_tips_values[0] << endl;
		}		

		if (renew_trunk) {
			cout << "renew trunk " << renew_trunk_values[0] << endl;
		}	
			
		if (trim_ends) {
			cout << "trim ends " << trim_ends_values[0] << " " << trim_ends_values[1] << endl;		
		}
	
		if (section_tree) {
			cout << "section tree " << section_tree_values[0] << " " << section_tree_values[1] << " " << section_tree_values[2] << endl;		
		}	
	
		if (time_slice) {
			cout << "time slice " << time_slice_values[0] << endl;
		}		
		if (prune_to_label) {
			cout << "prune to label " << prune_to_label_values[0] << endl;
		}
	
		if (prune_to_trunk) {
			cout << "prune to trunk" << endl;
		}	
		
		if (prune_to_time) {
			cout << "prune to time " << prune_to_time_values[0] << " " << prune_to_time_values[1] << endl;
		}		
		
		if (collapse_labels) {
			cout << "collapse labels" << endl;
		}	
		
		if (rotate) {
			cout << "rotate " << rotate_values[0] << endl;
		}	
		
		if (add_tail) {
			cout << "add tail " << add_tail_values[0] << endl;
		}		
		
		// TIP ORDERING
		if (ordering) {
			cout << "Tip ordering:" << endl;
			for (int i = 0; i < ordering_values.size(); i++) {
				cout << ordering_values[i] << " ";
			}
			cout << endl;
		}			
						
		cout << endl;
	
	}
	
	// TREE STRUCTURE
	if ( printtree() ) {
		cout << "Tree structure:" << endl;
		if (print_rule_tree) { cout << "print rule tree" << endl; }
		if (print_all_trees) { cout << "print all trees" << endl; }
		cout << endl;
	}	
	
	// SUMMARY STATISTICS
	if ( summary() ) { 
		cout << "Summary statistics:" << endl;
		if (summary_tmrca) { cout << "tmrca" << endl; }
		if (summary_length) { cout << "length" << endl; }
		if (summary_proportions) { cout << "proportions" << endl; }
		if (summary_coal_rates) { cout << "coal rates" << endl; }
		if (summary_mig_rates) { cout << "mig rates" << endl; }
		if (summary_sub_rates) { cout << "sub rates" << endl; }		
		if (summary_diversity) { cout << "diversity" << endl; }
		if (summary_fst) { cout << "fst" << endl; }	
		if (summary_tajima_d) { cout << "tajima d" << endl; }	
		cout << endl; 
	}

	// TIP STATISTICS
	if ( tips() ) { 
		cout << "Tip statistics:" << endl;
		if (tips_time_to_trunk) { cout << "time to trunk" << endl; }
		if (x_loc_history) { 
			cout << "x loc history " << x_loc_history_values[0] << " " << x_loc_history_values[1] << " " << x_loc_history_values[2] << endl; 
		}	
		if (y_loc_history) { 
			cout << "y loc history " << y_loc_history_values[0] << " " << y_loc_history_values[1] << " " << y_loc_history_values[2] << endl; 
		}			
		cout << endl; 
	}	

	// SKYLINE STATISTICS
	if ( skyline() ) {
		cout << "Skyline statistics:";
		cout << " " << skyline_values[0];
		cout << " " << skyline_values[1];
		cout << " " << skyline_values[2] << endl;
		if (skyline_tmrca) { cout << "tmrca" << endl; }
		if (skyline_length) { cout << "length" << endl; }
		if (skyline_proportions) { cout << "proportions" << endl; }
		if (skyline_coal_rates) { cout << "coal rates" << endl; }
		if (skyline_mig_rates) { cout << "mig rates" << endl; }
		if (skyline_diversity) { cout << "diversity" << endl; }
		if (skyline_fst) { cout << "fst" << endl; }	
		if (skyline_tajima_d) { cout << "tajima d" << endl; }
		if (skyline_timetofix) { cout << "time to fix" << endl; }
		if (skyline_xmean) { cout << "x mean" << endl; }
		if (skyline_ymean) { cout << "y mean" << endl; }
		if (skyline_xdrift) { cout << "x drift" << endl; }		
		if (skyline_ratemean) { cout << "rate mean" << endl; }	
		if (skyline_xtrunkdiff) { cout << "x trunk diff" << endl; }
		if (skyline_locsample) { cout << "loc sample" << endl; }
		if (skyline_locgrid) { cout << "loc grid" << endl; }	
		if (skyline_drift_rate_from_tips) { cout << "drift rate from tips" << endl; }	
		cout << endl;
	}
	
	// PAIR STATISTICS
	if ( pairs() ) {
		cout << "Pair statistics:" << endl;
		if (pairs_diversity) { cout << "pairwise diversity" << endl; }		
	}
	
}

bool Parameters::general() {
	bool check;
	if (burnin)
		check = true;
	else 
		check = false;
	return check;
}

bool Parameters::manip() {
	bool check;
	if (push_times_back || reduce_tips || renew_trunk || prune_to_trunk || prune_to_label || collapse_labels || trim_ends || section_tree || time_slice || rotate || add_tail || ordering)
		check = true;
	else 
		check = false;
	return check;
}

bool Parameters::printtree() {
	bool check;
	if (print_rule_tree || print_all_trees)
		check = true;
	else 
		check = false;
	return check;
}

bool Parameters::summary() {
	bool check;
	if (summary_tmrca || summary_length || summary_proportions || summary_coal_rates || summary_mig_rates || summary_sub_rates || summary_diversity || summary_fst || summary_tajima_d || summary_diffusion_coefficient)
		check = true;
	else 
		check = false;
	return check;
}

bool Parameters::tips() {
	bool check;
	if (tips_time_to_trunk || x_loc_history || y_loc_history || coord_history)
		check = true;
	else 
		check = false;
	return check;
}

bool Parameters::skyline() {
	bool check;
	if ( skyline_values.size() == 3 &&
		(skyline_tmrca || skyline_length || skyline_proportions || skyline_coal_rates || skyline_mig_rates || skyline_diversity || skyline_fst || skyline_tajima_d || skyline_timetofix || skyline_xmean || skyline_ymean || skyline_xdrift || skyline_ratemean || skyline_xtrunkdiff || skyline_locsample || skyline_locgrid || skyline_drift_rate_from_tips) )
		check = true;
	else 
		check = false;
	return check;
}

bool Parameters::pairs() {
	bool check;
	if (pairs_diversity)
		check = true;
	else
		check = false;
	return check;
}