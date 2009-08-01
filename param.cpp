/* param.cpp
Member function definitions for Parameter class
Parameter values are all global
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

#include "param.h"

Parameters::Parameters() {
	
	// default parameter values
	push_times_back = false;
	push_times_back_start = INF;
	push_times_back_stop = INF;
	
	prune_to_trunk = false;
	
	prune_to_label = false;
	prune_to_label_label = 0;
	
	trim_ends = false;
	trim_ends_start = INF;
	trim_ends_stop = INF;
	
	section_tree = false;
	section_tree_start = INF;
	section_tree_window = INF;
	section_tree_step = INF;
	
	time_slice = false;
	time_slice_time = INF;
	
	print_hp_tree = false;
	
	summary_tmrca = false;		
	summary_length = false;			
	summary_proportions = false;	
	summary_coal_rates = false;		
	summary_mig_rates = false;		
	summary_diversity = false;		
	summary_fst = false;				
	summary_tajima_d = false;	
	
	skyline_start = INF;
	skyline_stop = INF;
	skyline_step = INF;
	
	skyline_tmrca = false;	
	skyline_proportions = false;
	skyline_coal_rates = false;
	skyline_mig_rates = false;
	skyline_diversity = false;
	skyline_fst = false;
	skyline_tajima_d = false;	
	
	cout << "Checking parameters" << endl;
	
	// check to see file exists
	string paramString;
	ifstream paramFile ("in.param");
	if (paramFile.is_open()) {
		while (! paramFile.eof() ) {
			getline (paramFile,paramString);
			importLine(paramString);
		}
	}
	else { 
		cout << "Unable to find in.param" << endl; 
		cout << "Running with default parameters" << endl;
		
		print_hp_tree = true;		
		summary_coal_rates = true;		
		summary_mig_rates = true;
		
	}
	
}

/* Reads a string and attempts to extract parameters from it */
void Parameters::importLine(string line) {

	// READING LINE STRING
	string pstring = "";				// fill with a-z or _
	string vstring = "";				// fill with 0-9 or . or -
	vector<double> values;				// convert vstring to double and push here

	for (string::iterator is = line.begin(); is != line.end(); is++) {
	
		if (*is == '#')		// ignore rest of line after comment
			break; 
		else if (*is >= 'a' && *is <= 'z')
    		pstring += *is;
		else if ( (*is >= '0' && *is <= '9') || *is == '.' || *is == '-')
    		vstring += *is;
    	else if (vstring.size() > 0) {
    		values.push_back( atof(vstring.c_str()) );
    		vstring = "";
    	}
		
	}
	
	if (vstring.size() > 0) {
		values.push_back( atof(vstring.c_str()) );
	}	
		
		
	// SETTING PARAMETERS	
	if (pstring == "pushtimesback") { 
		if (values.size() == 1) {
			push_times_back = true; 
			push_times_back_stop = values[0];			
		}
		if (values.size() == 2) {
			push_times_back = true; 
			push_times_back_start = values[0];	
			push_times_back_start = values[1];		
		}		
	}
	
	if (pstring == "prunetotrunk") { 
		prune_to_trunk = true; 
	}
	
	if (pstring == "prunetolabel") { 
		if (values.size() == 1) {
			prune_to_label = true; 
			prune_to_label_label = (int) values[0];			
		}
	}	
	
	if (pstring == "trimends") { 
		if (values.size() == 2) {
			trim_ends = true; 
			trim_ends_start = values[0];
			trim_ends_stop = values[1];	
		}
	}	
	
	if (pstring == "sectiontree") { 
		if (values.size() == 3) {
			section_tree = true; 
			section_tree_start = values[0];
			section_tree_window = values[1];
			section_tree_step = values[2];
		}
	}	

	if (pstring == "timeslice") { 
		if (values.size() == 1) {
			time_slice = true; 
			time_slice_time = values[0];	
		}
	}		
	
	if (pstring == "printhptree") { 
		print_hp_tree = true; 
	}	

	if (pstring == "summarytmrca") { summary_tmrca = true; }
	if (pstring == "summarylength") { summary_length = true; }
	if (pstring == "summaryproportions") { summary_proportions = true; }
	if (pstring == "summarycoalrates") { summary_coal_rates = true; }
	if (pstring == "summarymigrates") { summary_mig_rates = true; }
	if (pstring == "summarydiversity") { summary_diversity = true; }
	if (pstring == "summaryfst") { summary_fst = true; }
	if (pstring == "summarytajimad") { summary_tajima_d = true; }	
	
	if (pstring == "skylinesettings") { 
		if (values.size() == 1) {
			skyline_step = values[0];			
		}
		if (values.size() == 3) {
			skyline_start = values[0]; 
			skyline_stop = values[1];	
			skyline_step = values[2];		
		}		
	}

	if (pstring == "skylinetmrca") { skyline_tmrca = true; }
	if (pstring == "skylineproportions") { skyline_proportions = true; }
	if (pstring == "skylinecoalrates") { skyline_coal_rates = true; }
	if (pstring == "skylinemigrates") { skyline_mig_rates = true; }
	if (pstring == "skylinediversity") { skyline_diversity = true; }
	if (pstring == "skylinefst") { skyline_fst = true; }
	if (pstring == "skylinetajimad") { skyline_tajima_d = true; }		
	
}

/* prints parameters */
void Parameters::print() {

	

}