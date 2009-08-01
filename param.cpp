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
	
	section = false;
	section_start = INF;
	section_window = INF;
	section_step = INF;
	
	time_slice = false;
	time_slice_time = INF;
	
	print_hp_tree = true;
	
	summary_tmrca = false;		
	summary_length = false;			
	summary_proportions = false;	
	summary_coal_rates = true;		
	summary_mig_rates = true;		
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
	}
	
}

/* Reads a string and attempts to extract parameters from it */
void Parameters::importLine(string line) {

	string pstring = "";				// fill with a-z or _
	string vstring = "";				// fill with 0-9 or . or -
	vector<double> values;				// convert vstring to double and push here

	// fill pstring and values vector
	for (string::iterator is = line.begin(); is != line.end(); is++) {
	
		if (*is == '#')		// ignore rest of line after comment
			break; 
		else if (*is >= 'a' && *is <= 'z')
    		pstring += *is;
		else if ( (*is >= '0' && *is <= '9') || *is == '.' || *is == '-')
    		vstring += *is;
    	else {
    		values.push_back( atof(vstring.c_str()) );
    	}
		
	}
	
	
}

/* prints parameters */
void Parameters::print() {

	

}