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