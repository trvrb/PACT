/* series.cpp
Copyright 2009 Trevor Bedford <bedfordt@umich.edu>
Member function definitions for Series class
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
using std::cout;
using std::endl;

#include <stdexcept>
using std::runtime_error;
using std::out_of_range;

#include <set>
using std::multiset;

#include <cmath>

#include "series.h"

Series::Series() {
	
}

/* inserts value into growing set */
void Series::insert(double n) {
	values.insert(n);
}

/* clears all stored values */
void Series::clear() {
	values.clear();
}

/* returns value at position n */
double Series::at(int n) {

	if (n > values.size() - 1 || values.size() == 0) {
		throw out_of_range("Series::at");
	}

	multiset<double>::iterator is = values.begin();
	for (int i = 0; i < n; i++) {
		is++;
	}
	return *is;

}

/* returns arithmetic mean of stored values */
double Series::mean() {

	if (values.size() == 0) {
		throw out_of_range("Series::mean");
	}

	double mean = 0.0;
	for (multiset<double>::iterator is = values.begin(); is != values.end(); is++ ) {
		mean += *is;
	}
	mean /= (double) values.size();
	return mean;
		
}

/* returns quantile rank (p) of stored values */
/* using the empirical distribution function with averaging */
double Series::quantile(double p) {

	if (values.size() == 0) {
		throw out_of_range("Series::quantile");
	}

	double q;						// quantile value
	
	int n = values.size();			// sample size
	int j = floor(n*p);				// integer part
	double g = n*p - floor(n*p);	// fractional part
	
	if (g > 0 || n == 1)			// samples are ordered 0 to N-1
		q = at(j);				
	else
		q = 0.5 * ( at(j-1) + at(j));

	return q;

}