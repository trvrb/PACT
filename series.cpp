/* series.cpp
Copyright 2009-2010 Trevor Bedford <bedfordt@umich.edu>
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
using std::isnan;
using std::isinf;
using std::sqrt;

#include "series.h"

Series::Series() {
	
}

/* inserts value into growing set */
void Series::insert(double n) {
	if (!isnan(n) && !isinf(n)) {
		values.insert(n);
	}
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

	double mean = 0.0;

	if (values.size() != 0) {
		for (multiset<double>::iterator is = values.begin(); is != values.end(); is++ ) {
			mean += *is;
		}
		mean /= (double) values.size();
	}
	else {
		mean /= mean;
	}
	
	return mean;
		
}

/* returns quantile rank (p) of stored values */
/* using the empirical distribution function with averaging */
double Series::quantile(double p) {

	double q = 0.0;						// quantile value

	if (values.size() != 0) {
		int n = values.size();			// sample size
		int j = floor(n*p);				// integer part
		double g = n*p - floor(n*p);	// fractional part
		
		if (g > 0 || n == 1)			// samples are ordered 0 to N-1
			q = at(j);				
		else
			q = 0.5 * ( at(j-1) + at(j));
	}
	else {
		q /= q; 
	}

	return q;

}

/* returns standard deviation of stored values */
double Series::sd() {

	double sd = 0.0;
	double seriesmean = mean();

	if (values.size() != 0) {
		for (multiset<double>::iterator is = values.begin(); is != values.end(); is++ ) {
			sd += (*is - seriesmean) * (*is - seriesmean);
		}
		sd /= (double) (values.size() - 1);
		sd = sqrt(sd);
	}
	else {
		sd /= sd;
	}
	
	return sd;
		
}

/* returns x standard deviations up or down from the mean */
double Series::sdrange(double x) {
	
	return mean() + x*sd();
		
}