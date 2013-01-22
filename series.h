/* series.h
Copyright 2009-2012 Trevor Bedford <t.bedford@ed.ac.uk>
Series class definition
This object holds a series of simple measurements, and associated operations.
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

#ifndef SERIES_H
#define SERIES_H

#include <set>
using std::multiset;

class Series {

public:
	Series();								// constructor
	
	void insert(double);					// inserts a value into the growing set
	void clear();							// clears all stored values
	
	double at(int);							
	double mean();							// returns arithmetic mean of stored values 
	double median();						// returns the median of stored values 	
	double quantile(double);				// returns the quantile rank of the stored values
	double sd();							// returns standard deviation of the stored values
	double sdrange(double);					// returns x standard deviations up or down from the mean
											
private:
	multiset<double> values;				// measurement values, never ordered
											// multiset insures that these are always kept sorted

};

#endif