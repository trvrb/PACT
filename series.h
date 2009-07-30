/* series.h
Series class definition
Extends vector
This object holds a series of simple measurements, and associated operations.
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
	double quantile(double);				// returns the quantile rank of the stored values
											
private:
	multiset<double> values;				// measurement values, never ordered
											// multiset insures that these are always kept sorted

};

#endif