/* measurement.h
Measurement class definition
This object represents a tree measurement, like TMRCA or migration rates
Mainly exists to ease working with seasonal averages
*/

#ifndef MEASURE_H
#define MEASURE_H

#include <map>
using std::map;
#include <set>
using std::set;
#include <vector>
using std::vector;

class Measurement {

public:
	Measurement();						// constructor
	
	void print();						// print space deliminated
	vector<double> get();				// returns vector of values	
	int size();							// returns size of vector
	void clear();						// clears data
	
	void increment(vector<double>);		// increment
	void divideBy(double);				// divide all values by the same amount
	void divideBy(vector<double>);		// divide each value by a particular amount
	
	double total();						// returns the sum of the values vector
	double mean();						// returns the mean of the values vector
																	
private:
	vector<double> values;				// vector of values
	int n;								// length of values

};

#endif