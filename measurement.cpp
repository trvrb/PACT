/* measurement.cpp
Member function definitions for Measurement class
*/

#include "measurement.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <sstream>

Measurement::Measurement() {
	
	n = 0;
	
}

void Measurement::print() {

	if (n > 0) {

		for (int i = 0; i < n - 1; i++) {
			cout << values[i] << " ";
		}
		cout << values[n-1] << endl;
	
	}
	
	else {
		cout << endl;
	}

}

vector<double> Measurement::get() { 
	return values; 
}

int Measurement::size() {
	return n;
}

void Measurement::clear() {
	n = 0;
	values.clear();
}

void Measurement::increment(vector<double> in) {

	if (n > 0) {			// already have data going
	
		for (int i = 0; i < n; i++) {
			values[i] += in[i];
		}
	
	}
	else {					// new data
	
		for (int i = 0; i < in.size(); i++) {
			values.push_back(in[i]);
		}
		n = values.size();
	}
	

}

void Measurement::divideBy(double in) {

	for (int i = 0; i < n; i++) {
		values[i] /= in;
	}
	
}

void Measurement::divideBy(vector<double> in) {

	for (int i = 0; i < n; i++) {
		values[i] /= in[i];
	}
	
}

double Measurement::total() {

	double v;
	for (int i = 0; i < n; i++) {
		v += values[i];
	}
	
	return v;
	
}

double Measurement::mean() {

	double v;
	for (int i = 0; i < n; i++) {
		v += values[i];
	}
	v /= (double) n;
	
	return v;
	
}
