/* stat.cpp
Member function definitions for Statistic class
*/

#include "stat.h"

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

Statistic::Statistic() {
	
	n = 0;
	
}

void Statistic::print() {

	for (int i = 0; i < n - 1; i++) {
		cout << values[i] << " ";
	}
	cout << values[n-1] << endl;

}

vector<double> Statistic::get() {

	return values;

}

void Statistic::increment(vector<double> in) {

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

void Statistic::divideBy(double in) {

	for (int i = 0; i < n; i++) {
		values[i] /= in;
	}
	
}

void Statistic::divideBy(vector<double> in) {

	for (int i = 0; i < n; i++) {
		values[i] /= in[i];
	}
	
}

double Statistic::total() {

	double v;
	for (int i = 0; i < n; i++) {
		v += values[i];
	}
	
	return v;
	
}

double Statistic::mean() {

	double v;
	for (int i = 0; i < n; i++) {
		v += values[i];
	}
	v /= (double) n;
	
	return v;
	
}
