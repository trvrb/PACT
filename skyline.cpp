/* skyline.cpp
Member function definitions for Skyline class
*/

#include "skyline.h"

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

#define INF pow(double(10),double(100)) // 10^100 (~infinity)

/* Constructor function to initialize private data */
Skyline::Skyline() {
	maxsteps = 0;
	maxiter = 0;
	iter = -1;
}

void Skyline::appendIteration(vector<double> times, vector<double> values) {

	skylineindices.push_back(times);
	skylinevalues.push_back(values);
	iter++;

	if (maxsteps < times.size()) {
		maxsteps = times.size();
		maxiter = iter;
	}
				
}

void Skyline::printStatistics() {

	vector<double> s;
	for (int j=0; j < maxsteps; j++) {
		s = stats(j);
		if (s[0] < INF)
			cout << skylineindices[maxiter][j] << "\t" << s[0] << "\t" << s[1] << "\t" << s[2] << endl;
		else
			cout << skylineindices[maxiter][j] << "\tnan\tnan\tnan" << endl;
	}

}

/* ignore non-numbers */
vector<double> Skyline::stats(int row) {

	double lower = 0.0;
	double median = 0.0;
	double upper = 0.0;
	vector<double> v(3);
	multiset<double> valueset;
	
	/* move through iters and add to a set if a value exists for this row */
	for (int i=0; i <= iter; i++) {
		int size = skylineindices[i].size();
		if (row < size && !isnan(skylinevalues[i][row])) {
			valueset.insert(skylinevalues[i][row]);
		}
	}
			
	int count = valueset.size();
	int pos;
	int get;
	
	if (count > 0) {
	
		/* 95% lower */
		pos = floor(count * 0.025);
		get = 0;
		for (set<double>::iterator it=valueset.begin(); it!=valueset.end(); it++) {
			if (get == pos) {
				lower = *it;	
			}
			get++;
		}
		
		/* median */
		if (count % 2 == 1) {  		// odd, take the center element
			pos = floor(count/2);
			get = 0;
			for (set<double>::iterator it=valueset.begin(); it!=valueset.end(); it++) {
				if (get == pos) {
					median = *it;	
				}
				get++;
			}
		}
		else {						// even, take the mean of the two center elements
			int posA = count/2;
			int posB = posA - 1;
			get = 0;
			for (set<double>::iterator it=valueset.begin(); it!=valueset.end(); it++) {
				if (get == posA || get == posB) {
					median += *it;	
				}
				get++;
			}
			median /= 2;
		}
		
		/* 95% upper */
		pos = floor(count * 0.975);
		get = 0;
		for (set<double>::iterator it=valueset.begin(); it!=valueset.end(); it++) {
			if (get == pos) {
				upper = *it;	
			}
			get++;
		}
		
		v[0] = lower;
		v[1] = median;
		v[2] = upper;
	
	}
	
	else {
		v[0] = INF;
		v[1] = INF;
		v[2] = INF;
	}
	
	return v;
	
}

void Skyline::printMeans() {

	double m;
	for (int j=0; j < maxsteps; j++) {
		m = mean(j);
		if (m < INF)
			cout << skylineindices[maxiter][j] << "\t" << m << endl;
		else
			cout << skylineindices[maxiter][j] << "\tnan" << endl;
	}

}

/* ignore non-numbers */
double Skyline::mean(int row) {

	double mean = 0.0;
	multiset<double> valueset;
	
	/* move through iters and add to a set if a value exists for this row */
	for (int i=0; i <= iter; i++) {
		int size = skylineindices[i].size();
		if (row < size && !isnan(skylinevalues[i][row])) {
			valueset.insert(skylinevalues[i][row]);
		}
	}
			
	int count = valueset.size();
	int n = 0;
	
	if (count > 0) {
	
		/* mean */
		for (set<double>::iterator it=valueset.begin(); it!=valueset.end(); it++) {
			mean += *it;
			n++;
		}
		
		mean /= n;
			
	}
	
	else {
		mean = INF;
	}
	
	return mean;
	
}