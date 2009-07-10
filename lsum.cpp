/* lsum.cpp
Member function definitions for LabelSummary class
*/

#include "lsum.h"

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
LabelSummary::LabelSummary(set<int> ls) {
	
//	labels = ls;
	
	/* initializing value maps to 0.0 */
//	for ( set<int>::const_iterator lit = labels.begin(); lit != labels.end(); lit++ ) {
//		weights[*lit] = 0.0;
//		counts[*lit] = 0;
//	}
	
}

void LabelSummary::incrWeights(map<int,double> values) {
	for( map<int,double>::const_iterator iter = values.begin(); iter != values.end(); iter++ ) {
		weights[(*iter).first] += (*iter).second;
	}
}

void LabelSummary::incrCounts(map<int,int> values) {
	for( map<int,int>::const_iterator iter = values.begin(); iter != values.end(); iter++ ) {
		counts[(*iter).first] += (*iter).second;
	}
}

void LabelSummary::printCounts() {
	for( map<int,int>::const_iterator iter = counts.begin(); iter != counts.end(); iter++ ) {
//		cout << "label = " << *lit << ", count = " << counts[*lit] << ", weight = " << weights[*lit] << ", rate = " << counts[*lit] / weights[*lit] << endl;
		cout << counts[(*iter).first] << " ";
	}
	cout << endl;
}

void LabelSummary::printRates() {
	for( map<int,int>::const_iterator iter = counts.begin(); iter != counts.end(); iter++ ) {
//		cout << "label = " << *lit << ", count = " << counts[*lit] << ", weight = " << weights[*lit] << ", rate = " << counts[*lit] / weights[*lit] << endl;
		cout << counts[(*iter).first] / weights[(*iter).first] << " ";
	}
	cout << endl;
}
