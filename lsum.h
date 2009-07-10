/* lsum.h
LabelSummary class definition
This object is for summarizing statistics across a tree.  Stores a vector of values associated with each label.
*/

#ifndef LSUM_H
#define LSUM_H

#include <map>
using std::map;
#include <set>
using std::set;
#include <vector>
using std::vector;

class LabelSummary {

public:
	LabelSummary(set<int>);							// constructor, takes a set of labels
													// values initialized to 0.0
													
	void incrWeights(map<int,double>);				// increment weights
	void incrCounts(map<int,int>);					// increment counts

	void printCounts();								// print counts	
	void printRates();								// print counts / weights
					
private:
	set<int> labels;								// set of labels
	map<int,double> weights;						// maps labels to values
	map<int,int> counts;							// maps labels to values
	
};

#endif