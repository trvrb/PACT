/* io.h
IO class definition
This object reads a BEAST or Migrate treefile and performs calculations on the resulting vector of
CoalescentTrees.
*/

#ifndef IO_H
#define IO_H

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "coaltree.h"

class IO {

public:
	IO(string);								// constructor takes a input file
	
	void printHPTree(string);				// print highest posterior tree to output file
																			
private:
	vector<CoalescentTree> treelist;		// vector of coalescent trees
	vector<double> problist;				// vector of assocatied probabilities

};

#endif