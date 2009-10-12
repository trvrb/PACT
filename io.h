/* io.h
Copyright 2009 Trevor Bedford <bedfordt@umich.edu>
IO class definition
This object reads a BEAST or Migrate treefile and performs calculations on the resulting vector of
CoalescentTrees.

Currently, fills a vector with CoalescentTrees.  This is memory-intensive, but faster and more elegant 
than doing multiple readthroughs of the tree file.  CoalescentTrees vector takes about 10X the memory of
the corresponding tree file.

Uses Parameter object to figure out which operations to perform.
*/

#ifndef IO_H
#define IO_H

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "coaltree.h"
#include "param.h"

class IO {

public:
	IO();									// constructor takes a input file
	
	void treeManip();						// perform tree manipulation operations
	void printTree();						// print highest posterior tree to .rules
	void printStatistics();					// print coalescent statistics to .stats
	void printTips();						// print tip statistics to .tips
	void printSkylines();					// print skyline values to .skylines
																			
private:

	Parameters param;						// parameters object, read from in.param

	string inputFile;						// complete name of input tree file
	string outputPrefix;					// prefix for output files .rules and .stats

	vector<CoalescentTree> treelist;		// vector of coalescent trees
	vector<double> problist;				// vector of assocatied probabilities

};

#endif