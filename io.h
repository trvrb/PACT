/* io.h
Copyright 2009-2010 Trevor Bedford <bedfordt@umich.edu>
IO class definition
This object reads a BEAST or Migrate treefile and performs calculations on the resulting vector of
CoalescentTrees.

Currently, fills a vector with CoalescentTrees.  This is memory-intensive, but faster and more elegant 
than doing multiple readthroughs of the tree file.  CoalescentTrees vector takes about 10X the memory of
the corresponding tree file.

Uses Parameter object to figure out which operations to perform.
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
	void printPairs();						// pair pair statistics to .pairs
																			
private:

	Parameters param;						// parameters object, read from in.param

	string inputFile;						// complete name of input tree file
	string outputPrefix;					// prefix for output files .rules and .stats

	vector<CoalescentTree> treelist;		// vector of coalescent trees
	vector<double> problist;				// vector of assocatied probabilities

};

#endif