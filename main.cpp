/* 	Takes an evolutionary tree in parentheses notation and outputs the same tree
	as a Mathematica compatible rule list

 	Input tree from command line, and output tree to standard out
 	
	Issues to be aware of:
	Tree tips must contain at least a single non-digit, otherwise they may be confused with the internal numeric ordering

*/

/*	Refactoring 7/9/09
	PaddedTree is kludgy, want a padTree function, otherwise everything should work identically
	Is there a way to split coal_tree?
		Seems like a good split would be printing / information functions vs. manipulation functions.
		Also, there are a number of low level functions that would be nice to add.  Things like move_branch.
		This will be especially useful if I want to build this into an MCMC.
		What sort of primitives would be necessary to make this work?
	Build in command-line options.
	Tests?  Error checking in general.
	Run time forwards.
	Always labeled.  Just sometimes labelled with all 0s.
*/

/*	
	Times runs forward starting at 0, unless othewise indicated.	
*/

/*
	Refactoring 7/24/09
	I have a lot of simple for loops to inspect nodes.  This assumes that nodes are labeled consecutively from 
	0 to nodeCount.  This is not always the case (specifically after pruning the tree).  Need to rewrite these operations
	to work with a set of nodes.  Or, might be easier to cleanup tree and maps to be consecutive.
*/

using namespace std;

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <map>
#include <vector>
#include <set>

// Template class for standard library style tree
#include "tree.hh"

// Class for coalescent nodes within a tree object
#include "node.h"
#include "node.cpp"

// Extension of the tree class to deal specifically with coalescent trees
#include "coaltree.h"
#include "coaltree.cpp"

// Measurements done on a CoalescentTree
#include "measurement.h"
#include "measurement.cpp"

// Collects a series of measurement objects, usually from multiple trees
#include "series.h"
#include "series.cpp"

// Soon to be deprecated
#include "lsum.h"
#include "lsum.cpp"
#include "skyline.h"
#include "skyline.cpp"

#define INF pow(double(10),double(100)) // 10^100 (~infinity)

int main(int argc, char* argv[]) {				// arguments passed from the command line
		
//	double a = atof(argv[1]);
//	double b = atof(argv[2]);
		
//	Skyline skl;
	
	string parenString;
	ifstream inputFile ("trees.txt");
	if (inputFile.is_open()) {
		while (! inputFile.eof() ) {
			getline (inputFile,parenString);
			
			if (parenString.size() > 1) {
	
				CoalescentTree ct(parenString);
				ct.pushTimesBack(0.2,9.98356);

		//		ct.pruneToTrunk();					
		//		ct.pruneToLabel(2);
		//		ct.trimEnds(2003,2004);
		//		ct.timeSlice(5);
		
				ct.printRuleList();
		//		ct.printTree();
		
				Measurement ms;
				ms.increment(ct.getLabelCoalCounts());
				ms.divideBy(ct.getLabelCoalOpp());
				ms.print();
		
			}
			
		}
		inputFile.close();
	}
	else
		cout << "Unable to open tree file, expecting trees.txt";
	
	return 0;
}


