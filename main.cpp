/* 	(P)osterior (A)nalysis of (C)oalescent (T)rees 
	Copyright 2009 Trevor Bedford <bedfordt@umich.edu>

	This program is designed to interpret and manipulate labeled phylogenetic trees.  Statistics 
	regarding the structured coalescent may be calculated.

*/

/*	Agenda:
	Fix padTree
	Include migration events in PrintParen
	Is there a way to split coal_tree?
		Seems like a good split would be printing / information functions vs. manipulation functions.
		Also, there are a number of low level functions that would be nice to add.  Things like move_branch.
		What sort of primitives would be necessary to make this work?
	Build command-line options.
	Tests?  Error checking in general.
*/

/* Tree class: Copyright 2001-2006 Kasper Peeters <kasper.peeters@aei.mpg.de>

*/

#include <iostream>
#include <fstream>
using std::ifstream;
using std::cout;
using std::endl;

#include <string>
using std::string;

#include <vector>
using std::vector;

// Template class for standard library style tree
#include "tree.hh"

// Class for coalescent nodes within a tree object
#include "node.h"

// Extension of the tree class to deal specifically with coalescent trees
#include "coaltree.h"

// Collects a series of measurements, usually from multiple trees
#include "series.h"

int main(int argc, char* argv[]) {				// arguments passed from the command line
		
//	double a = atof(argv[1]);
//	double b = atof(argv[2]);
			
	string parenString;
	ifstream inputFile ("trees.txt");
	if (inputFile.is_open()) {
		while (! inputFile.eof() ) {
			getline (inputFile,parenString);
			
			if (parenString.size() > 1) {
	
				CoalescentTree ct(parenString);
				ct.pushTimesBack(2002,2007);

		//		ct.pruneToTrunk();					
		//		ct.pruneToLabel(2);
		//		ct.trimEnds(2002.75,2003.23);
		//		ct.timeSlice(2004);
		//		ct.section(2001.75,0.5,1);
	
		//		ct.printTree();
	
		//		ct.printRuleList();
		
		//		for (int i = 1; i <= 6; i++) {
		//			cout << ct.getCoalWeight(i) << endl;
		//		}
	
				Series s;
				s.insert(1);
				s.insert(2);
				s.insert(3);
				s.insert(4);
				s.insert(5);		
				cout << s.quantile(0.5) << endl;
					
			}
			
		}
		inputFile.close();
	}
	else
		cout << "Unable to open tree file, expecting trees.txt";
	
	return 0;
}


