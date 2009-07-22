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
*/

/*	
	Times runs forward starting at 0, unless othewise indicated.	


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

#include "tree.hh"
#include "coal_tree.h"
#include "coal_tree.cpp"
#include "lsum.h"
#include "lsum.cpp"
#include "skyline.h"
#include "skyline.cpp"

#define INF pow(double(10),double(100)) // 10^100 (~infinity)

int main(int argc, char* argv[]) {				// arguments passed from the command line
		
//	double a = atof(argv[1]);
//	double b = atof(argv[2]);
//	cout << "a = " << a << endl;
//	cout << "b = " << b << endl;	
	
	Skyline skl;
	
	string parenString;
	ifstream inputFile ("trees.txt");
	if (inputFile.is_open()) {
		while (! inputFile.eof() ) {
			getline (inputFile,parenString);
			
			if (parenString.size() > 1) {
	
				CoalescentTree ct(parenString, "beast");
	//			LabelSummary lsum(ct.getLabelSet());

/*
				for (double i = a; i > 0; i -= 1.0 ) {
				
					double j;
					if (i - 0.1666 > 0)
						j = i - 0.1666;
					else
						j = 0.0;
						
					CoalescentTree temptree = ct;
					temptree.trimEnds(j,i);
					lsum.incrWeights(temptree.getRevMigWeights());
					lsum.incrCounts(temptree.getRevMigCounts());
					
				
				}
*/

		//		lsum.incrWeights(ct.getRevMigWeights());
		//		lsum.incrCounts(ct.getRevMigCounts());

				ct.pushTimesBack(2004);

		//		ct.pruneToTrunk();					
		//		ct.pruneToLabel(1);
		//		ct.trimEnds(a,b);
		
				ct.printRuleList();
		//		ct.printTimeTree();
		//		ct.printMigTotal();
		//		ct.printMigRates();
		//		ct.printCoalRates();
		//		ct.printTrunkRates();

		//		ct.printLabelProportions();
		
		//		lsum.printRates();
		
		//		ct.setStepSize(0.001);
		//		ct.NeSkyline();
		//		ct.labelSkyline(a);													
		//		skl.appendIteration(ct.getSkylineIndex(),ct.getSkylineValue());
			}
			
		}
		inputFile.close();
	}
	else
		cout << "Unable to open tree file, expecting trees.txt";


//	skl.printMeans();
	
	return 0;
}


