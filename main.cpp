/* 	(P)osterior (A)nalysis of (C)oalescent (T)rees 
	Copyright 2009 Trevor Bedford <bedfordt@umich.edu>
	
	tree.hh: Copyright 2001-2006 Kasper Peeters <kasper.peeters@aei.mpg.de>

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

/*
	This file is part of PACT.

    PACT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PACT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PACT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::ios;

#include <stdexcept>
using std::runtime_error;
using std::out_of_range;

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

	Series s;
	vector<CoalescentTree> treelist;
	vector<double> problist;
	
	string line;
	int pos;
	ifstream inputFile ("trees.txt");
	if (inputFile.is_open()) {
		while (! inputFile.eof() ) {
			getline (inputFile,line);
			
			// Catching log probabilities of trees
			string annoString;
			
			// migrate annotation
			annoString = "log likelihood = ";
			pos = line.find(annoString);

			if (pos >= 0) {
				string thisString;
				thisString = line.substr(pos+annoString.size());
				thisString.erase(thisString.find(']'));
				double ll = atof(thisString.c_str());
				problist.push_back(ll);
			}
			
			// beast annotation
			annoString = "[&lnP=";
			pos = line.find(annoString);

			if (pos >= 0) {
				string thisString;
				thisString = line.substr(pos+annoString.size());
				thisString.erase(thisString.find(']'));
				double ll = atof(thisString.c_str());
				problist.push_back(ll);				
			}
			
			// find first occurance of '(' in line
			pos = line.find('(');
			if (pos >= 0) {
			
				string paren = line.substr(pos);
				
				CoalescentTree ct(paren);
				
				ct.pushTimesBack(1968,2003);
		//		ct.pruneToTrunk();					
		//		ct.pruneToLabel(2);
		//		ct.trimEnds(2002.75,2003.23);
		//		ct.timeSlice(2004);
		//		ct.section(2001.75,0.5,1);
				
				treelist.push_back(ct);
	
		//		ct.printTree();

		//		ct.printRuleList();
		
		//		s.insert( ct.getCoalRate() );
		
				cout << "tree " << treelist.size() << " read" << endl;
					
			}
			
		}
		inputFile.close();
	}
	else {
		throw runtime_error("input file not understood");
	}

	// TREE OUTPUT /////////////////
	// output tree with highest likelihood
	// if likelihoods don't exist, output final tree
	
	int index;
	if (problist.size() == treelist.size()) {
		
		double ll = problist[0];
		index = 0;
		for (int i = 1; i < problist.size(); i++) {
			if (problist[i] > ll) {
				ll = problist[i];
				index = i;
			}
		}
		
	}
	else {
		index = treelist.size() - 1;
	}
	
	treelist[index].printRuleList("out.rules");



	return 0;
}


