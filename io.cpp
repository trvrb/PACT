/* io.cpp
Member function definitions for IO class
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

#include "io.h"
#include "coaltree.h"

IO::IO(string inputFile) {

	ifstream inStream;
	inStream.open( inputFile.c_str(),ios::out);

	string line;
	int pos;
	if (inStream.is_open()) {
		while (! inStream.eof() ) {
			getline (inStream,line);
			
			if (line.size() > 0) {
				if (line[0] != '#') {
			
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
						treelist.push_back(ct);
						cout << "tree " << treelist.size() << " read" << endl;
							
					}

				}
			}			
		}
		inStream.close();
	}
	else {
		throw runtime_error("input file not found");
	}
		
}

void IO::printHPTree(string outputFile) {

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
	
	treelist[index].printRuleList(outputFile);

}