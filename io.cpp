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
#include "series.h"

IO::IO() {

	/* these will eventually be taken care of with a Parameters class */
	inputFile = "in.trees";
	outputPrefix = "out";

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
						ct.pushTimesBack(2007,2007);
					//	ct.section(2001.75,0.5,1);
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

/* go through problist and treelist and print highest posterior probability tree */
void IO::printHPTree() {

	string outputFile = outputPrefix + ".rules";
	cout << "Printing tree with highest posterior probability to " << outputFile << endl;

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

/* go through treelist and summarize coalescent statistics */
void IO::printStatistics() {

	string outputFile = outputPrefix + ".stats";
	cout << "Printing coalescent statistics to " << outputFile << endl;

	/* initializing output stream */
	ofstream outStream;
	outStream.open( outputFile.c_str(),ios::out);
	
	outStream << "statistic\t95% lower\tmean\t95% upper" << endl; 
	
	int n = treelist[0].getMaxLabel();

	// COALESCENCE /////////////////////
	for (int label = 1; label <= n; label++) {
	
		Series s;
		for (int i = 0; i < treelist.size(); i++) {
			double n = treelist[i].getCoalRate(label);
			s.insert(n);
		}
		outStream << "coal_" << label << "\t";
		outStream << s.quantile(0.025) << "\t" << s.mean() << "\t" << s.quantile(0.975) << endl;
	
	}
	
	// MIGRATION ///////////////////////
	for (int from = 1; from <= n; from++) {
		for (int to = 1; to <= n; to++) {	
			if (from != to) {

				Series s;
				for (int i = 0; i < treelist.size(); i++) {
					double n = treelist[i].getMigRate(from,to);
					s.insert(n);
				}
				outStream << "mig_" << from << "_" << to << "\t";
				outStream << s.quantile(0.025) << "\t" << s.mean() << "\t" << s.quantile(0.975) << endl;

			}
		}	
	}
	
	// TRUNK PROPORTIONS //////////////
	for (int label = 1; label <= n; label++) {

	Series s;
	for (int i = 0; i < treelist.size(); i++) {
		CoalescentTree ct = treelist[i];
		ct.pruneToTrunk();
		double n = ct.getLabelPro(label);
		s.insert(n);
	}
	outStream << "pro_" << label << "\t";
	outStream << s.quantile(0.025) << "\t" << s.mean() << "\t" << s.quantile(0.975) << endl;
	
	}
	
	
	outStream.close();

}