/* coaltree.h
Copyright 2009-2010 Trevor Bedford <bedfordt@umich.edu>
CoalescentTree class definition
This object stores and manipulates coalescent trees, rooted bifurcating trees with nodes mapped to time points
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

#ifndef CTREE_H
#define CTREE_H

#include <set>
using std::set;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "tree.hh"
#include "node.h"
#include "rng.h"

class CoalescentTree {

public:
	CoalescentTree(string);					// constructor, takes a parentheses string as input
											// starts with most recent sample set at time = 0
											// sharing a most recent sample time ensures skyline calculations 
											// will work properly

	// TREE MANIPULATION
	void pushTimesBack(double);				// push dates to agree with a most recent sample date at t
	void pushTimesBack(double,double);		// oldest sample and most recent sample	
	void reduceTips(double);					// reduces tree to ancestors of a subset of tips
	void renewTrunk(double);				// renews definition of trunk, working back from all recent tips
	void renewTrunkRandom(double);			// renews definition of trunk, working back from a random tip
	void pruneToTrunk();					// reduces CoalescentTree object to trunk
	void pruneToLabel(string);				// reduces CoalescentTree object to only include a particular set of tips
	void pruneToName(string);				// reduces CoalescentTree object to only include a particular tip	
	void pruneToTime(double,double);		// reduces CoalescentTree object to only include tips in a certain time frame
	void collapseLabels();					// sets all labels in CoalescentTree to 1
	void trimEnds(double,double);			// reduces CoalescentTree object to only those nodes between
											// time start and time stop	
	void sectionTree(double,double,double);	// break tree up into sections
	void timeSlice(double);					// reduces CoalescentTree to all the ancestors of time slice
	void trunkSlice(double);				// removes from CoalescentTree all descendents of trunk at slice of time
	void leafSlice(double, double);			// reduces CoalescentTree to ancestors of leafs in a window of time
	void padTree();							// TODO: fix this
											// pads CoalescentTre with additional nodes at each coalescent event
											// included mainly for compatibility with TreePlot	
	void rotateLoc(double);					// rotate X&Y locations around origin

	// TREE STRUCTURE
	void printTree();						// print indented tree with coalescent times			
	void printRuleList(string);				// print to file name in Mathematica rule list format
											// used with Graphics primitives
	void printParen();						// TODO: migration events
											// print parentheses tree										

	// BASIC STATISTICS
	double getPresentTime();				// returns most recent time in tree
	double getRootTime();					// returns most ancient time in tree
	double getTMRCA();						// span of time in tree
//	int getMaxLabel();						// returns the highest label present
	int getLeafCount();						// returns the count of leaf nodes in tree
	int getNodeCount();						// returns the total number of nodes in tree

	// LABEL STATISTICS		
	double getLength();						// return total tree length
	double getLength(string);				// return length with this label
	double getLabelPro(string);				// return proportion of tree with label
	double getTrunkPro();					// proportion of tree that can trace its history from present day samples
	set<string> getLabelSet();				// return labelset
	
	// COALESCENT STATISTICS
	// problem with weight calculation for sectioned data
	int getCoalCount();						// total count of coalescent events on tree
	int getCoalCount(string);				// count of coalescent events involving label on tree	
	double getCoalWeight();					// total opportunity for coalescence on tree
	double getCoalWeight(string);			// total opportunity for coalescence on tree
	double getCoalRate();
	double getCoalRate(string);

	// MIGRATION STATISTICS
	int getMigCount();
	int getMigCount(string,string);	
	double getMigRate();			
	double getMigRate(string,string);
	
	// DIVERSITY STATISTICS	
	double getDiversity();					// return mean of (2 * time to common ancestor) for every pair of leaf nodes
	double getDiversity(string);			// diversity only involving a particular label
	double getDiversityWithin();			// diversity where both samples have the same label
	double getDiversityBetween();			// diversity where both samples have different labels
	double getFst();						// Fst = (divBetween - divWithin) / divBetween
	double getTajimaD();					// return D = pi - S/a1, where pi is diversity, S is the total tree length, 
											// and a1 is a normalization factor

	// LOCATION STATISTICS
	double getMeanX();						// return mean X location across all tips in the tree
	double getMeanY();						// return mean Y location across all tips in the tree
	vector<double> getTipsX();				// returns a vector of double for X position of every tip in tree
	vector<double> getTipsY();				// returns a vector of double for Y position of every tip in tree	
	
	// RATE STATISTICS
	double getMeanRate();					// return mean rate across all tips in the tree	

	// TIP STATISTICS
	vector<string> getTipNames();			// returns vector of tip names
	double getTime(string);
	string getLabel(string);
	double timeToTrunk(string);				// time it takes for a named tip to coalesce with the trunk
	
									
private:
	RNG rgen;								// random number generator
	tree<Node> nodetree;					// linked tree containing Node objects	
	set<string> labelset;					// set of all label names
										
	// HELPER FUNCTIONS
	string initialDigits(string);			// return initial digits in a string, 34ATZ -> 34, 3454 -> 0
	void reduce();							// goes through tree and removes inconsequential nodes	
	void peelBack();						// removes excess root from tree
	void adjustCoords();					// sets coords in Nodes to allow tree drawing
	int getMaxNumber();						// return larger number in tree
	int renumber(int);						// renumbers tree in preorder traversal starting from int 
											// returning 1 greater than the max in the tree
	tree<Node>::iterator findNode(int);		// return iterator to a Node in nodetree based upon matching number
	tree<Node>::iterator findNode(string);	// return iterator to a Node in nodetree based upon matching name	
											// if not found, returns iterator to end of tree
												
	tree<Node>::iterator commonAncestor(tree<Node>::iterator,tree<Node>::iterator);
	
};

#endif