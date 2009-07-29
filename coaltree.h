/* coaltree.h
CoalescentTree class definition
This object stores and manipulates coalescent trees, rooted bifurcating trees with nodes mapped to time points
*/

#ifndef CTREE_H
#define CTREE_H

#include <vector>
#include "node.h"
#include "tree.hh"

class CoalescentTree {

public:
	CoalescentTree(string);					// constructor, takes a parentheses string as input
											// can only be "migrate" or "beast"

	// TREE MANIPULATION
	void pushTimesBack(double);				// push dates to agree with a most recent sample date at t
	void pushTimesBack(double,double);		// oldest sample and most recent sample	
	void pruneToTrunk();					// reduces CoalescentTree object to trunk
	void pruneToLabel(int);					// reduces CoalescentTree object to only include a particular set of tips
	void trimEnds(double,double);			// reduces CoalescentTree object to only those nodes between
											// time start and time stop	
	void section(double,double,double);		// break tree up into section
	void timeSlice(double);					// reduces CoalescentTree to all the ancestors of time slice
	void padTree();							// TODO: fix this
											// pads CoalescentTre with additional nodes at each coalescent event
											// included mainly for compatibility with TreePlot	

	// TREE STRUCTURE
	void printTree();						// print indented tree with coalescent times
	void printRuleList();					// print tree in Mathematica rule list format with times included
											// used with Graphics primitives
	void printParen();						// TODO: migration events
											// print parentheses tree										

	// BASIC STATISTICS
	double getPresentTime();				// returns most recent time in tree
	double getRootTime();					// returns most ancient time in tree
	double getTMRCA();						// span of time in tree
	int getMaxLabel();						// returns the highest label present
	int getLeafCount();						// returns the count of leaf nodes in tree
	int getNodeCount();						// returns the total number of nodes in tree


	// LABEL STATISTICS		
	double getLength();						// return total tree length
	double getLength(int);					// return length with this label
	vector<double> getLengths();	
	vector<double> getLabelPro();			// proportion of tree with each label
	double getTrunkPro();					// proportion of tree that can trace its history from present day samples
	
	// COALESCENT STATISTICS
	// problem with weight calculation for sectioned data
	int getCoalCount();						// total count of coalescent events on tree
	int getCoalCount(int);					// count of coalescent events involving label on tree	
	double getCoalWeight();					// total opportunity for coalescence on tree
	double getCoalWeight(int);				// total opportunity for coalescence on tree
	double getCoalRate();
	double getCoalRate(int);
	vector<double> getCoalRates();		

	// MIGRATION STATISTICS
	int getMigCount();
	int getMigCount(int,int);	
	double getMigRate();					// normalized by 'from' length
	double getMigRate(int,int);
	vector<double> getMigRates();
	
	// DIVERSITY STATISTICS	
	double getDiversity();					// return mean of (2 * time to common ancestor) for every pair of leaf nodes
	double getDiversity(int);				// diversity only involving a particular label
	double getDiversityWithin();			// diversity where both samples have the same label
	double getDiversityBetween();			// diversity where both samples have different labels
	double getFst();						// Fst = (divBetween - divWithin) / divBetween
	double getTajimaD();					// return D = pi - S/a1, where pi is diversity, S is the total tree length, 
											// and a1 is a normalization factor
									
private:
	tree<Node> nodetree;					// linked tree containing Node objects		
										
	// HELPER FUNCTIONS
	void reduce();							// goes through tree and removes inconsequential nodes	
	int getMaxNumber();						// return larger number in tree
	int renumber(int);						// renumbers tree in preorder traversal starting from int 
											// returning 1 greater than the max in the tree
	tree<Node>::iterator findNode(int);		// return iterator to a Node in nodetree based upon matching number
											// if not found, returns iterator to end of tree
												
	tree<Node>::iterator commonAncestor(tree<Node>::iterator,tree<Node>::iterator);
	
};

#endif