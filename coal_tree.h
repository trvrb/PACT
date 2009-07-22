/* coal_tree.h
CoalescentTree class definition
This object stores and manipulates coalescent trees, rooted bifurcating trees with nodes mapped to time points
*/

#ifndef COALTREE_H
#define COALTREE_H

#include <map>
using std::map;
#include <set>
using std::set;
#include <vector>
using std::vector;
#include "tree.hh"

class CoalescentTree {

public:
	CoalescentTree(string,string);		// constructor, takes a parentheses string as input as well as an options string
										// can only be "migrate" or "beast"

	void pruneToTrunk();				// reduces CoalescentTree object to trunk
										// updates all relevant variables
	void pruneToLabel(int);				// reduces CoalescentTree object to only include a particular set of tips
										// updates all relevant variables
	void trimEnds(double,double);		// reduces CoalescentTree object to only those nodes between
										// time start and time stop
			
	void printTree();					// print indented tree
	void printTimeTree();				// print indented tree with coalescent times
	void printParenTree();				// print parentheses tree
	void printPaddedRuleList();			// print tree in Mathematica suitable format, 1->2,1->3, etc...
										// padded with extra nodes at coalescent time points
										// used with TreePlot
	void printRuleList();				// print tree in Mathematica rule list format with times included
										// used with Graphics primitives
										
	void printTrunkRatio();				// print proportion of tree that can trace its history forward
										// to present day samples
	void printLabelProportions();		// print list of proportions with label i										
	void printMigTotal();				// print overall migration rate across tree 
	void printMigRates();				// print list of posterior migration rates 
	void printCoalRates();				// print list of coalescent rates for each label
	void printTrunkRates();				// print list of coalescent rates for each label 
										// only considering coalescent events that include the trunk
	
	map<int,double> getCoalWeights();	// return map of coalescent weights for each label	
	map<int,int> getCoalCounts();		// return map of coalescent counts for each label
	map<int,double> getMigWeights();	// return map of forward migration weights for each label
	map<int,int> getMigCounts();		// return map of forward migration counts for each label
	map<int,double> getRevMigWeights();	// return map of backward migration weights for each label
	map<int,int> getRevMigCounts();		// return map of backward migration counts for each label	

	set<int> getLabelSet();				// returns labelSet
																	
	void NeSkyline();					// updates skyline with effective coalesent timescale Ne*tau
	void subRateSkyline();				// updates skyline with substitution rate on phylogeny	
										// takes mean rate across all concurrent lineages
	void divSkyline();					// updates skyline with diversity (pi) estimates
	void tmrcaSkyline();				// updates skyline with TMRCA estimates	
	void tajimaSkyline();				// updates skyline with estimates of Tajima's D
	void tcSkyline();					// updates skyline with time for each sampled lineage to coalesce with phylogeny trunk
	void labelSkyline(int);				// updates skyline with proportion of branches with label	
										// takes mean proportion across all concurrent lineages									
											
	vector<double> getSkylineIndex();	// go through skyline object and print medians
	vector<double> getSkylineValue();	// go through skyline object and print medians
																	
	void setStepSize(double);			// sets the window at which to take parameter values
	void pushTimesBack(double);			// push dates to agree with a most recent sample date at t
		
private:
	int leafCount;						// number of leafs in the tree (sample size n)
	int nodeCount;						// number of nodes in the tree
	set<int> leafSet;					// set of leaf nodes
	set<int> trunkSet;					// set of trunk nodes
	set<int> labelSet;					// set of distinct node labels
	
	double prunetime;					// time back from present with which to consider nodes as trunk
	double mp;							// multiplier to convert branch lengths into years
										// set externally
	
	tree<int> ctree;					// coalescent tree, has n leaf nodes (labelled 1..n) and n-1 internal nodes
										// root is labelled 0, other internal nodes continue with n+1
										// each node contains its label
	map<int,double>	bmap;				// maps a node to the length of the branch leading into the node
	map<int,double>	rmap;				// maps a node to the rate of the branch leading into the node
	map<int,double>	tmap;				// maps a node to the age of the node
	map<int,int>	lmap;				// maps a node to a numeric label
	map<int,int>	ancMap;				// maps a node to its parent	
	set<double> tlist;					// sorted list of times starting with root		
	bool blCheck;						// whether branch lengths exist in input tree
	bool brCheck;						// whether branch rates exist in input tree
	bool nlCheck;						// whether node labels exist in input tree	
	double stepsize;					// length of time between samples of parameter values
	
	void constructPaddedTree();			// takes a coalescent tree and pads with internal nodes at every event
	tree<int> paddedTree;				// padded with extra nodes at coalescent time points
	int paddedNodeCount;
	map<int,double>	paddedbmap;
	map<int,double>	paddedrmap;
	map<int,double> paddedtmap;
	map<int,int> 	paddedlmap;	
	map<int,int>	paddedAncMap;		// takes a node and maps it to its parent
		
	vector<double> skylineindex;
	vector<double> skylinevalue;
		
	tree<int> extractSubtree(set<int>);	// takes a set of node labels and walks up the coalescent tree
										// returning a tree object, tmap still works with this object											
	double getTreeLength(tree<int> &);	// takes a tree object and returns the total tree length, works with tmap
	
};

#endif