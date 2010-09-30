/* coaltree.cpp
Copyright 2009-2010 Trevor Bedford <bedfordt@umich.edu>
Member function definitions for CoalescentTree class
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

#include <iostream>
#include <sstream>
#include <fstream>
using std::ofstream;
using std::stringstream;
using std::cout;
using std::endl;
using std::flush;
using std::ios;
using std::fixed;

#include <string>
using std::string;

#include <set>
using std::set;

#include <vector>
using std::vector;

#include <stdexcept>
using std::runtime_error;
using std::out_of_range;

#include <cstdlib>
using std::atof;
using std::atoi;

#include <cmath>

#include "coaltree.h"
#include "node.h"
#include "tree.hh"

/* Constructor function to initialize private data */
/* Takes NEWICK parentheses tree as string input */
CoalescentTree::CoalescentTree(string paren) {

	string::iterator is;
	tree<Node>:: iterator it, jt;
	
	// make sure that parentheses are matched, counting '(' and counting ')'
	int leftcount = 0;
	int rightcount = 0;
	for ( is=paren.begin(); is < paren.end(); ++is ) {
		if (*is == '(') { leftcount++; }
		if (*is == ')') { rightcount++; }
	}
	if (leftcount != rightcount) {
		throw runtime_error("unmatched parentheses in in.trees");
	}

	// STARTING TREE /////////////////
	// starting point as single root node
	Node rootNode = Node(0);
	it = nodetree.set_head(rootNode);
	
	// WALK THROUGH NEWICK STRING ////
	// collect a substring, stop at ( ) , :
	
	string nameOrLength = "";
	string bracketed = "";
	int nodeCount = 1;
	bool lengthCheck = false;
	bool bracketCheck = false;
	bool braceCheck = false;
	
	// fill 'nameOrLength' with names and branch lengths
	// 'bracketed' placeholder for bracketed strings
	for (is = paren.begin(); is < paren.end(); ++is ) {
									
		// OUTSIDE OF BRACKETS
		// branch tree, update names, updates branch lengths
		if (!bracketCheck) {
		
			// filling nameOrLength
			if ( (*is >= '0' && *is <= '9') || (*is >= 'A' && *is <= 'Z') || (*is >= 'a' && *is <= 'z') || *is == '.' || *is == '-' || *is == '_' ) {
				nameOrLength += *is;
			}	
		
			// : --> name node, keep pointer where it is, prime loop to update a length, set as tip
			if (*is == ':') {
			
				if (nameOrLength.length() > 0) {
					(*it).setName(nameOrLength);
					(*it).setLeaf(true);
					(*it).setLabel(initialDigits(nameOrLength));
					if(initialDigits(nameOrLength) != "0") {
						labelset.insert(initialDigits(nameOrLength));
					}
					nameOrLength = "";
				}
				
				lengthCheck = true;
				
			}	
							
			if ( (*is == '[' || *is == '(' || *is == ')' || *is == ',') && nameOrLength.length() > 0) {
			
				//  update node length, if lengthCheck is flagged
				if (lengthCheck) {
					(*it).setLength(atof(nameOrLength.c_str()));
					lengthCheck = false;			
				}
				
				//  update node name, assuming branch lengths are absent, set as tip
				else {
					(*it).setName(nameOrLength);
					(*it).setLeaf(true);
					(*it).setLabel(initialDigits(nameOrLength));
					if(initialDigits(nameOrLength) != "0") {
						labelset.insert(initialDigits(nameOrLength));
					}
				}
				
				nameOrLength = "";
				
			}			
			
			// ( --> add child node, move pointer to this child node
			if (*is == '(') {
				Node thisNode(nodeCount);
				it = nodetree.append_child(it,thisNode);
				nodeCount++;
			}
	
			// , --> add sister node, move pointer to this sister node
			if (*is == ',') {
				Node thisNode(nodeCount);
				it = nodetree.insert_after(it,thisNode);
				nodeCount++;
			}
			
			// ) --> move pointer to parent node, need to inherit state when moving up the tree
			if (*is == ')') {
				string childLabel = (*it).getLabel();
				it = nodetree.parent(it);
				(*it).setLabel(childLabel);
			}		
		
		}
		
		// prime bracketed
		if (*is == '[') { 
			bracketCheck = true; 
			bracketed = "";
		}		
		
		// INSIDE OF BRACKETS
		// update labels, add migration events
		if (bracketCheck) {
		
			// fill bracketed
			if (bracketCheck && *is != '[' && *is != ']') {
				bracketed += *is;
			}

			if (*is == '{') { braceCheck = true; }
			if (*is == '}') { braceCheck = false; }			
			
		
			if (*is == ']' || (*is == ',' && braceCheck == false)) {
			
	//			cout << bracketed << endl;
			
				// fill 4 strings delimited by ' ', ':' and '=' 
				string::iterator ib;	
				string stringOne = "";
				string stringTwo = "";
				string stringThree = "";
				string stringFour = "";
				int counter = 1;
				
				ib = bracketed.begin();
				while (ib != bracketed.end()) {
					if (*ib != '&' && *ib != '{' && *ib != '}' && *ib != '"') {		// ignore these completely
						if (*ib != ' ' && *ib != '=' && *ib != ':' && *ib != ',') {				// ignore these and increment counter
							if (counter == 1) { stringOne += *ib; }
							if (counter == 2) { stringTwo += *ib; }
							if (counter == 3) { stringThree += *ib; }
							if (counter == 4) { stringFour += *ib; }
						}
						else {
							++counter;
						}
					}
				++ib;
				}
				
				// MIGRATION
				// insert an additional node up the tree
				if (stringOne == "M") {
				
					int fromInt = atoi(stringTwo.c_str()) + 1;
					int toInt = atoi(stringThree.c_str()) + 1;
					double migLength = atof(stringFour.c_str());
					
					stringstream fromStream;
					fromStream << fromInt;
					string from = fromStream.str();
					
					stringstream toStream;
					toStream << toInt;
					string to = toStream.str();					
										
					// push current node back by distance equal to migLength
					double newLength = (*it).getLength() - migLength;
					(*it).setLength(migLength);
			//		(*it).setLabel(to); 				// this should happen automatically
		
					// create new intermediate node
					Node migNode(nodeCount);
					migNode.setLabel(from);
					labelset.insert(from);
					migNode.setLength(newLength);
					nodeCount++;
					
					// wrap this new node so that it inherits the old node
					it = nodetree.wrap(it,migNode);						
					
				}
				
				// LOCATION
				// label current node
				if (stringOne == "states") {
					string loc = stringTwo.c_str();
					(*it).setLabel(loc);
					labelset.insert(loc);
				}
				
				// ANTIGENIC
				if (stringOne == "antigenic") {
					double xloc = atof(stringTwo.c_str());
					(*it).setX(xloc);
					double yloc = atof(stringThree.c_str());
					(*it).setY(yloc);					
				}				
				
				
				bracketed = "";
				
				if (*is == ']') { 
					bracketCheck = false;
				}
				
			}
		
		}
			
	}
	
	// adding branch length to the parent node's time to get the node's time
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		jt = nodetree.parent(it);
		if (nodetree.is_valid(jt)) {
			double t = (*jt).getTime() + (*it).getLength();
			(*it).setTime(t);
		}
	}	
	  			  		
	/* go through tree and append to trunk set */
	/* only the last 1/100 of the time span is considered */
	double presentTime = getPresentTime();
	double trunkTime = presentTime / (double) 100;
	it = nodetree.begin();
	(*it).setTrunk(true);
	while(it != nodetree.end()) {
		/* find nodes at present */
		if ((*it).getTime() > presentTime - trunkTime) {
			jt = it;
			/* move up tree adding nodes to trunk set */
			while (nodetree.is_valid(jt)) {
				(*jt).setTrunk(true);
				jt = nodetree.parent(jt);
			}
		}
		++it;
	}
	
	/* pushing the most recent sample up to time = 0 */
	pushTimesBack(0);
		
}

/* return initial digits in a string, incremented by 1, return 0 on failure 34ATZ -> 35, 3454 -> 0 */
string CoalescentTree::initialDigits(string name) {

	// label is the first digit characters of node string
	int initial = -1;
	
	// if string contains a letter
	bool containsLetter = false;
	for (int i = 0; i < name.length(); ++i) {
		if ( (name[i] >= 'A' && name[i] <= 'Z') || (name[i] >= 'a' && name[i] <= 'z') ) {
			containsLetter = true;
		}
	}
					
	if (containsLetter) {
		// construct substring of initial numbers
		string labelString = "";
		int count = 0;
		while (name[count] >= '0' && name[count] <= '9') {
			labelString += name[count];
			++count;
		}
		// set label to this substring
		initial = atoi(labelString.c_str());
	}
	
	stringstream outStream;
	outStream << initial + 1;
	string out = outStream.str();
	
	return out;

}

/* push dates to agree with a most recent sample date at endTime */
void CoalescentTree::pushTimesBack(double endTime) {
		
	// need to adjust times by this amount
	double diff = endTime - getPresentTime();
		
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); ++it) {
		double t = (*it).getTime();
		(*it).setTime(t + diff);
	}	
		  	
}

/* push dates to agree with a most recent sample date at endTime and oldest sample date is startTime */
/* will fail if used on contempory samples */
void CoalescentTree::pushTimesBack(double startTime, double endTime) {
	
	tree<Node>::iterator it, jt;
	double presentTime = getPresentTime();
	
	if (startTime < endTime) {
	
		// STRETCH OR SHRINK //////////////	 
			 
		// find oldest sample	
		double oldestSample = presentTime;
		for (tree<Node>::leaf_iterator lit = nodetree.begin_leaf(); lit != nodetree.end_leaf(); ++lit) {
			if ((*lit).getTime() < oldestSample) { 
				oldestSample = (*lit).getTime(); 
			}
		}
		
		double mp = (endTime - startTime) / (presentTime - oldestSample);
		
		// go through tree and multiple lengths by mp	
		for (it = nodetree.begin(); it != nodetree.end(); ++it) {
			double l = (*it).getLength();
			(*it).setLength(l * mp);
		}	
		
		// update times in tree
		for (it = nodetree.begin(); it != nodetree.end(); ++it) {
			jt = nodetree.parent(it);
			if (nodetree.is_valid(jt)) {
				double t = (*jt).getTime() + (*it).getLength();
				(*it).setTime(t);
			}
		}	
	
	}

	// PUSH BACK /////////////////////

	// need to adjust times by this amount
	double diff = endTime - getPresentTime();
		
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); ++it) {
		double t = (*it).getTime();
		(*it).setTime(t + diff);
	}		
		 
}

/* old version of renewTrunk.  This peels back from all current nodes. */
//void CoalescentTree::renewTrunk(double t) {
//
//	/* go through tree and append to trunk set */
//	double presentTime = getPresentTime();
//	tree<Node>::iterator it, jt;
//	
//	for(it = nodetree.begin(); it != nodetree.end(); ++it) {
//		(*it).setTrunk(false);
//	}
//	
//	it = nodetree.begin();
//	(*it).setTrunk(true);
//	while(it != nodetree.end()) {
//		/* find nodes at present */
//		if ((*it).getTime() > presentTime - t) {
//			jt = it;
//			/* move up tree adding nodes to trunk set */
//			while (nodetree.is_valid(jt)) {
//				(*jt).setTrunk(true);
//				jt = nodetree.parent(jt);
//			}
//		}
//		++it;
//	}	
//				
//}

/* reduces a tree to a random subset of samples */
void CoalescentTree::reduceTips(double pro) {

	/* start by finding all tips with this label */
	set<int> tipset; 
	tree<Node>::iterator it, jt;
	
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		if ( rgen.uniform(0,1) < pro && (*it).getLeaf() ) {
		
			/* move up tree adding nodes to label set */
			jt = it;
			while (nodetree.is_valid(jt)) {
				tipset.insert( (*jt).getNumber() );
				jt = nodetree.parent(jt);
			}
		
		}
	}
			
	/* erase other nodes from the tree */
	it = nodetree.begin();
	while(it != nodetree.end()) {
		if (tipset.end() == tipset.find( (*it).getNumber() )) {
			it = nodetree.erase(it);
		}
		else {
    		++it;
    	}
    }
        
   	peelBack();     
	reduce();
				
}

/* reduces a tree to just its trunk, takes a single random sample from "current" tips and works backward from this */
void CoalescentTree::renewTrunk(double t) {

	/* go through tree and append to trunk set */
	double presentTime = getPresentTime();
	tree<Node>::iterator it, jt;
	
	/* count tips and set every node as non-trunk */
	int count = 0;	
	for(it = nodetree.begin(); it != nodetree.end(); ++it) {
		(*it).setTrunk(false);
		if ((*it).getTime() > presentTime - t && (*it).getLeaf()) {
			count++;
		}
	}

	int selection = rgen.uniform(0,count);
	count = 0;
	
	it = nodetree.begin();
	(*it).setTrunk(true);
	while(it != nodetree.end()) {
		/* find nodes at present */
		if ((*it).getTime() > presentTime - t && (*it).getLeaf()) {
			if (selection == count) {
				jt = it;
				/* move up tree adding nodes to trunk set */
				while (nodetree.is_valid(jt)) {
					(*jt).setTrunk(true);
					jt = nodetree.parent(jt);
				}
				break;
			}
			count++;
		}
		++it;
	}	
	
}

/* reduces a tree to just its trunk, takes most recent sample and works backward from this */
void CoalescentTree::pruneToTrunk() {
	
	/* erase other nodes from the tree */	
	tree<Node>::iterator it;
	it = nodetree.begin();
	while(it != nodetree.end()) {
		if ( !(*it).getTrunk() ) {
			it = nodetree.erase(it);
		}
		else {
    		++it;
    	}
    }
            
	reduce();
			
}


/* reduces a tree to samples with a single label */
void CoalescentTree::pruneToLabel(string label) {

	/* start by finding all tips with this label */
	set<int> nodeset; 
	tree<Node>::iterator it, jt;
	
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		if ( (*it).getLabel() == label && (*it).getLeaf() ) {
		
			/* move up tree adding nodes to label set */
			jt = it;
			while (nodetree.is_valid(jt)) {
				nodeset.insert( (*jt).getNumber() );
				jt = nodetree.parent(jt);
			}
		
		}
	}
			
	/* erase other nodes from the tree */
	it = nodetree.begin();
	while(it != nodetree.end()) {
		if (nodeset.end() == nodeset.find( (*it).getNumber() )) {
			it = nodetree.erase(it);
		}
		else {
    		++it;
    	}
    }
        
//	peelBack();     
	reduce();
				
}

/* reduces a tree to samples within a specific time frame */
void CoalescentTree::pruneToTime(double start, double stop) {

	/* start by finding all tips within this time frame */
	set<int> nodeset; 
	tree<Node>::iterator it, jt;
	
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		if ( (*it).getTime() > start && (*it).getTime() < stop && (*it).getLeaf() ) {
		
			/* move up tree adding nodes to label set */
			jt = it;
			while (nodetree.is_valid(jt)) {
				nodeset.insert( (*jt).getNumber() );
				jt = nodetree.parent(jt);
			}
		
		}
	}
			
	/* erase other nodes from the tree */
	it = nodetree.begin();
	while(it != nodetree.end()) {
		if (nodeset.end() == nodeset.find( (*it).getNumber() )) {
			it = nodetree.erase(it);
		}
		else {
    		++it;
    	}
    }
        
//	peelBack();     
	reduce();
				
}


/* sets all labels in tree to 1 */
void CoalescentTree::collapseLabels() {

	labelset.clear();
	labelset.insert("1");

	tree<Node>::iterator it, jt;
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		(*it).setLabel("1");
	}
							
}

/* trims a tree at its edges 

   		   |-------	 			 |-----
from  ------			 	to	--
   		   |----------		 	 |-----

*/
void CoalescentTree::trimEnds(double start, double stop) {
			
	/* erase nodes from the tree where neither the node nor its parent are between start and stop */
	tree<Node>::iterator it, jt;
	it = nodetree.begin();
	while(it != nodetree.end()) {	
	
		jt = nodetree.parent(it);
	
		if (nodetree.is_valid(jt)) {
	
			/* if node > stop and parent < stop, erase children and prune node back to stop */
			/* this pruning causes an internal node to become an leaf node */
			if ((*it).getTime() > stop && (*jt).getTime() < stop) {
			
				(*it).setTime( stop );
				(*it).setLength( (*it).getTime() - (*jt).getTime() );
				//(*it).setLeaf(false);
				(*it).setLeaf(true);
				nodetree.erase_children(it);
				it = nodetree.begin();
			
			}
			
			/* if node > start and parent < start, push parent up to start */
			/* and reparent anc[node] to be a child of root */
			/* neither node nore anc[node] can be root */
			else if ((*it).getTime() > start && (*jt).getTime() < start) {
			
				(*jt).setTime(start);
				(*jt).setLength(0.0);
				(*jt).setInclude(false);
				nodetree.move_after(nodetree.begin(),jt);
				it = nodetree.begin();
			
			}
			
			else {
				++it;
			}
		
		}
		
		else {
    		++it;
    	}
    }
        
    /* second pass for nodes < start */
    it = nodetree.begin();   
	while(it != nodetree.end()) {	
		if ((*it).getTime() < start) {
			it = nodetree.erase(it);
		}
		else {
    		++it;
    	}
    }
        
	// go through tree and update lengths based on times
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		jt = nodetree.parent(it);
		if (nodetree.is_valid(jt)) {
			(*it).setLength( (*it).getTime() - (*jt).getTime() );
		}
	}	               
               
	reduce();
					
}

/* cuts up tree into multiple sections */
void CoalescentTree::sectionTree(double start, double window, double step) {

	tree<Node>::iterator it, jt;
	tree<Node> holdtree = nodetree;
	int current = 1;

	double rootTime = getRootTime();
	double presentTime = getPresentTime();

	/* newtree holds the growing tree structure */
	tree<Node> newtree;
	Node tempNode(-1);	
	newtree.set_head(tempNode);

	/* move window forward in time, make sure there are nodes in this window */
	for (double t = start; t < presentTime; t += step) {
		if (t > rootTime) {
		
			// operations all affect nodetree
			nodetree = holdtree;
			trimEnds(t,t + window);	
			current = renumber(current);			// need unique node numbers
			
			//	printTree();
					
			// need to move multiple sibling branches
			// need four sibling iterators, to1 to2 from1 from2
			tree<Node>::sibling_iterator to1, to2, from1, from2;
	
			from1 = nodetree.begin();
			from2 = nodetree.begin();
			while ( nodetree.is_valid(nodetree.next_sibling(from2)) ){
				from2 = nodetree.next_sibling(from2);
			}

			to1 = newtree.begin();
			to2 = newtree.begin();
			while ( newtree.is_valid(newtree.next_sibling(to2)) ){
				to2 = newtree.next_sibling(to2);
			}
					
			newtree.merge(to1,to2,from1,from2,true);
	
		}
	}

	nodetree = newtree;
	
}

/* Reduces tree to just the ancestors of a single slice in time */
/* Used to calcuate diversity, TMRCA and Tajima's D at a particular time */
void CoalescentTree::timeSlice(double slice) {

	/* desire only nodes spanning the time slice */
	/* find these nodes and add them and their ancestors to a set */
	set<int> sliceset; 
	tree<Node>::iterator it, jt, kt;
	it = nodetree.begin();
	while(it != nodetree.end()) {	
	
		jt = nodetree.parent(it);
	
		/* if node > slice and parent < slice, erase children and prune node back to stop */
		/* this pruning causes an internal node to become a leaf node */
		if ((*it).getTime() > slice && (*jt).getTime() <= slice) {
		
			// adjusting node
			(*it).setTime( slice );
			(*it).setLength( (*it).getTime() - (*jt).getTime() );
			(*it).setLeaf(true);
			nodetree.erase_children(it);
			
			// move up tree adding nodes to sliceset
			jt = it;
			while (nodetree.is_valid(jt)) {
				sliceset.insert( (*jt).getNumber() );
				jt = nodetree.parent(jt);
			}
			
			it = nodetree.begin();
		
		}
		
		else {
    		++it;
    	}
    	
    }
    
	/* erase other nodes from the tree */
	it = nodetree.begin();
	while(it != nodetree.end()) {
		if (sliceset.end() == sliceset.find( (*it).getNumber() )) {
			it = nodetree.erase(it);
		}
		else {
	   		++it;
    	}
	}    
    
	peelBack();
	reduce();

}

/* Removes descendents of trunk from tree at a single slice in time */
/* Used to calcuate time to fixation */
void CoalescentTree::trunkSlice(double slice) {

	/* desire only nodes spanning the time slice */
	/* these nodes must be trunk nodes */
	tree<Node>::iterator it, jt, kt;
	it = nodetree.begin();
	while(it != nodetree.end()) {	
	
		jt = nodetree.parent(it);
	
		/* if node > slice AND parent < slice AND node is trunk AND parent is trunk */
		/* erase children and prune node back to stop */
		/* this pruning causes an internal node to become a leaf node */
		if ((*it).getTime() > slice && (*jt).getTime() <= slice && (*it).getTrunk() && (*jt).getTrunk()) {
		
			// adjusting node
			(*it).setTime( slice );
			(*it).setLength( (*it).getTime() - (*jt).getTime() );
			(*it).setLeaf(true);
			nodetree.erase_children(it);
		
			it = nodetree.begin();
		
		}
		
		else {
    		++it;
    	}
    	
    }
        
}

/* Reduces tree to just the ancestors of leaf nodes existing in a particular window of time */
/* Used to calcuate diversity, TMRCA and Tajima's D for a window of time */
void CoalescentTree::leafSlice(double start, double stop) {

	/* desire only nodes spanning the time slice */
	/* find these nodes and add them and their ancestors to a set */
	set<int> sliceset; 
	tree<Node>::iterator it, jt, kt;
	it = nodetree.begin();
	while(it != nodetree.end()) {	
	
		/* if node > start AND node < stop && node is leaf */
		/* add ancestors to sliceset */
		/* there will be no children to erase */
		if ((*it).getTime() > start && (*it).getTime() <= stop && (*it).getLeaf()) {
		
			// move up tree adding nodes to sliceset
			jt = it;
			while (nodetree.is_valid(jt)) {
				sliceset.insert( (*jt).getNumber() );
				jt = nodetree.parent(jt);
			}
				
		}
		
		++it;
    	
    }
    
	/* erase other nodes from the tree */
	it = nodetree.begin();
	while(it != nodetree.end()) {
		if (sliceset.end() == sliceset.find( (*it).getNumber() )) {
			it = nodetree.erase(it);
		}
		else {
    		++it;
    	}
    }    
    
	peelBack();
	reduce();

}


/* padded with extra nodes at coalescent time points */
/* causing problems with migration tree */
void CoalescentTree::padTree() { 

	int current = getMaxNumber() + 1;

	tree<Node>::iterator it, end, iterTemp, iterN;
	
	/* construct set of coalescent times */
	set<double>::const_iterator is;
	set<double> tset;
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		tset.insert( (*it).getTime() );
	}
	
	/* pad tree with extra nodes, make sure there is a node at each time slice correspoding to coalescent event */
	it = nodetree.begin();
	while(it != nodetree.end()) {
	
		/* finding what the correct depth of the node should be */
		int newDepth = -1;
		for (is = tset.begin(); is != tset.find( (*it).getTime() ); ++is) {
    		newDepth++;
    	}
    	
    	is++;
	
		if (newDepth > nodetree.depth(it)) {
		
			/* padding with number of nodes equal to the difference in depth levels */
			for(int i = 0; i < newDepth - nodetree.depth(it); i++) {

				Node newNode(current);
				newNode.setLabel( (*it).getLabel() );
				newNode.setTime( *is );
				newNode.setLength( *is - (*it).getTime() );
				
				nodetree.wrap(it,newNode);
	
				current++;
				it = nodetree.begin();
				
			}
	
		}
		
		++it;
	
	}
			  	
}

/* Print indented tree */
void CoalescentTree::printTree() { 

	int rootdepth = nodetree.depth(nodetree.begin());
		
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); ++it) {
		for(int i=0; i<nodetree.depth(it)-rootdepth; ++i) 
			cout << "  ";
		cout << (*it).getNumber();
		if ((*it).getName() != "") { 
			cout << " " << (*it).getName();
		}
		cout << " (" << (*it).getTime() << ")";
		cout << " [" << (*it).getLabel() << "]";			
		cout << " {" << (*it).getLength() << "}";	
		cout << " <" << (*it).getX() << "," << (*it).getY() << ">";			
		if ( !(*it).getInclude()) { 
			cout << " *";
		}		
		cout << endl << flush;
	}
		
}

/* print tree in Mathematica suitable format
Output is:
	leaf list
	trunk list
	tree rules
	label rules
	coordinate rules 
	tip name rules
*/	
void CoalescentTree::printRuleList(string outputFile) {

//	printTree();

	/* initializing output stream */
	ofstream outStream;
	outStream.open( outputFile.c_str(),ios::app);

	/* setting up y-axis ordering, x-axis is date */
	adjustCoords();
	
	tree<Node>::iterator it, jt;
	
	/* print leaf nodes */
	/* a node may be a leaf on the current tree, but not a leaf on the original tree */
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		if ((*it).getLeaf()) {
			outStream << (*it).getNumber() << " ";
		}
	}
	outStream << endl;	
	
	/* print trunk nodes */
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		if ((*it).getTrunk()) {
			outStream << (*it).getNumber() << " ";
		}
	}
	outStream << endl;	
			
	/* print the tree in rule list (Mathematica-ready) format */
	/* print only upward links */
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {		// increment past root
		jt = nodetree.parent(it);
		if (nodetree.is_valid(jt)) {
			outStream << (*it).getNumber() << "->" << (*jt).getNumber() << " ";
		}
	}
	outStream << endl;
	
	
	/* print mapping of nodes to labels */
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {	
		outStream << (*it).getNumber() << "->\"" << (*it).getLabel() << "\" ";
	}
	outStream << endl;
	
	/* print mapping of nodes to coordinates */
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		outStream << (*it).getNumber() << "->{" << fixed << (*it).getTime() << "," << (*it).getCoord() << "} ";	
	}
	outStream << endl;
		  	  	
	/* print mapping of nodes to names */
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {	
		if ((*it).getName() != "")
			outStream << (*it).getNumber() << "->\"" << (*it).getName() << "\" ";
	}
	outStream << endl;  	  	
	  	  	
	outStream.close();
	  	  	
}

/* Print parentheses tree */
void CoalescentTree::printParen() { 

	tree<Node>::post_order_iterator it;
	it = nodetree.begin_post();
   	
	int currentDepth = nodetree.depth(it);
	for (int i = 0; i < currentDepth; i++) { 
		cout << "("; 
	} 
	cout << (*it).getNumber() << ":" << (*it).getLength(); 
	++it;
	
	/* need to add a '(' whenever the depth increases and a ')' whenever the depth decreases */
	/* only print leaf nodes */
	while(it != nodetree.end_post()) {
		if (nodetree.depth(it) > currentDepth) { 
			cout << ", ("; 
			for (int i = 0; i < nodetree.depth(it) - currentDepth - 1; i++) { 
				cout << "("; 
			}
			if (nodetree.number_of_children(it) == 0) { 
				cout << (*it).getNumber() << ":" << (*it).getLength(); 
			}
		}
		if (nodetree.depth(it) == currentDepth) { 
			if (nodetree.number_of_children(it) == 0) { 
				cout << ", " << (*it).getNumber() << ":" << (*it).getLength(); ; 
			}
		}
		if (nodetree.depth(it) < currentDepth) {
			if (nodetree.number_of_children(it) == 0) { 
				cout << (*it).getNumber() << ":" << (*it).getLength(); 
				cout << ")";		
			}
			else {
				cout << ")";	
				cout << ":" << (*it).getLength();
			}
		}
		currentDepth = nodetree.depth(it);
		++it;
		
	}
	
	cout << endl;

	
}


/* most recent node in tree, will always be a leaf */
double CoalescentTree::getPresentTime() {
	
	double t = (*nodetree.begin()).getTime();
	for (tree<Node>::leaf_iterator it = nodetree.begin_leaf(); it != nodetree.end_leaf(); ++it) {
		if ((*it).getTime() > t) { 
			t = (*it).getTime(); 
		}
	}
	return t;

}

/* most ancient node in tree */
double CoalescentTree::getRootTime() {

	double t = (*nodetree.begin()).getTime();
	for (tree<Node>::leaf_iterator it = nodetree.begin_leaf(); it != nodetree.end_leaf(); ++it) {
		if ((*it).getTime() < t) { 
			t = (*it).getTime(); 
		}
	}
	return t;
}

/* amount of time it takes for all samples to coalesce */
double CoalescentTree::getTMRCA() {
	
	int leafcount = 0;
	for (tree<Node>::leaf_iterator it = nodetree.begin_leaf(); it != nodetree.end_leaf(); ++it) {
		leafcount++;
	}
	
	double tmrca = 0.0;
	if (leafcount > 1) {
		tmrca = getPresentTime() - getRootTime();
	}
	else {
		tmrca /= tmrca;
	}
	
	return tmrca;
	
}

/* number of labels 1 to n */
//int CoalescentTree::getMaxLabel() {
//
//	double n = 0;
//	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); ++it ) {
//		if ( (*it).getLabel() > n ) {
//			n = (*it).getLabel();
//		}
//	}	
//	return n;
//
//}

/* number of leaf nodes */
int CoalescentTree::getLeafCount() {

	double n = 0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); ++it ) {
		if ( (*it).getLeaf() ) {
			n++;
		}
	}	
	return n;

}

/* total number of nodes */
int CoalescentTree::getNodeCount() {
	return nodetree.size();
}

/* total length of the tree */
double CoalescentTree::getLength() {

	double length = 0.0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); ++it ) {
		if ( (*it).getInclude() ) {
			length += (*it).getLength();
		}
	}	
	return length;

}

/* length of the tree with label l */
double CoalescentTree::getLength(string l) {

	double length = 0.0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); ++it ) {
		if ( (*it).getInclude() && (*it).getLabel() == l ) {
			length += (*it).getLength();
		}
	}	
	return length;

}

/* get proportion of tree with label */
double CoalescentTree::getLabelPro(string l) { 

	return getLength(l) / getLength();
	
}

/* proportion of tree that can trace its history forward to present day samples */
/* trunk traced back from the last 1/100 of the time width */
double CoalescentTree::getTrunkPro() { 

	double totalLength = getLength();

	double trunkLength = 0.0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); ++it ) {
		if ( (*it).getInclude() && (*it).getTrunk()) {
			trunkLength += (*it).getLength();
		}
	}	
	
	return trunkLength / totalLength;
    	
}

set<string> CoalescentTree::getLabelSet() {
	return labelset;
}

/* returns the count of coalescent events */
int CoalescentTree::getCoalCount() {

	/* count coalescent events, these are nodes with two children */
	int count = 0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); ++it) {
		if ((*it).getInclude() && nodetree.number_of_children(it) == 2) {		
			count++;
		}
	}
	return count;

}

/* returns the count of coalescent events with label */
int CoalescentTree::getCoalCount(string l) {

	/* count coalescent events, these are nodes with two children */
	int count = 0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); ++it) {
		if ((*it).getInclude() && nodetree.number_of_children(it) == 2 && (*it).getLabel() == l ) {		
			count++;
		}
	}
	return count;

}

/* returns the opportunity for coalescence over the whole tree */
/* running this will padTree() may be faster and more accurate */
double CoalescentTree::getCoalWeight() {

	// setting step to be 1/1000 of the total length of the tree
	double start = getRootTime();
	double stop = getPresentTime();
	double step = (stop - start) / (double) 1000;
	
	// step through tree counting concurrent lineages
	double weight = 0.0;
	for (double t = start; t <= stop; t += step) {
	
		int lineages = 0;
		tree<Node>::iterator it, jt;
		for (it = nodetree.begin(); it != nodetree.end(); ++it) {
			jt = nodetree.parent(it);
			if ( (*it).getInclude() && nodetree.is_valid(jt) && (*it).getTime() >= t && (*jt).getTime() < t) {		
				lineages++;
			}
		}
		
		if (lineages > 0) {
			weight += ( ( lineages * (lineages - 1) ) / 2 ) * step;
		}
		
	}	
	
	return weight;

}

/* returns the opportunity for coalescence for label */
double CoalescentTree::getCoalWeight(string l) {

	// setting step to be 1/1000 of the total length of the tree
	double start = getRootTime();
	double stop = getPresentTime();
	double step = (stop - start) / (double) 1000;
		
	// step through tree counting concurrent lineages
	double weight = 0.0;
	for (double t = start; t <= stop; t += step) {
	
		int lineages = 0;
		tree<Node>::iterator it, jt;
		for (it = nodetree.begin(); it != nodetree.end(); ++it) {
			jt = nodetree.parent(it);
			if ( (*it).getInclude() && nodetree.is_valid(jt) && (*it).getTime() >= t && (*jt).getTime() < t && (*it).getLabel() == l ) {		
				lineages++;
			}
		}
				
		if (lineages > 0) {
			weight += ( ( lineages * (lineages - 1) ) / 2 ) * step;
		}
		
	}	
	
	return weight;

}

double CoalescentTree::getCoalRate() {
	return getCoalCount() / getCoalWeight();
}

double CoalescentTree::getCoalRate(string l) {
	return getCoalCount(l) / getCoalWeight(l);
}

/* returns the count of migration events over entire tree */
int CoalescentTree::getMigCount() {

	/* count migration events, these are nodes in which the parent label differs from child label */
	tree<Node>::iterator it, jt;
	int count = 0;
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		jt = nodetree.parent(it);
		if (nodetree.is_valid(jt)) {		
			if ( (*it).getInclude() && (*jt).getInclude() && (*it).getLabel() != (*jt).getLabel() ) {		
				count++;
			}
		}
	}
	return count;

}

/* returns the count of migration events from label to label */
int CoalescentTree::getMigCount(string from, string to) {

	/* count migration events, these are nodes in which the parent label differs from child label */
	tree<Node>::iterator it, jt;
	int count = 0;
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		jt = nodetree.parent(it);
		if (nodetree.is_valid(jt)) {
			if ( (*it).getInclude() && (*jt).getInclude() && (*it).getLabel() == to && (*jt).getLabel() == from ) {		
				count++;
			}
		}
	}
	return count;

}

/* returns the overall rate of migration */
double CoalescentTree::getMigRate() {
	return getMigCount() / getLength();
}

/* returns the rate of migration from label to label */
/* this is important to check on */
/* currently, this is set up as calculating the rate from working backwards in time */
/* i.e. the migration rate from 1->2 is measured from the count going backwards on 2->1 divided */
/* by the backward opportunity of 2 */
/* getMigCount(from,to) / getLength(to) */
/* this needs attention */
/* seems to match with empirical estimates with getMigCount(from,to) / getLength() */
/* seems wrong however */
double CoalescentTree::getMigRate(string from, string to) {
	return getMigCount(from,to) / getLength(to);
}

/* return mean of (2 * time to common ancestor) for every pair of leaf nodes */
double CoalescentTree::getDiversity() {

	double div = 0.0;
	int count = 0;

	/* iterating over every pair of leaf nodes */
	tree<Node>::leaf_iterator it, jt, kt;
	for (it = nodetree.begin_leaf(); it != nodetree.end_leaf(); ++it) {
		for (jt = it; jt != nodetree.end_leaf(); ++jt) {
			if ((*it).getInclude() && (*jt).getInclude() && it != jt) {
	
				/* find common ancestor and calculate time from it to jt via common ancestor */
				kt = commonAncestor(it,jt);
				div += ( (*it).getTime() - (*kt).getTime() ) + ( (*jt).getTime() - (*kt).getTime() );
				count++;
			
			}
		}
	}
	
	div /= (double) count;
	return div;
	
}

/* return mean of (2 * time to common ancestor) for pairs of leaf nodes with labels a and b */
double CoalescentTree::getDiversity(string l) {

	double div = 0.0;
	int count = 0;

	/* iterating over every pair of leaf nodes */
	tree<Node>::leaf_iterator it, jt, kt;
	for (it = nodetree.begin_leaf(); it != nodetree.end_leaf(); ++it) {
		for (jt = it; jt != nodetree.end_leaf(); ++jt) {
			if ((*it).getInclude() && (*jt).getInclude() && it != jt && (*it).getLabel() == l && (*jt).getLabel() == l ) {
	
				/* find common ancestor and calculate time from it to jt via common ancestor */
				kt = commonAncestor(it,jt);
				div += ( (*it).getTime() - (*kt).getTime() ) + ( (*jt).getTime() - (*kt).getTime() );
				count++;
			
			}
		}
	}
	
	div /= (double) count;
	return div;
	
}

/* return mean of (2 * time to common ancestor) for pairs of leaf nodes with identical labels */
double CoalescentTree::getDiversityWithin() {

	double div = 0.0;
	int count = 0;

	/* iterating over every pair of leaf nodes */
	tree<Node>::leaf_iterator it, jt, kt;
	for (it = nodetree.begin_leaf(); it != nodetree.end_leaf(); ++it) {
		for (jt = it; jt != nodetree.end_leaf(); ++jt) {
			if ((*it).getInclude() && (*jt).getInclude() && it != jt && (*it).getLabel() == (*jt).getLabel() ) {
	
				/* find common ancestor and calculate time from it to jt via common ancestor */
				kt = commonAncestor(it,jt);
				div += ( (*it).getTime() - (*kt).getTime() ) + ( (*jt).getTime() - (*kt).getTime() );
				count++;
			
			}
		}
	}
	
	div /= (double) count;
	return div;
	
}

/* return mean of (2 * time to common ancestor) for pairs of leaf nodes with different labels */
double CoalescentTree::getDiversityBetween() {

	double div = 0.0;
	int count = 0;

	/* iterating over every pair of leaf nodes */
	tree<Node>::leaf_iterator it, jt, kt;
	for (it = nodetree.begin_leaf(); it != nodetree.end_leaf(); ++it) {
		for (jt = it; jt != nodetree.end_leaf(); ++jt) {
			if ((*it).getInclude() && (*jt).getInclude() && it != jt && (*it).getLabel() != (*jt).getLabel() ) {
	
				/* find common ancestor and calculate time from it to jt via common ancestor */
				kt = commonAncestor(it,jt);
				div += ( (*it).getTime() - (*kt).getTime() ) + ( (*jt).getTime() - (*kt).getTime() );
				count++;
			
			}
		}
	}
	
	div /= (double) count;
	return div;
	
}

/* returns population subdivision Fst = (divBetween - divWithin) / divBetween */
double CoalescentTree::getFst() {

	double divWithin = getDiversityWithin();
	double divBetween = getDiversityBetween();
	double fst = (divBetween - divWithin) / divBetween;
	return fst;

}

/* return D = pi - S/a1, where pi is diversity, S is the total tree length, and a1 is a normalization factor */
/* expect D = 0 under neutrality */
double CoalescentTree::getTajimaD() {

	double div = getDiversity();
	double S = getLength();
	
	double a1 = 0.0;
	double a2 = 0.0;	
	int n = getLeafCount();
	for (int i = 1; i < n; i++) {
		a1 += 1 / (double) i;
		a2 += 1 / (double) (i*i);		
	}
				
	double e1 = (1.0/a1) * ((double)(n+1) / (3*(n-1)) - (1.0/a1));	
	double e2 = (1.0 / (a1*a1 + a2) ) * ( (double)(2*(n*n+n+3)) / (9*n*(n-1)) - (double)(n+2) / (n*a1) + a2/(a1*a1) );
	double denom = sqrt(e1*S + e2*S*(S-1));
	double tajima = (div - S/a1) / denom;	
	
	return tajima;

}


/* returns vector of tip names */
vector<string> CoalescentTree::getTipNames() {

	vector<string> names;
	for (tree<Node>::leaf_iterator lit = nodetree.begin_leaf(); lit != nodetree.end_leaf(); ++lit) {
		names.push_back( (*lit).getName() );
	}
	return names;

}

double CoalescentTree::getTime(string name) {
	tree<Node>::iterator it = findNode(name);	
	return (*it).getTime();
}

string CoalescentTree::getLabel(string name) {
	tree<Node>::iterator it = findNode(name);	
	return (*it).getLabel();
}

/* time it takes for a named tip to coalesce with the trunk */
double CoalescentTree::timeToTrunk(string name) {

	tree<Node>::iterator it, jt;
	it = findNode(name);
	jt = it;
	
	/* walk back from this node until trunk is reached */
	while ( !(*jt).getTrunk() ) {
		jt = nodetree.parent(jt);
	}
	
	return (*it).getTime() - (*jt).getTime();	

}

/* removes extraneous nodes from tree */
void CoalescentTree::reduce() {

	tree<Node>::iterator it, jt, kt;

	/* removing pointless nodes, ie nodes that have no coalecent
	events or migration events associated with them */
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		jt = nodetree.parent(it);
		if (nodetree.is_valid(jt)) {
			if (nodetree.number_of_children(it) == 1) {								// no coalescence	
				kt = nodetree.child(it,0);
				if ((*kt).getLabel() == (*it).getLabel()) { 						// mo migration
	//				cout << "it = " << *it << ", kt = " << *kt << endl;
					(*kt).setLength( (*kt).getLength() + (*it).getLength() );	
					nodetree.reparent(jt,it);										// push child node up to be sibling of node
					nodetree.erase(it);												// erase node									
					it = nodetree.begin();
				}
			}
		}
	}

}

/* peels back trunk. works from root forward, stopping when first split is reached */
void CoalescentTree::peelBack() {

	tree<Node>::iterator it, jt, kt;

	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		jt = nodetree.parent(it);
		if ( nodetree.is_valid(jt) && nodetree.number_of_children(it) == 1) {	
			kt = nodetree.child(it,0);
			(*kt).setLength( (*kt).getLength() + (*it).getLength() );	
			nodetree.reparent(jt,it);								// push child node up to be sibling of node
			nodetree.erase(it);										// erase node									
			it = nodetree.begin();
		}
		if (nodetree.number_of_children(it) == 2) {
			break;
		}
	}

	// adjust root    
	it = nodetree.begin();
	if (nodetree.number_of_children(it) == 1) {
		nodetree.move_after(nodetree.begin(),++nodetree.begin());
		nodetree.erase(nodetree.begin());
		(*nodetree.begin()).setLength(0.0);
	}
	
}	
	
void CoalescentTree::adjustCoords() {

	tree<Node>::iterator it, jt;

	/* reorder tree so that the bottom node of two sister nodes always has the most recent child more children */
	/* this combined with preorder traversal will insure the trunk follows a rough diagonal */
	it = nodetree.begin();
	while(it != nodetree.end()) {
		jt = nodetree.next_sibling(it);
		if (nodetree.is_valid(jt)) {
			int cit = nodetree.size(it);
			int cjt = nodetree.size(jt);
			if (cit > cjt) {
				nodetree.swap(jt,it);
				it = nodetree.begin();
			}
		}
		++it;
	}

	/* set coords of tips according to preorder traversal */
  	int count = 0;
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
	//	if (nodetree.depth(it) == 0) { count = 0; }		// this resets count for each subtree
		if ( (*it).getLeaf() ) {
			(*it).setCoord(count);	
			count++;
		}
	}
	
	/* revise coords of internal nodes according to postorder traversal */
  	tree<Node>::post_order_iterator post_it, post_jt, post_kt;
  	for (post_it = nodetree.begin_post(); post_it != nodetree.end_post(); post_it++) {
  		if (nodetree.number_of_children(post_it) == 1) {
  			post_jt = nodetree.child(post_it,0);
  			(*post_it).setCoord((*post_jt).getCoord());	
  		}  		  	
  		if (nodetree.number_of_children(post_it) == 2) {
  			post_jt = nodetree.child(post_it,0);
  			post_kt = nodetree.child(post_it,1);
  			double avg = ( (*post_jt).getCoord() + (*post_kt).getCoord() ) / (double) 2;
  			(*post_it).setCoord(avg);	
  		}
	}	

}	
	
/* returns maximium node associated with a node in the tree */
int CoalescentTree::getMaxNumber() {

	int n = 0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); ++it) {
		if ((*it).getNumber() > n) {
			n = (*it).getNumber();
		}
	}
	return n;

}

/* renumber tree via preorder traversal starting from n */
int CoalescentTree::renumber(int n) {

	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); ++it) {
		(*it).setNumber(n);
		n++;
	}
	return n;

}

/* given a number, returns iterator to associated node, or if not found, returns iterator to end of tree */
tree<Node>::iterator CoalescentTree::findNode(int n) {
	
	tree<Node>::iterator it;
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		if ((*it).getNumber() == n)
			break;
	}
	return it;
	
}

/* given a name, returns iterator to associated node, or if not found, returns iterator to end of tree */
tree<Node>::iterator CoalescentTree::findNode(string name) {
	
	tree<Node>::iterator it;
	for (it = nodetree.begin(); it != nodetree.end(); ++it) {
		if ((*it).getName() == name)
			break;
	}
	return it;
	
}

/* given two iterators, returns an iterator to their most recent common ancestor */
tree<Node>::iterator CoalescentTree::commonAncestor(tree<Node>::iterator ia, tree<Node>::iterator ib) {

	/* make a set */
	set<int> nodeSet;
	
	/* walk down from first node to root, appending to nodeSet */
	tree<Node>::iterator it;
	it = ia;
	while (nodetree.is_valid(it)) {
		nodeSet.insert( (*it).getNumber() );
		it = nodetree.parent(it);
	}
	
	/* walk down from second node, stopping when a member of nodeSet is encountered */
	it = ib;	
	while (nodetree.is_valid(it)) {
		if (nodeSet.end() == nodeSet.find( (*it).getNumber() )) {
			it = nodetree.parent(it);
		}
		else {
			break;
		}
	}

//	cout << "a = " << (*ia).getNumber() << ", b = " << (*ib).getNumber() << ", anc = " << (*it).getNumber() << endl;
	
	return it;
	
}
