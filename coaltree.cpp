/* coaltree.cpp
Member function definitions for CoalescentTree class
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

#include <string>
using std::string;

#include <set>
using std::set;

#include <vector>
using std::vector;

#include <cmath>

#include "coaltree.h"
#include "node.h"
#include "tree.hh"

/* Constructor function to initialize private data */
/* Takes NEWICK parentheses tree as string input */
CoalescentTree::CoalescentTree(string paren) {

	string::iterator is;
	tree<Node>:: iterator it, jt;

	// STRIP PAREN STRING ///////////
	// strip spaces from paren string
	// strip & and following character, replace following : with =
	// assumes migration events follow the format [&M 5 3:8.49916e-05]
	is = paren.begin() ;
	bool mig = false;
	while (is < paren.end()) {
		if (*is == ' ')
			is = paren.erase(is);
		else if (*is == '&') {
			is = paren.erase(is);
			is = paren.erase(is);
			mig = true;
		}
		else if (*is == ':' && mig) {
			*is = '=';
			mig = false;
		}
		else
			is++;
	}

	// GATHER TIPS ////////////////
	// read in node names, filling tips vector
	// names can only be 0-9 A-Z a-z
	// exported tree renames tips with consecutive numbering starting at 1
	// go through paren string and collect tips, at the same time replace names with matching numbers in paren string
	// if first character of name is a number and the second character is a letter, assume first character is label

	vector<Node> tipsList;
	int current = 1;
	int stringPos = 0;
	string thisString = "";
	
	while (stringPos < paren.length()) {
		
		char thisChar = paren[stringPos];
		
		if ( (thisChar >= 'A' && thisChar <= 'Z') || (thisChar >= 'a' && thisChar <= 'z') || (thisChar >= '0' && thisChar <= '9') ) {
			thisString += thisChar;
		} 	  	
				
		else if (thisChar == ':' && thisString.length() > 0) {
							
			/* nodetree update */	
			Node thisNode(current);
			thisNode.setName(thisString);
			
			// label is the first character of node string, incremented by 1
			// only attempt this if first character is number and second character is letter
			// otherwise set to 1
			if ( (thisString[0] >= '0' && thisString[0] <= '9') &&
					( (thisString[1] >= 'A' && thisString[1] <= 'Z') || (thisString[1] >= 'a' && thisString[1] <= '1') ) ) {
				thisNode.setLabel(atoi(thisString.substr(0,1).c_str()) + 1);
			}
			
			thisNode.setLeaf(true);
			
			tipsList.push_back(thisNode);

			/* replace name with number */	
			stringstream out;
			out << current;
			paren = paren.substr(0,stringPos - thisString.size()) + out.str() + paren.substr(stringPos,paren.length());
			
			/* move counter back */
			/* need to take into acount the length of the digits */
			stringPos -= thisString.size() - (out.str()).length() + 1;		
			
			thisString = "";
			current++;
		
		}
		
		else {
			thisString = "";
		}
		
		stringPos++;
		
	}
		
	// STARTING TREE /////////////////
	// construct starting point for tree (multifurcation from root)
	it = nodetree.set_head(tipsList[0]);
	for(int i = 1; i < tipsList.size(); i++) {
		it = nodetree.insert_after(it, tipsList[i]);
   	}
  	
  	// CONSTRUCT TREE /////////////////////
	// read parentheses string from left to right, stop when a close parenthesis is encountered
	// push the left and right nodes onto their own branch
	// replace parenthesis string with their parent node ((1,2),3)  --->  (4,3)
		
	// end when all parentheses have been eliminated
	while (paren.at(0) == '(') {
	
//		cout << paren << endl;
			
		int left, right, from, to, openParen, closeParen, openMig, closeMig;	
		double leftLength, rightLength, migLength;
		stringPos = 0;
		thisString = "";
					
		for ( is=paren.begin(); is < paren.end(); is++ ) {

			if (*is == '(') {
				openParen = stringPos;
				openMig = stringPos;
			}
			if ( (*is >= '0' && *is <= '9') || (*is >= 'A' && *is <= 'Z') || (*is >= 'a' && *is <= 'z') || *is == '.' || *is == '-' ) {
				thisString += *is;
			}	
			else {
			
				bool stringExists;
				if (thisString.length() > 0)
					stringExists = true;
				else
					stringExists = false;
				
				// branch length
				if (( *is == '[' || *is == ',' || *is == ')' ) && stringExists) {		
					leftLength = rightLength;
					rightLength = atof(thisString.c_str());
				}					
				
				// node number
				if (*is == ':' && stringExists) {	
					left = right;
					right = atoi(thisString.c_str());		
				}

				if (*is == ',') {
					openMig = stringPos;
				}		
				
				// MIGRATION EVENTS ////////////////////
				
				/* need to extend ctree here */
				/* can only deal with migration events that effect a tip node */
				/* this section is only called when brackets follow a tip node */
				if (*is == '=' && stringExists) {
				
					/* grabbing migration event */
					string labelString = thisString;
					labelString.erase(0,1);
					from = atoi(labelString.c_str()) + 1;
					labelString = thisString;
					labelString.erase(1,1);		
					to = atoi(labelString.c_str()) + 1;
					
				}				
				
				if (*is == ']') {
				
					closeMig = stringPos; 
				
					migLength = atof(thisString.c_str());
		//			cout << tempN << " mig from " << from << " to " << to << ", at " << migLength << endl;
					
					// push child node back by distance equal to migLength
					it = findNode(right);	
					(*it).setLength( rightLength - migLength );
		
					// create new intermediate node
					Node migNode(current);
					migNode.setLabel(to);
					migNode.setLength(migLength);
					
					// wrap this new node so that it inherits the old node
					nodetree.wrap(it,migNode);		
							
					/* replace parenthesis with new node label */	
					/* code is set up to deal with the situation of two labels before a parenthesis */
					stringstream out;
					out << current << ":" << migLength;
					paren = paren.substr(0,openMig+1) + out.str() + paren.substr(closeMig+1,paren.length());
					
					current++;			
					break;
				
				}
								
				thisString = "";
				
			}				
			
			// COALESCENT EVENTS //////////////////
			
			if (*is == ')') { 
						
				closeParen = stringPos; 
				
				// append a new node
				// append this new node with two branches (left node and right node)
								
				tree<Node>:: iterator iterLeft, iterRight, iterNew;
				iterLeft = findNode(left);
				iterRight = findNode(right);	

				(*iterLeft).setLength(leftLength);
				(*iterRight).setLength(rightLength);				
				
				Node newNode(current);
				newNode.setLabel( (*iterLeft).getLabel() );
			
				iterNew = nodetree.wrap(iterLeft,newNode);		
				nodetree.move_after(iterLeft,iterRight);
				
				/* replace parenthesis with new node label */		
				stringstream out;
				out << current;
				paren = paren.substr(0,openParen) + out.str() + paren.substr(closeParen+1,paren.length());
				
				current++;			
				break;
					
				
			}
	
		stringPos++;
			
		}
	
	}
		
	
	// adding branch length to the parent node's time to get the node's time
	for (it = nodetree.begin(); it != nodetree.end(); it++) {
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
		it++;
	}
	
	
	/* pushing the most recent sample up to time = 0 */
	pushTimesBack(0);
	
}

/* push dates to agree with a most recent sample date at endTime */
void CoalescentTree::pushTimesBack(double endTime) {
		
	// need to adjust times by this amount
	double diff = endTime - getPresentTime();
		
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); it++) {
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
		for (tree<Node>::leaf_iterator lit = nodetree.begin_leaf(); lit != nodetree.end_leaf(); lit++) {
			if ((*lit).getTime() < oldestSample) { 
				oldestSample = (*lit).getTime(); 
			}
		}
		
		double mp = (endTime - startTime) / (presentTime - oldestSample);
		
		// go through tree and multiple lengths by mp	
		for (it = nodetree.begin(); it != nodetree.end(); it++) {
			double l = (*it).getLength();
			(*it).setLength(l * mp);
		}	
		
		// update times in tree
		for (it = nodetree.begin(); it != nodetree.end(); it++) {
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
		
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); it++) {
		double t = (*it).getTime();
		(*it).setTime(t + diff);
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
    		it++;
    	}
    }
            
	reduce();
			
}


/* reduces a tree to samples with a single label */
void CoalescentTree::pruneToLabel(int label) {

	/* start by finding all tips with this label */
	set<int> labelset; 
	tree<Node>::iterator it, jt;
	
	for (it = nodetree.begin(); it != nodetree.end(); it++) {
		if ( (*it).getLabel() == label && (*it).getLeaf() ) {
		
			/* move up tree adding nodes to label set */
			jt = it;
			while (nodetree.is_valid(jt)) {
				labelset.insert( (*jt).getNumber() );
				jt = nodetree.parent(jt);
			}
		
		}
	}
			
	/* erase other nodes from the tree */
	it = nodetree.begin();
	while(it != nodetree.end()) {
		if (labelset.end() == labelset.find( (*it).getNumber() )) {
			it = nodetree.erase(it);
		}
		else {
    		it++;
    	}
    }
        
	reduce();
				
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
			/* this pruning causes a leaf node to become an internal node */
			if ((*it).getTime() > stop && (*jt).getTime() <= stop) {
			
				(*it).setTime( stop );
				(*it).setLength( (*it).getTime() - (*jt).getTime() );
				(*it).setLeaf(false);
				nodetree.erase_children(it);
				it = nodetree.begin();
			
			}
			
			/* if node > start and parent < start, push parent up to start */
			/* and reparent anc[node] to be a child of root */
			/* neither node nore anc[node] can be root */
			else if ((*it).getTime() > start && (*jt).getTime() < start && nodetree.depth(it) > 1) {	
				(*jt).setTime(start);
				(*jt).setLength(0.0);
				(*jt).setInclude(false);
				Node tempNode(-1);
				nodetree.move_ontop(nodetree.append_child(nodetree.begin(),tempNode),jt);
				it = nodetree.begin();
			}
			
			else {
				it++;
			}
		
		}
		
		else {
    		it++;
    	}
    }
    
    /* second pass for nodes < start */
    it = nodetree.begin();
    (*it).setTime(start);
    (*it).setLength(0.0);
	(*it).setInclude(false);    
	while(it != nodetree.end()) {	
		if ((*it).getTime() < start && nodetree.depth(it) > 0) {
			it = nodetree.erase(it);
		}
		else {
    		it++;
    	}
    }
               
	// go through tree and update lengths based on times
	for (it = nodetree.begin(); it != nodetree.end(); it++) {
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
			
			it = newtree.begin();
			jt = nodetree.begin();	
			newtree.replace(it,jt);			
		
			newtree.insert(newtree.begin(),tempNode);
			
		}
	}

	newtree.erase(newtree.begin());
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
		/* this pruning causes a leaf node to become an internal node */
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
    		it++;
    	}
    	
    }
    
	/* erase other nodes from the tree */
	it = nodetree.begin();
	while(it != nodetree.end()) {
		if (sliceset.end() == sliceset.find( (*it).getNumber() )) {
			it = nodetree.erase(it);
		}
		else {
    		it++;
    	}
    }    
    
    /* peel back trunk */
	for (it = nodetree.begin(); it != nodetree.end(); it++) {
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
	nodetree.move_after(nodetree.begin(),++nodetree.begin());
	nodetree.erase(nodetree.begin());
	(*nodetree.begin()).setLength(0.0);

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
	for (it = nodetree.begin(); it != nodetree.end(); it++) {
		tset.insert( (*it).getTime() );
	}
	
	/* pad tree with extra nodes, make sure there is a node at each time slice correspoding to coalescent event */
	it = nodetree.begin();
	while(it != nodetree.end()) {
	
		/* finding what the correct depth of the node should be */
		int newDepth = -1;
		for (is = tset.begin(); is != tset.find( (*it).getTime() ); is++) {
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
		
		it++;
	
	}
			  	
}

/* Print indented tree */
void CoalescentTree::printTree() { 

	int rootdepth = nodetree.depth(nodetree.begin());
		
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); it++) {
		for(int i=0; i<nodetree.depth(it)-rootdepth; ++i) 
			cout << "  ";
		cout << (*it).getNumber();
		if ((*it).getName() != "") { 
			cout << " " << (*it).getName();
		}
		cout << " (" << (*it).getTime() << ")";
		cout << " [" << (*it).getLabel() << "]";			
		cout << " {" << (*it).getLength() << "}";		
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
*/	
void CoalescentTree::printRuleList(string outputFile) {

	/* initializing output stream */
	ofstream outStream;
	outStream.open( outputFile.c_str(),ios::out);

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
		it++;
	}
	
	/* print leaf nodes */
	/* a node may be a leaf on the current tree, but not a leaf on the original tree */
	for (it = nodetree.begin(); it != nodetree.end(); it++) {
		if ((*it).getLeaf()) {
			outStream << (*it).getNumber() << " ";
		}
	}
	outStream << endl;	
	
	/* print trunk nodes */
	for (it = nodetree.begin(); it != nodetree.end(); it++) {
		if ((*it).getTrunk()) {
			outStream << (*it).getNumber() << " ";
		}
	}
	outStream << endl;	
			
	/* print the tree in rule list (Mathematica-ready) format */
	/* print only upward links */
	for (it = nodetree.begin(); it != nodetree.end(); it++) {		// increment past root
		jt = nodetree.parent(it);
		if (nodetree.is_valid(jt)) {
			outStream << (*it).getNumber() << "->" << (*jt).getNumber() << " ";
		}
	}
	outStream << endl;
	
	
	/* print mapping of labels in Mathematica format */
	for (it = nodetree.begin(); it != nodetree.end(); it++) {	
		outStream << (*it).getNumber() << "->" << (*it).getLabel() << " ";
	}
	outStream << endl;
	
	/* print mapping of nodes to coordinates */
  	int count = 0;
	for (it = nodetree.begin(); it != nodetree.end(); it++) {
		if (nodetree.depth(it) == 0) {
			count = 0;
		}
		outStream << (*it).getNumber() << "->{" << (*it).getTime() << "," << count << "} ";	
		count++;
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
	it++;
	
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
		it++;
		
	}
	
	cout << endl;

	
}


/* most recent node in tree, will always be a leaf */
double CoalescentTree::getPresentTime() {
	
	double t = (*nodetree.begin()).getTime();
	for (tree<Node>::leaf_iterator it = nodetree.begin_leaf(); it != nodetree.end_leaf(); it++) {
		if ((*it).getTime() > t) { 
			t = (*it).getTime(); 
		}
	}
	return t;

}

/* most ancient node in tree */
double CoalescentTree::getRootTime() {

	double t = (*nodetree.begin()).getTime();
	for (tree<Node>::leaf_iterator it = nodetree.begin_leaf(); it != nodetree.end_leaf(); it++) {
		if ((*it).getTime() < t) { 
			t = (*it).getTime(); 
		}
	}
	return t;
}

/* amount of time it takes for all samples to coalesce */
double CoalescentTree::getTMRCA() {
	return getPresentTime() - getRootTime();
}

/* number of labels 1 to n */
int CoalescentTree::getMaxLabel() {

	double n = 0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); it++ ) {
		if ( (*it).getLabel() > n ) {
			n = (*it).getLabel();
		}
	}	
	return n;

}

/* number of leaf nodes */
int CoalescentTree::getLeafCount() {

	double n = 0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); it++ ) {
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
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); it++ ) {
		if ( (*it).getInclude() ) {
			length += (*it).getLength();
		}
	}	
	return length;

}

/* length of the tree with label l */
double CoalescentTree::getLength(int l) {

	double length = 0.0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); it++ ) {
		if ( (*it).getInclude() && (*it).getLabel() == l ) {
			length += (*it).getLength();
		}
	}	
	return length;

}

/* get proportion of tree with label */
double CoalescentTree::getLabelPro(int l) { 
		
	return getLength(l) / getLength();
	
}

/* proportion of tree that can trace its history forward to present day samples */
/* trunk traced back from the last 1/100 of the time width */
double CoalescentTree::getTrunkPro() { 

	double totalLength = getLength();

	double trunkLength = 0.0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); it++ ) {
		if ( (*it).getInclude() && (*it).getTrunk()) {
			trunkLength += (*it).getLength();
		}
	}	
	
	return trunkLength / totalLength;
    	
}

/* returns the count of coalescent events */
int CoalescentTree::getCoalCount() {

	/* count coalescent events, these are nodes with two children */
	int count = 0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); it++) {
		if ((*it).getInclude() && nodetree.number_of_children(it) == 2) {		
			count++;
		}
	}
	return count;

}

/* returns the count of coalescent events with label */
int CoalescentTree::getCoalCount(int l) {

	/* count coalescent events, these are nodes with two children */
	int count = 0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); it++) {
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
		for (it = nodetree.begin(); it != nodetree.end(); it++) {
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
double CoalescentTree::getCoalWeight(int l) {

	// setting step to be 1/1000 of the total length of the tree
	double start = getRootTime();
	double stop = getPresentTime();
	double step = (stop - start) / (double) 1000;
		
	// step through tree counting concurrent lineages
	double weight = 0.0;
	for (double t = start; t <= stop; t += step) {
	
		int lineages = 0;
		tree<Node>::iterator it, jt;
		for (it = nodetree.begin(); it != nodetree.end(); it++) {
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

double CoalescentTree::getCoalRate(int l) {
	return getCoalCount(l) / getCoalWeight(l);
}

/* returns the count of migration events over entire tree */
int CoalescentTree::getMigCount() {

	/* count migration events, these are nodes in which the parent label differs from child label */
	tree<Node>::iterator it, jt;
	int count = 0;
	for (it = nodetree.begin(); it != nodetree.end(); it++) {
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
int CoalescentTree::getMigCount(int from, int to) {

	/* count migration events, these are nodes in which the parent label differs from child label */
	tree<Node>::iterator it, jt;
	int count = 0;
	for (it = nodetree.begin(); it != nodetree.end(); it++) {
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
double CoalescentTree::getMigRate(int from, int to) {
	return getMigCount(from,to) / getLength(to);
}

/* return mean of (2 * time to common ancestor) for every pair of leaf nodes */
double CoalescentTree::getDiversity() {

	double div = 0.0;
	int count = 0;

	/* iterating over every pair of leaf nodes */
	tree<Node>::leaf_iterator it, jt, kt;
	for (it = nodetree.begin_leaf(); it != nodetree.end_leaf(); it++) {
		for (jt = it; jt != nodetree.end_leaf(); jt++) {
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
double CoalescentTree::getDiversity(int l) {

	double div = 0.0;
	int count = 0;

	/* iterating over every pair of leaf nodes */
	tree<Node>::leaf_iterator it, jt, kt;
	for (it = nodetree.begin_leaf(); it != nodetree.end_leaf(); it++) {
		for (jt = it; jt != nodetree.end_leaf(); jt++) {
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
	for (it = nodetree.begin_leaf(); it != nodetree.end_leaf(); it++) {
		for (jt = it; jt != nodetree.end_leaf(); jt++) {
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
	for (it = nodetree.begin_leaf(); it != nodetree.end_leaf(); it++) {
		for (jt = it; jt != nodetree.end_leaf(); jt++) {
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


/* removes extraneous nodes from tree */
void CoalescentTree::reduce() {

	tree<Node>::iterator it, jt, kt;

	/* removing pointless nodes, ie nodes that have no coalecent
	events or migration events associated with them */
	for (it = nodetree.begin(); it != nodetree.end(); it++) {
		jt = nodetree.parent(it);
		if (nodetree.is_valid(jt) && nodetree.number_of_children(it) == 1) {	// no coalescence
			kt = nodetree.child(it,0);
			if ( (*it).getLabel() == (*kt).getLabel() ) {						// no migration
//				cout << "it = " << *it << ", kt = " << *kt << endl;
				(*kt).setLength( (*kt).getLength() + (*it).getLength() );	
 				nodetree.reparent(jt,it);										// push child node up to be sibling of node
 				nodetree.erase(it);												// erase node									
				it = nodetree.begin();
			}
		}
	}

}

/* returns maximium node associated with a node in the tree */
int CoalescentTree::getMaxNumber() {

	int n = 0;
	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); it++) {
		if ((*it).getNumber() > n) {
			n = (*it).getNumber();
		}
	}
	return n;

}

/* renumber tree via preorder traversal starting from n */
int CoalescentTree::renumber(int n) {

	for (tree<Node>::iterator it = nodetree.begin(); it != nodetree.end(); it++) {
		(*it).setNumber(n);
		n++;
	}
	return n;

}

/* given a number, returns iterator to associated node, or if not found, returns iterator to end of tree */
tree<Node>::iterator CoalescentTree::findNode(int n) {
	
	tree<Node>::iterator it;
	
	for (it = nodetree.begin(); it != nodetree.end(); it++) {
		if ((*it).getNumber() == n)
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
