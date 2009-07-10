/* coal_tree.cpp
Member function definitions for CoalescentTree class
*/

#include "coal_tree.h"
#include "tree.hh"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <sstream>

#define INF pow(double(10),double(100)) // 10^100 (~infinity)

/* Constructor function to initialize private data */
/* Takes NEWICK parentheses tree as string input */
CoalescentTree::CoalescentTree(string paren, string options) {

	prunetime = 0.01;
	stepsize = 0.1;	
	mp = 1.0;			// 1.927 for mig within coal, 1.47 for mig outside coal

	nlCheck = false;	
	
	string::iterator iterStr;

	// strip spaces from paren string
	// also & and next character
	iterStr=paren.begin() ;
	while (iterStr < paren.end()) {
		if (*iterStr == ' ')
			iterStr = paren.erase(iterStr);
		else if (*iterStr == '&') {
			iterStr = paren.erase(iterStr);
			iterStr = paren.erase(iterStr);
		}
		else
			iterStr++;
	}

	// check paren string for :, indicates branch lengths present
	blCheck = false;
	for ( iterStr=paren.begin() ; iterStr < paren.end(); iterStr++ ) {
		if (*iterStr == ':') {
			blCheck = true;
		}
	}
	
	// check paren string for [, indicates rates lengths present
	brCheck = false;
	for ( iterStr=paren.begin() ; iterStr < paren.end(); iterStr++ ) {
		if (*iterStr == '[') {
			brCheck = true;
		}
	}	

	if (options == "migrate") {
		nlCheck = true;
		brCheck = false;
	}

	// read in node names, filling tips
	// tips maps a node string to an int
	// names can only be 0-9 A-Z
	// exported tree renames tips with consecutive numbering starting at 1
	// ignore spaces
	// in BEAST trees I can just directly use the node ids
	
	map<string,int> tips; 
	leafCount = 0;
	
	string tempNode = "";
	bool nl = false;
	for ( iterStr=paren.begin() ; iterStr < paren.end(); iterStr++ ) {
		if (*iterStr >= '0' && *iterStr <= '9' || *iterStr >= 'A' && *iterStr <= 'Z') {
    		tempNode += *iterStr;
    	}
    	else {
    	    if ( (brCheck && *iterStr == '[' && tempNode.length()>0) ||
    	     (!brCheck && !nlCheck && blCheck && *iterStr == ':' && tempNode.length()>0) ||
    	     (!brCheck && !blCheck && tempNode.length()>0) ||
    	     (!brCheck && nlCheck && blCheck && *iterStr == ':' && tempNode.length()>0 && !nl) ) {
    	    	leafCount++;
    	    	if (options == "migrate") {
    				tips.insert( make_pair(tempNode,leafCount ) );
    				leafSet.insert(leafCount);	
    			}
   				if (options == "beast") {
    				tips.insert( make_pair(tempNode,atoi(tempNode.c_str()) ) );
    				leafSet.insert(atoi(tempNode.c_str()));	
    			}		
    			if (nlCheck) {													// filling lmap and labelSet if labels exist
    				lmap[ leafCount ] = atoi(tempNode.substr(0,1).c_str()) + 1;	// label is the first character of node string
    				labelSet.insert(atoi(tempNode.substr(0,1).c_str()) + 1);	// incremented by 1
    			}						
    			else {															// if not labelled set all labels to 1
     				lmap[ leafCount ] = 1;
    				labelSet.insert(1);   			
    			}
    		}
    		if (!brCheck && nlCheck && blCheck && *iterStr == '[')
    			nl = true;
    		if (!brCheck && nlCheck && blCheck && *iterStr == ']')
    			nl = false;		
    		tempNode = "";
    	}	
	}

	nodeCount = leafCount;
	int currentNode = nodeCount + 1;

	// construct starting point for tree (multifurcation from root)
	ctree.set_head(0);
	tmap[0] = 0.0;
	
	map<string,int>::iterator iter;
	for( iter = tips.begin(); iter != tips.end(); iter++ ) {
    	ctree.begin() = ctree.insert_after(ctree.begin(), iter->second);
		tmap[iter->second] = 0.0;
   	}
  	
	tlist.insert( 0.0 );



	// read parentheses string from left to right, stop when a close parenthesis is encountered
	// push the left and right nodes onto their own branch
	// replace parenthesis string with their parent node ((1,2),3)  --->  (4,3)
	// fill bmap
	
	tree<int>::iterator it, jt, end, iterLeft, iterRight, iterN, iterTemp;
	int leftNode, rightNode, openParen, closeParen, openMig, closeMig, stringPos;	
	vector<int> nodeList;
	vector<double> blList;
	string tempBl, tempBr;

	/* mapping of nodes to the label of their parents */
	map<int,int> lmaptemp;
	int leftLabel,rightLabel;	
	
	// end when all parentheses have been eliminated
	while (paren.at(0) == '(') {
	
		it = ctree.begin_post();
		end = ctree.end_post();
		stringPos = 0;
	
		// will need to keep track of this across iterations
		int thisNode, from, to;
		double migBranch;
	
		// fill nodeList until a close parenthesis is hit
		tempNode = "";
		tempBl = "";
		tempBr = "";
		nl = false; 
		
		for ( iterStr=paren.begin(); iterStr < paren.end(); iterStr++ ) {

			if (*iterStr == '(') {
				openParen = stringPos;
				openMig = stringPos;
			}
			if ( (*iterStr >= '0' && *iterStr <= '9') || (*iterStr >= 'A' && *iterStr <= 'Z') || *iterStr == '.' || *iterStr == 'e' || *iterStr == '-' ) {
				tempNode += *iterStr;
				tempBl += *iterStr;
				tempBr += *iterStr;
			}	
			else {
				if (brCheck) {	
					if (*iterStr == '[' && tempNode.length()>0) {
						if (tips.end() == tips.find(tempNode))				// not found, direct insert
							nodeList.push_back(atoi(tempNode.c_str()));
						else
							nodeList.push_back(tips[tempNode]);				// found, insert mapped int
					}
					if (*iterStr == ']' && tempNode.length()>0)
						rmap[nodeList.back()] = atof(tempBr.c_str());
					if (*iterStr != ':' && tempNode.length()>0)
						bmap[nodeList.back()] = atof(tempBl.c_str());
				}
				if (!brCheck && nlCheck && blCheck) {
					if (*iterStr == '[' && tempNode.length()>0 && !nl) {
						bmap[nodeList.back()] = atof(tempBl.c_str());
//						cout << "  bl = " << tempBl.c_str() << endl;
						nl = true;
					}
					if (*iterStr == ':' && tempNode.length()>0 && !nl) {
						if (tips.end() == tips.find(tempNode))				// not found, direct insert
							nodeList.push_back(atoi(tempNode.c_str()));
						else
							nodeList.push_back(tips[tempNode]);				// found, insert mapped int
//						cout << "node = " << tempNode.c_str() << endl;
					}
					if (*iterStr != ':' && tempNode.length()>0 && !nl) {
						bmap[nodeList.back()] = atof(tempBl.c_str());
//						cout << "  bl = " << tempBl.c_str() << endl;
					}
					if (*iterStr == ',') {
						openMig = stringPos;
					}	
					if (*iterStr == '[' && tempNode.length()==0) {
						nl = true;
					}						
					if (*iterStr == ']') {
					
						closeMig = stringPos; 
					
						migBranch = atof(tempNode.c_str());
					//	cout << thisNode << " mig from " << from << " to " << to << ", at " << migBranch << endl;
					
						/* wrapping a new node */
						tree<int>::iterator iterTo, iterFrom;
						it = ctree.begin();
						while(it!=end) {
							if (*it == thisNode) { 
								iterFrom = it; 
								break;
							}
							it++;
						}
						iterTo = ctree.wrap(iterFrom,currentNode);
						
						lmap[currentNode] = to;
						bmap[currentNode] = migBranch;
			
						/* this sometimes results in a negative branch length */
						/* however, I'm fairly certain this is a problem with Migrate */
						/* rather than a problem with my code */
						bmap[thisNode] = bmap[thisNode] - migBranch;		// this places the migration event 
																			// on the coalescent branch
				
					//	bmap[thisNode] = bmap[thisNode];					// this makes a new branch for the 
																			// migration event 					
						
						/* replace parenthesis with new node label */	
						/* code is set up to deal with the situation of two labels before a parenthesis */
						stringstream out;
						out << currentNode << ":" << bmap[currentNode];
						paren = paren.substr(0,openMig+1) + out.str() + paren.substr(closeMig+1,paren.length());
					//	cout << paren << endl;
						
						currentNode++;
						nodeCount++;
						iterStr=paren.begin();	// reset count
						break;
					
						nl = false;
					}
					
					/* this deals with what goes on between the brackets */
					/* need to extend ctree here */
					/* can only deal with migration events that effect a tip node */
					/* this section is only called when brackets follow a tip node */
					if (*iterStr == ':' && tempNode.length()>0 && nl) {
					
						/* grabbing migration event */
						string labelString = tempNode;
						thisNode = *--nodeList.end();
						labelString.erase(0,1);
						from = atoi(labelString.c_str()) + 1;
						labelString = tempNode;
						labelString.erase(1,1);		
						to = atoi(labelString.c_str()) + 1;
						lmaptemp[thisNode] = to;
						
					}
				}				
				if (!brCheck && !nlCheck && blCheck) {
					if (*iterStr == ':' && tempNode.length()>0) {
						if (tips.end() == tips.find(tempNode))				// not found, direct insert
							nodeList.push_back(atoi(tempNode.c_str()));
						else
							nodeList.push_back(tips[tempNode]);				// found, insert mapped int
					}
					if (*iterStr != ':' && tempNode.length()>0)
						bmap[nodeList.back()] = atof(tempBl.c_str());
				}
				if (!brCheck && !blCheck && tempNode.length()>0) {
					if (tips.end() == tips.find(tempNode))				// not found, direct insert
						nodeList.push_back(atoi(tempNode.c_str()));
					else
						nodeList.push_back(tips[tempNode]);				// found, insert mapped int
				}
				tempNode = "";
				tempBl = "";
				tempBr = "";
			}				
			if (*iterStr == ')') { 
				
				closeParen = stringPos; 
				leftNode = *----nodeList.end();
				rightNode = *--nodeList.end();
				
				// iterate over the tree
				// once right node is reached, append a new node
				// append this new node with two branches (left node and right node)
							
				while(it!=end) {
					if (*it == leftNode) { iterLeft = it; }
					if (*it == rightNode) { iterRight = it; }
					it++;
				}
								
				/* append a new node after root, append this node with two branches */
				/* replace with subtrees from before*/
				iterN = ctree.insert_after(iterRight,currentNode);
				iterTemp = ctree.append_child(iterN,0);
				ctree.replace(iterTemp,iterLeft);
				iterTemp = ctree.append_child(iterN,0);
				ctree.replace(iterTemp,iterRight);	
				
				/* update the label map with new node */
				leftLabel = lmap[leftNode];
				rightLabel = lmap[rightNode];
			
				if (leftLabel == rightLabel)
					lmap[currentNode] = leftLabel;
				else
					cout << "left and right labels don't match" << endl;
				
				lmaptemp.clear();
				
				/* erase copied nodes */
				ctree.erase(iterLeft);
				ctree.erase(iterRight);
						
				/* replace parenthesis with new node label */		
				stringstream out;
				out << currentNode;
				paren = paren.substr(0,openParen) + out.str() + paren.substr(closeParen+1,paren.length());
			//	cout << paren << endl;
				
				currentNode++;
				nodeCount++;
				iterStr=paren.begin();	// reset count
				break;
					
				
			}
	
		stringPos++;
			
		}
	
	}
	
	
	// removing null root and replacing node label of root with 0
	ctree.erase(ctree.begin_post());	
	*ctree.begin() = 0;
	lmap[0] = leftLabel;
		
	// go through bmap and multiply by mp	
	tree<int>::pre_order_iterator pre_it;
	for (pre_it = ++ctree.begin(); pre_it != ctree.end(); pre_it++) {
		bmap[*pre_it] *= mp;
	}	
	
	// use tree and bmap to get tmap
	double t;
	for (pre_it = ++ctree.begin(); pre_it != ctree.end(); pre_it++) {
		t = tmap[*ctree.parent(pre_it)] - bmap[*pre_it];
		tmap[*pre_it] = t;
	}

	// push ages up so most recent time = 0
	double mint = 0;
	for (pre_it = ctree.begin(); pre_it != ctree.end(); pre_it++) {
		if (tmap[*pre_it] < mint ) {
			mint = tmap[*pre_it];
		}	
	}
	
	// fill tlist
	for (pre_it = ctree.begin(); pre_it != ctree.end(); pre_it++) {
		tmap[*pre_it] -= mint;	
		tlist.insert(tmap[*pre_it]);
	}	
	
	/* construct ancMap */
	ancMap.clear();
	it = ctree.begin();
	end = ctree.end();
	while (it != end) {
		if (*it != 0) {
			ancMap[*it] = *ctree.parent(it);
		}
		it++;
	}
  	
	constructPaddedTree();  	
		
	/* go through tree and append to trunk set */
	it = ctree.begin();
	end = ctree.end();
	while(it!=end) {
		/* find leaf nodes at present */
		if (tmap[*it] < prunetime && *it <= leafCount) {
			jt = it;
			/* move up tree adding nodes to trunk set */
			while (*jt != 0) {
				trunkSet.insert(*jt);
				jt = ctree.parent(jt);
			}
		}
		it++;
	}
	trunkSet.insert(0);

}

/* push dates to agree with a most recent sample date at t */
void CoalescentTree::pushTimesBack(double time) {

	tree<int>::iterator it;

	// find most recent node
	double mint = *tlist.begin();
	for (it = ctree.begin(); it != ctree.end(); it++) {
		if (tmap[*it] < mint ) {
			mint = tmap[*it];
		}	
	}
	
	// need to adjust times by this amount
	double diff = time - mint;
		
	// modifying tmap and tlist
	tlist.clear();
	for (it = ctree.begin(); it != ctree.end(); it++) {
		tmap[*it] += diff;
		tlist.insert(tmap[*it]);
		
	}	
		  	
	constructPaddedTree();

}

/* padded with extra nodes at coalescent time points */
/* causing problems with migration tree */
void CoalescentTree::constructPaddedTree() { 

	/* don't want to modify the original tree */

	paddedTree = ctree;
  	paddedtmap = tmap;
  	paddedrmap = rmap;
  	paddedlmap = lmap;

	tree<int>::pre_order_iterator it, end, iterTemp, iterN;

	it = paddedTree.begin();
	end = paddedTree.end();
   
	if(!paddedTree.is_valid(it)) return;
	
	set<double>::const_iterator depth_it;
	int newDepth;
	int currentNode = nodeCount;

	/* pad tree with extra nodes, make sure there is a node at each time slice correspoding to coalescent event */
	while(it!=end) {
	
		/* finding what the correct depth of the node should be */
		newDepth = -1;
		for ( depth_it=tlist.end(); depth_it != tlist.find(paddedtmap[*it]); depth_it-- )
    		newDepth++;
    	
    	depth_it++;
	
		if (newDepth > paddedTree.depth(it)) {
		
			/* padding with number of nodes equal to the difference in depth levels */
			for(int i=0; i<newDepth-paddedTree.depth(it); i++) {
	
				double rtemp = paddedrmap[*it];
				int ltemp = paddedlmap[*it];
	
				/* need to make a new sibling for a subtree with less than the correct depth */
				iterN = paddedTree.insert(it, currentNode);
				paddedtmap[currentNode] = *depth_it;
				
				/* create a temporary child for this new node */
				iterTemp = paddedTree.append_child(iterN,0);
				
				/* replace this child with the subtree it */
				paddedTree.replace(iterTemp,it);	
			
				/* erase copied nodes */
				paddedTree.erase(it);
				
				it = iterN;
				currentNode++;
				depth_it++;
				
				paddedrmap[*it] = rtemp;
				paddedlmap[*it] = ltemp;
				
			}
	
		}
		
		it++;
	
	}
	
	/* construct paddedAncMap */
	paddedNodeCount = 0;
	paddedAncMap.clear();
	it = paddedTree.begin();
	end = paddedTree.end();
	while (it != end) {
		if (*it != 0) {
			paddedAncMap[*it] = *paddedTree.parent(it);
			paddedbmap[*it] = paddedtmap[*paddedTree.parent(it)] - paddedtmap[*it];
		}
		paddedNodeCount++;
		it++;
	}
	
	  	
}

/* reduces a tree to just its trunk, takes most recent sample and works backward from this */
void CoalescentTree::pruneToTrunk() {
	
	/* erase other nodes from the tree */
	tree<int>::iterator it, end;	
	it = ctree.begin();
	end = ctree.end();
	while(it!=end) {
		if (trunkSet.end() == trunkSet.find(*it)) {
			it = ctree.erase(it);
	//		tmap.erase(*it);			// cannot do this for some reason
	//		bmap.erase(*it);
	//		rmap.erase(*it);
		}
		else {
    		it++;
    	}
    }
        
    /* update other private data */
	/* maps don't need updating */
    nodeCount = ctree.size();
    tlist.clear();
    leafCount = 0;
	for (it=ctree.begin(); it!=ctree.end(); it++) {
		tlist.insert(tmap[*it]);
		if (ctree.number_of_children(it)==0)
			leafCount++;
	}
		
	constructPaddedTree(); 	
	
}

/* reduces a tree to samples with a single label */
void CoalescentTree::pruneToLabel(int label) {

	/* start by finding all tips with this label */
	set<int> labelset; 
	tree<int>::iterator it, jt, end;
	int newLeafs = 0;
	
	it = ctree.begin();
	end = ctree.end();
   
	if(!ctree.is_valid(it)) return;

	while(it!=end) {
		/* find leaf nodes at present */
		if (lmap[*it] == label) {     
			newLeafs++;
			jt = it;
			/* move up tree adding nodes to trunk set */
			while (*jt != 0) {
				labelset.insert(*jt);
				jt = ctree.parent(jt);
			}
		}
		it++;
	}
	
	/* make sure root is included */
	labelset.insert(0);
	
	/* erase other nodes from the tree */
	it = ctree.begin();
	end = ctree.end();
	while(it!=end) {
		if (labelset.end() == labelset.find(*it)) {
			it = ctree.erase(it);
	//		tmap.erase(*it);
	//		bmap.erase(*it);
	//		rmap.erase(*it);
		}
		else {
    		it++;
    	}
    }
        
    /* update other private data */
	/* maps don't need updating */
    leafCount = newLeafs;
    nodeCount = ctree.size();
    tlist.clear();
	for (it=ctree.begin(); it!=ctree.end(); it++)
		tlist.insert(tmap[*it]);
		
	constructPaddedTree(); 	
	
}

/* trims a tree at its edges */
/*
   		   |-------	 			 |-----
from  ------			 	to	--
   		   |----------		 	 |-----
*/
/* redo maps */
void CoalescentTree::trimEnds(double start, double stop) {
		
	/* erase nodes from the tree where neither the node nor its parent are between start and stop */
	tree<int>::iterator it, end, root, itertemp;	
	it = ctree.begin();
	end = ctree.end();
	root = ctree.begin();
	while(it!=end) {	
		/* if node < start and anc[node] > start, erase children and prune node back to start */
		if (tmap[*it] < start && tmap[ancMap[*it]] >= start) {
			tmap[*it] = start;
			bmap[*it] = tmap[ancMap[*it]] - tmap[*it];
			ctree.erase_children(it);
			it++;
		}
		/* if node < stop and anc[node] > stop, push anc[node] up to stop */
		/* and reparent anc[node] to be a child of root */
		/* neither node nore anc[node] can be root */
		else if (tmap[*it] < stop && tmap[ancMap[*it]] > stop && *it != 0 && ancMap[*it] != 0) {
			tmap[ancMap[*it]] = stop;
			bmap[ancMap[*it]] = 0;
			itertemp = ctree.append_child(root);
			ctree.move_ontop(itertemp, ctree.parent(it));
			it = ctree.begin();
		}
		else {
    		it++;
    	}
    }
    
    /* second pass for nodes < stop */
    it = ctree.begin();
	end = ctree.end();
	while(it!=end) {	
		if (tmap[*it] > stop && *it != 0) {
			it = ctree.erase(it);
		}
		else {
    		it++;
    	}
    }
    tmap[0] = stop;
        
/*
	int leafCount;						
	int nodeCount;						
	set<int> leafSet;					
	set<int> trunkSet;					
	set<int> labelSet;															
	map<int,double>	bmap;				
	map<int,double>	rmap;				
	map<int,double>	tmap;				
	map<int,int>	lmap;				
	map<int,int>	ancMap;	
*/

	/* want to reduce maps to just the subset dealing with these nodes */
	/* set of remaining nodes */
	set<int> nodeSet;
	for (it = ctree.begin(); it != ctree.end(); it++) {
		nodeSet.insert(*it);
	}
	/* go through and erase nonfunctionaling map elements */
	for (int i = 0; i < nodeCount; i++) {
		if (nodeSet.count(i) == 0) {
			leafSet.erase(i);
			trunkSet.erase(i);
			bmap.erase(i);
			rmap.erase(i);
			tmap.erase(i);
			lmap.erase(i);
			ancMap.erase(i);
		}
	}

    /* update other private data */
	/* maps don't need updating */
    nodeCount = ctree.size();
    tlist.clear();
    leafCount = 0;
	for (it=ctree.begin(); it!=ctree.end(); it++) {
		tlist.insert(tmap[*it]);
		if (ctree.number_of_children(it)==0)
			leafCount++;
	}
	
	/* this isn't working for some reason */
	constructPaddedTree(); 	
	
}


/* Print indented tree */
void CoalescentTree::printTree() { 

	tree<int>::pre_order_iterator it;
	tree<int>::pre_order_iterator end;

	it = ctree.begin();
	end = ctree.end();
   
	if(!ctree.is_valid(it)) return;
	
	int rootdepth=ctree.depth(it);
	cout << "-----" << endl;
	while(it!=end) {
		for(int i=0; i<ctree.depth(it)-rootdepth; ++i) 
			cout << "  ";
		cout << (*it) << endl << flush;
		it++;
	}
	cout << "-----" << endl;
}

/* Print parentheses tree */
void CoalescentTree::printParenTree() { 

	tree<int>::post_order_iterator it;
	tree<int>::post_order_iterator end;

	it = ctree.begin_post();
	end = ctree.end_post();
   
	if(!ctree.is_valid(it)) return;
	
	int rootdepth = ctree.depth(it);
	int currentdepth = ctree.depth(it);
	for (int i = 0; i < currentdepth; i++) { cout << "("; } 
	cout << (*it) << "";
	it++;
	
	/* basically, need to add a '(' whenever the depth increases and a ')' whenever the depth decreases */
	/* only print leaf nodes */
	while(it!=end) {
		if (ctree.depth(it) > currentdepth) { 
			cout << ", ("; 
			for (int i = 0; i < ctree.depth(it) - currentdepth - 1; i++) { cout << "("; }
			if (*it <= leafCount && *it > 0) { cout << (*it) << ""; }
		}
		if (ctree.depth(it) == currentdepth) { 
			if (*it <= leafCount && *it > 0) { cout << ", " << (*it) << ""; }
		}
		if (ctree.depth(it) < currentdepth) {
			if (*it <= leafCount && *it > 0) { cout << (*it) << ""; }
			cout << ")"; 
		}
		currentdepth = ctree.depth(it);
		it++;
		
	}
	
	cout << endl;

	
}

/* Print tree with times */
void CoalescentTree::printTimeTree() { 

	tree<int>::pre_order_iterator it;
	tree<int>::pre_order_iterator end;

	it = ctree.begin();
	end = ctree.end();
   
	if(!ctree.is_valid(it)) return;
		
	int rootdepth=ctree.depth(it);
	cout << "-----" << endl;
	while(it!=end) {
		for(int i=0; i<ctree.depth(it)-rootdepth; ++i) 
			cout << "  ";
		cout << (*it) << " (" << tmap[*it] << ")";
		if (brCheck)
			cout << " [" << rmap[*it] << "]";
		if (nlCheck)
			cout << " [" << lmap[*it] << "]";			
		cout << endl << flush;
		it++;
	}
	cout << "-----" << endl;
	cout << "leafs = " << leafCount << ", nodes = " << nodeCount << endl;
	
}

/* print tree in Mathematica suitable format, 1->2,1->3, etc... */
/* padded with extra nodes at coalescent time points */
/* print mapping of nodes to rates if rates exist */
void CoalescentTree::printPaddedRuleList() { 

	/* print tip count */
	cout << leafCount << endl;
	
	/* print nodes that exist at present */
	for (int i=1; i<leafCount; i++) {
		if (paddedtmap[i] < 0.01)				// should specify this better
			cout << i << " ";
	}
	cout << endl;
	
	/* print trunk nodes */
	
	set<int> trunk; 
	
	tree<int>::iterator it, jt, end;
	it = paddedTree.begin();
	end = paddedTree.end();
   
	while(it!=end) {
		/* find leaf nodes at present */
		if (paddedtmap[*it] < 0.01) {
			jt = it;
			/* move up tree adding nodes to trunk set */
			while (*jt != 0) {
				trunk.insert(*jt);
				jt = paddedTree.parent(jt);
			}
		}
		it++;
	}
	
	for ( set<int>::const_iterator lit=trunk.begin(); lit != trunk.end(); lit++ ) {
		cout << *lit << " ";
	}
	cout << endl;
	
	/* print the tree in rule list (Mathematica-ready) format */
	/* print only upward links */
	tree<int>::pre_order_iterator iterTemp, iterN;
	it = paddedTree.begin();
	end = paddedTree.end();
	
	it++;
	cout << *it << "->" << *paddedTree.parent(it);
	it++;
	while(it!=end) {
		cout << " " << *it << "->" << *paddedTree.parent(it);
		it++;
	}
	cout << endl;
	
	/* print list of coalescent times in Mathematica-ready format */
	cout << *(tlist.begin());
  	for( set<double>::const_iterator iter = ++tlist.begin(); iter != tlist.end(); iter++ ) {
		cout << " " << *iter;
    }
	cout << endl;
	
	/* print mapping of rates in Mathematica format */
	if (brCheck) {
		it = paddedTree.begin();
		it++;
		cout << *it << "->" << paddedrmap[*it];
		it++;
		while(it!=end) {
			cout << " " << *it << "->" << paddedrmap[*it];
			it++;
		}
		cout << endl;
	}
	
	/* print mapping of labels in Mathematica format */
	if (nlCheck) {
		it = paddedTree.begin();
		it++;
		cout << *it << "->" << paddedlmap[*it];
		it++;
		while(it!=end) {
			cout << " " << *it << "->" << paddedlmap[*it];
			it++;
		}
		cout << endl;
	}	
	  	
  	
}


/* print tree in Mathematica suitable format, 1->2,1->3, etc... */
/* padded with extra nodes at coalescent time points */
/* print mapping of nodes to rates if rates exist */
void CoalescentTree::printRuleList() { 

	tree<int>::iterator it, jt, end;

	/* reorder tree so that the bottom node of two sister nodes always has more children */
	it = ctree.begin();
	end = ctree.end();	
	while(it!=end) {
		jt = ctree.next_sibling(it);
		if (ctree.is_valid(jt)) {
			int cit = ctree.size(it);
			int cjt = ctree.size(jt);
			if (cit > cjt) {
				ctree.swap(jt,it);
				it = ctree.begin();
			}
		}
		it++;
	}

	/* print tip count */
	cout << leafCount << endl;
		
	/* print trunk nodes */
	for ( set<int>::const_iterator lit=trunkSet.begin(); lit != trunkSet.end(); lit++ ) {
		cout << *lit << " ";
	}
	cout << endl;
	
	/* print the tree in rule list (Mathematica-ready) format */
	/* print only upward links */
	tree<int>::pre_order_iterator iterTemp, iterN;
	it = ctree.begin();
	end = ctree.end();
	
	it++;
	cout << *it << "->" << *ctree.parent(it);
	it++;
	while(it!=end) {
		cout << " " << *it << "->" << *ctree.parent(it);
		it++;
	}
	cout << endl;
	
	
	/* print mapping of labels in Mathematica format */
	if (nlCheck) {
		it = ctree.begin();
		it++;
		cout << *it << "->" << lmap[*it];
		it++;
		while(it!=end) {
			cout << " " << *it << "->" << lmap[*it];
			it++;
		}
		cout << endl;
	}	
	
	/* print mapping of nodes to coordinates */
  	it = ctree.begin();
  	int count = 0;
  	cout << *it << "->{" << -1 * tmap[*it] << "," << count << "}";
	it++;
	count++;
	while(it!=end) {
		cout << " " << *it << "->{" << -1 * tmap[*it] << "," << count << "}";
		it++;
		count++;
	}
	cout << endl;	
	  	  	
}

/* print proportion of tree that can trace its history forward to present day samples */
void CoalescentTree::printTrunkRatio() { 

	double totalLength = 0.0;
	double trunkLength = 0.0;
	
	set<int> trunk; 
	tree<int>::iterator it, jt, end;

	/* sum branch lengths */
	for (int i=0; i<nodeCount; i++) {
		totalLength += bmap[i];
	}
	cout << "total tree length: " << totalLength << endl;
	
	/* print nodes that exist at present */
	cout << "leaf nodes at present:" << endl;
	
	it = ctree.begin();
	end = ctree.end();
   
	if(!ctree.is_valid(it)) return;

	while(it!=end) {
		/* find leaf nodes at present */
		if (tmap[*it] < 1 && *it <= leafCount) {
			jt = it;
			/* move up tree adding nodes to trunk set */
			while (*jt != 0) {
				trunk.insert(*jt);
				jt = ctree.parent(jt);
			}
		}
		it++;
	}
	
	/* move through trunk nodes and total branch lengths */
	for ( set<int>::const_iterator t=trunk.begin(); t != trunk.end(); t++ ) {
		trunkLength += bmap[*t];
    }
    
    cout << "trunk length: " << trunkLength << endl;
    cout << "ratio: " << trunkLength / totalLength << endl;
    	
}

/* print proportion of tree with each label */
void CoalescentTree::printLabelProportions() { 

	/* calculation is to get the total length in years of each label */

	double totalLength = 0;
	map <int,double> labelLengths;
	set<int>::const_iterator lit, ljt;

	for ( lit=labelSet.begin(); lit != labelSet.end(); lit++ ) {
   
		/* sum branch lengths */
		/* have to go through tree structure, nodes can be non-consecutive */
		double length = 0.0;
		for (tree<int>::iterator it = ctree.begin(); it != ctree.end(); it++ ) {
			if (*lit == lmap[*it]) {
				length += bmap[*it];
			}
		}
		totalLength += length;
//		cout << "label: " << *lit << ", length: " << length << endl;
		labelLengths.insert( make_pair(*lit,length ) );
		
	}
//	cout << "total length: " << totalLength << endl;
	
	for ( lit=labelSet.begin(); lit != labelSet.end(); lit++ ) {
		cout << labelLengths[*lit] / totalLength << " ";
	}
	cout << endl;
	
}

/* print total rate of migration across tree, measured as events/year */
void CoalescentTree::printMigTotal() { 

	/* calculation is to get the total length in years of each label */
	double totalLength = 0;
	set<int>::const_iterator lit, ljt;
	for ( lit=labelSet.begin(); lit != labelSet.end(); lit++ ) {
   
		/* sum branch lengths */
		/* have to go through tree structure, nodes can be non-consecutive */
		for (tree<int>::iterator it = ctree.begin(); it != ctree.end(); it++ ) {
			if (*lit == lmap[*it]) {
				totalLength += bmap[*it];
			}
		}		
	}
//	cout << "total length: " << totalLength << endl;
	
	/* count migration events */	
	int count = 0;
	for ( ljt=labelSet.begin(); ljt != labelSet.end(); ljt++ ) {		// to
		for ( lit=labelSet.begin(); lit != labelSet.end(); lit++ ) {	// from
			if (*lit != *ljt) {
				
				/* go through tree and find situations where parent matches lit and child matches ljt */
				/* using ancMap for speed */
				for (tree<int>::iterator it = ctree.begin(); it != ctree.end(); it++ ) {	
					if ( lmap[*it] == *ljt && lmap[ancMap[*it]] == *lit 
					&& *it != 0 && ancMap[*it] != 0) {						// exclude root
						count++;
//						cout << "from " << i << " [" << *lit << "]" << " to " << paddedAncMap[i] << " [" << *ljt << "]" << endl;
					}
				}
				
			}
		}
	}
	
	/* rate is equal to events / opportunity */
	cout << count / totalLength << endl;
	
}


/* print matrix of migration rates */
void CoalescentTree::printMigRates() { 

	/* calculation is to get the total length in years of each label */

	double totalLength = 0;
	map <int,double> labelLengths;
	set<int>::const_iterator lit, ljt;

	for ( lit=labelSet.begin(); lit != labelSet.end(); lit++ ) {
   
		/* sum branch lengths */
		/* have to go through tree structure, nodes can be non-consecutive */
		double length = 0.0;
		for (tree<int>::iterator it = ctree.begin(); it != ctree.end(); it++ ) {
			if (*lit == lmap[*it]) {
//				cout << *it << " " << bmap[*it] << endl;
				length += bmap[*it];
			}
		}
		totalLength += length;
//		cout << "label: " << *lit << ", length: " << length << endl;
		labelLengths.insert( make_pair(*lit,length ) );
		
	}
//	cout << "total length: " << totalLength << endl;
	
	/* count migration events */
	/* ordering matches that used by Migrate */
	
	for ( ljt=labelSet.begin(); ljt != labelSet.end(); ljt++ ) {		// to
		for ( lit=labelSet.begin(); lit != labelSet.end(); lit++ ) {	// from
			if (*lit != *ljt) {
				
				/* go through tree and find situations where parent matches lit and child matches ljt */
				/* using ancMap for speed */
				int count = 0;
				for (tree<int>::iterator it = ctree.begin(); it != ctree.end(); it++ ) {	
					if ( lmap[*it] == *ljt && lmap[ancMap[*it]] == *lit 
					&& *it != 0 && ancMap[*it] != 0) {						// exclude root
						count++;
//						cout << "from " << i << " [" << *lit << "]" << " to " << paddedAncMap[i] << " [" << *ljt << "]" << endl;
					}
				}
		
				/* rate is equal to events / opportunity */
				/* divide by length of lit */
				
//				cout << "from " << *lit << " to " << *ljt << ": " << "count = " << count << ", rate = ";
				cout << count / labelLengths[*lit] << " ";
//				cout << endl;
		
			}
		}
	}
	cout << endl;
	
}


map<int,double> CoalescentTree::getMigWeights() { 

	/* calculation is to get the total length in years of each label */

	double totalLength = 0;
	map <int,double> labelLengths;
	set<int>::const_iterator lit, ljt;
	
	map<int,double> migWeights;

	for ( lit=labelSet.begin(); lit != labelSet.end(); lit++ ) {
   
		/* sum branch lengths */
		/* have to go through tree structure, nodes can be non-consecutive */
		double length = 0.0;
		for (tree<int>::iterator it = ctree.begin(); it != ctree.end(); it++ ) {
			if (*lit == lmap[*it]) {
				if (bmap[*it] > 0) {			// correcting for the occasional negative branch length
					length += bmap[*it];
				}
			}
		}
		totalLength += length;
//		cout << "label: " << *lit << ", length: " << length << endl;
		labelLengths.insert( make_pair(*lit,length ) );
		
	}
//	cout << "total length: " << totalLength << endl;
		
	for ( ljt=labelSet.begin(); ljt != labelSet.end(); ljt++ ) {		// to
		for ( lit=labelSet.begin(); lit != labelSet.end(); lit++ ) {	// from
			if (*lit != *ljt) {
								
				migWeights[*ljt * 10 + *lit] = labelLengths[*lit];
		
			}
		}
	}
	
	return migWeights;
	
}

map<int,int> CoalescentTree::getMigCounts() { 

	/* calculation is to get the total length in years of each label */

	set<int>::const_iterator lit, ljt;	
	map<int,int> migCounts;
		
	for ( ljt=labelSet.begin(); ljt != labelSet.end(); ljt++ ) {		// to
		for ( lit=labelSet.begin(); lit != labelSet.end(); lit++ ) {	// from
			if (*lit != *ljt) {
				
				/* go through tree and find situations where parent matches lit and child matches ljt */
				/* using ancMap for speed */
				int count = 0;
				for (tree<int>::iterator it = ctree.begin(); it != ctree.end(); it++ ) {	
					if ( lmap[*it] == *ljt && lmap[ancMap[*it]] == *lit 
					&& *it != 0 && ancMap[*it] != 0					// exclude root
					&& bmap[*it] > 0 && bmap[ancMap[*it]] > 0) {	// correcting for the occasional negative branch length			
						count++;
//						cout << "from " << i << " [" << *lit << "]" << " to " << paddedAncMap[i] << " [" << *ljt << "]" << endl;
					}
				}
		
				/* rate is equal to events / opportunity */
				/* divide by length of lit */
				
				migCounts[*ljt * 10 + *lit] = count;
		
			}
		}
	}
	
	return migCounts;
	
}

map<int,double> CoalescentTree::getRevMigWeights() { 

	/* calculation is to get the total length in years of each label */

	double totalLength = 0;
	map <int,double> labelLengths;
	set<int>::const_iterator lit, ljt;
	
	map<int,double> migWeights;

	for ( lit=labelSet.begin(); lit != labelSet.end(); lit++ ) {
   
		/* sum branch lengths */
		/* have to go through tree structure, nodes can be non-consecutive */
		double length = 0.0;
		for (tree<int>::iterator it = ctree.begin(); it != ctree.end(); it++ ) {
			if (*lit == lmap[*it]) {
				if (bmap[*it] > 0) {			// correcting for the occasional negative branch length
					length += bmap[*it];
				}
			}
		}
		totalLength += length;
//		cout << "label: " << *lit << ", length: " << length << endl;
		labelLengths.insert( make_pair(*lit,length ) );
		
	}
//	cout << "total length: " << totalLength << endl;
		
	for ( ljt=labelSet.begin(); ljt != labelSet.end(); ljt++ ) {		// to
		for ( lit=labelSet.begin(); lit != labelSet.end(); lit++ ) {	// from
			if (*lit != *ljt) {
								
				migWeights[*ljt * 10 + *lit] = labelLengths[*ljt];		// forward uses from deme
																		// backward uses to deme
		
			}
		}
	}
	
	return migWeights;
	
}

map<int,int> CoalescentTree::getRevMigCounts() { 

	set<int>::const_iterator lit, ljt;	
	map<int,int> migCounts;
		
	for ( ljt=labelSet.begin(); ljt != labelSet.end(); ljt++ ) {		// to
		for ( lit=labelSet.begin(); lit != labelSet.end(); lit++ ) {	// from
			if (*lit != *ljt) {
				
				/* go through tree and find situations where parent matches lit and child matches ljt */
				/* using ancMap for speed */
				int count = 0;
				for (tree<int>::iterator it = ctree.begin(); it != ctree.end(); it++ ) {	
					if ( lmap[*it] == *lit	 						// from deme is child
					&& lmap[ancMap[*it]] == *ljt 					// to deme is parent
					&& *it != 0 && ancMap[*it] != 0					// exclude root
					&& bmap[*it] > 0 && bmap[ancMap[*it]] > 0) {	// correcting for the occasional negative branch length			
						count++;
//						cout << "from " << i << " [" << *lit << "]" << " to " << paddedAncMap[i] << " [" << *ljt << "]" << endl;
					}
				}
		
				/* rate is equal to events / opportunity */
				/* divide by length of lit */
				
				migCounts[*ljt * 10 + *lit] = count;
		
			}
		}
	}
	
	return migCounts;
	
}


/* go through tree and print the coalescent rate for each label */
void CoalescentTree::printCoalRates() {

	double step = 0.0001;
	tree<int>::iterator it;
	
	/* go through tree and find end points */
	double start, stop;
	start = tmap[*ctree.begin()];
	stop = tmap[*ctree.begin()];
	for (it = ctree.begin(); it != ctree.end(); it++) {
		if (start > tmap[*it]) {
			start = tmap[*it];
		}
	}
	
//	cout << "start = " << start << ", stop = " << stop << endl;
	
	set<int>::const_iterator lit;
	for ( lit=labelSet.begin(); lit != labelSet.end(); lit++ ) {
	
//		cout << "label = " << *lit << endl;
	
		/* go through tree and count concurrent lineages */
		/* weight number by n(n-1)/2 to get coalescent opportunity */
		double labelWeight = 0.0;
		for (double t = start; t <= stop; t += step) {
		
//			cout << t << endl;
			int count = 0;
			for (it = ctree.begin(); it != ctree.end(); it++) {
				if (tmap[*it] < t && tmap[ancMap[*it]] >= t && lmap[*it] == *lit) {		
//					cout << *it << " ";
					count++;
				}
			}
			double weight;
			if (count > 0)
				weight = ( ( count * (count - 1) ) / 2 ) * step;
			else
				weight = 0;
//			cout << "count = " << count << ", weight = " << weight << endl;
			labelWeight += weight;
		
		}
		
		
		/* count coalescent events, these are nodes with two children */
		int count = 0;
		for (it = ctree.begin(); it != ctree.end(); it++) {
			if (ctree.number_of_children(it) == 2 && lmap[*it] == *lit) {		
				count++;
			}
		}
//		cout << "label = " << *lit << ", weight = " << labelWeight << ", count = " << count << ", rate = " << count / labelWeight << ", timescale = " << labelWeight / count << endl;
		cout << count / labelWeight << " ";

	}
	cout << endl;
	
}

/* return map of coalescent weights for each label */
map<int,double> CoalescentTree::getCoalWeights() {

	double step = 0.0001;
	tree<int>::iterator it;
	map<int,double> labelWeights;
	
	/* go through tree and find end points */
	double start, stop;
	start = tmap[*ctree.begin()];
	stop = tmap[*ctree.begin()];
	for (it = ctree.begin(); it != ctree.end(); it++) {
		if (start > tmap[*it]) {
			start = tmap[*it];
		}
	}
		
	for ( set<int>::const_iterator lit = labelSet.begin(); lit != labelSet.end(); lit++ ) {
		
		labelWeights.insert(make_pair(*lit,0.0));
		
		/* go through tree and count concurrent lineages */
		/* weight number by n(n-1)/2 to get coalescent opportunity */
		for (double t = start; t <= stop; t += step) {
		
//			cout << t << endl;
			int count = 0;
			for (it = ctree.begin(); it != ctree.end(); it++) {
				if (tmap[*it] < t && tmap[ancMap[*it]] >= t && lmap[*it] == *lit) {		
//					cout << *it << " ";
					count++;
				}
			}
			double weight;
			if (count > 0)
				weight = ( ( count * (count - 1) ) / 2 ) * step;
			else
				weight = 0;
//			cout << "count = " << count << ", weight = " << weight << endl;
			labelWeights[*lit] += weight;
		
		}
			
	}
	
	return labelWeights;
}

/* return map of coalescent counts for each label */
/* making it a double to talk better with LabelSummary */
map<int,int> CoalescentTree::getCoalCounts() {

	tree<int>::iterator it;
	map<int,int> labelCounts;
	
	/* go through tree and find end points */
	double start, stop;
	start = tmap[*ctree.begin()];
	stop = tmap[*ctree.begin()];
	for (it = ctree.begin(); it != ctree.end(); it++) {
		if (start > tmap[*it]) {
			start = tmap[*it];
		}
	}
		
	for ( set<int>::const_iterator lit = labelSet.begin(); lit != labelSet.end(); lit++ ) {
		
		labelCounts.insert(make_pair(*lit,0));
				
		/* count coalescent events, these are nodes with two children */
		for (it = ctree.begin(); it != ctree.end(); it++) {
			if (ctree.number_of_children(it) == 2 && lmap[*it] == *lit) {		
				labelCounts[*lit]++;
			}
		}
			
	}
	
	return labelCounts;
}



set<int> CoalescentTree::getLabelSet() {
	return labelSet;
}

/* Bayesian skyline for effective coalesent timescale Ne*tau */
void CoalescentTree::NeSkyline() { 

	vector<int> lineages (tlist.size());
	tree<int>::iterator it, end;
	it = paddedTree.begin();
	end = paddedTree.end(); 		
	int rootdepth=paddedTree.depth(it);
	
	while(it!=end) {
		lineages.at(paddedTree.depth(it)-rootdepth)++;
		it++;
	}
	
	int loc = paddedTree.max_depth();
	double ne;
	int count;
	double step = stepsize;
  	for( set<double>::const_iterator iter = tlist.begin(); iter != --tlist.end(); iter++ ) {
  		double a = *iter;
  		iter++;
  		double b = *iter;
  		iter--;
  		if (a < b - 0.00000001) {
  			if (lineages.at(loc-1) < lineages.at(loc))
				ne = (b - a) * 0.5 * lineages.at(loc) * (lineages.at(loc) - 1);
	//		cout << "{" << a << "," << b << "} " << lineages.at(loc) << " " << ne << endl;
			while (step > a && step < b) {
	//			cout << step << "\t" << ne << endl;
				skylineindex.push_back(step);
				skylinevalue.push_back(ne);
				step += stepsize;
			}
	//		cout << step << "\t" << ne << endl;
			skylineindex.push_back(step);
			skylinevalue.push_back(ne);
			step += stepsize;
		}
		loc--;
    }
    
}

/* Skyline for rate of substitution, take mean rate for concurrent lineages */
void CoalescentTree::subRateSkyline() { 

	tree<int>::iterator it, jt, end;
	vector<int> lineages (tlist.size());
	vector<double> rates (tlist.size());	
	it = paddedTree.end();
	end = paddedTree.end(); 		
	int rootdepth=paddedTree.depth(it);
	
	/* lineages contains a running tally of mean rate at this depth */
	while(it!=paddedTree.begin()) {
		lineages.at(paddedTree.depth(it)-rootdepth)++;
		rates.at(paddedTree.depth(it)-rootdepth) += paddedrmap[*it];
//		cout << *it << " " << paddedrmap[*it] << endl;
		it--;
	}
	
	int loc = paddedTree.max_depth();
	double rate;
	int count;
	double step = stepsize;
  	for( set<double>::const_iterator iter = tlist.begin(); iter != --tlist.end(); iter++ ) {
  		double a = *iter;
  		iter++;
  		double b = *iter;
  		iter--;
  		if (a < b - 0.00000001) {
			rate = rates.at(loc) / lineages.at(loc);
//			cout << "{" << a << "," << b << "} " << lineages.at(loc) << " " << rate << endl;
			while (step > a && step < b) {
//				cout << step << "\t" << rate << endl;
				skylineindex.push_back(step);
				skylinevalue.push_back(rate);
				step += stepsize;
			}
//			cout << step << "\t" << rate << endl;
			skylineindex.push_back(step);
			skylinevalue.push_back(rate);
			step += stepsize;
		}
		loc--;
    }
    
}


/* Skyline for diversity (pi) */
/* At every event, take concurrent lineages and make all pairwise comparisons */
/* calculating total distance separating samples, E[div] = 2 Ne*tau */
void CoalescentTree::divSkyline() { 

	/* have a vector[time points] of sets[node labels] */
	tree<int>::iterator it, jt, end;
	vector< set<int> > lineages (tlist.size());
	it = paddedTree.end();
	end = paddedTree.end(); 		
	int rootdepth=paddedTree.depth(it);
	
	/* go through tree and add to set */
	while(it!=paddedTree.begin()) {
		(lineages.at(paddedTree.depth(it)-rootdepth)).insert(*it);
		it--;
	}
	
	int loc = paddedTree.max_depth();
	double rate;
	int count;
	double step = stepsize;
  	for( set<double>::const_iterator iter = tlist.begin(); iter != --tlist.end(); iter++ ) {
  		set<int> nodes = lineages.at(loc);
  		double a = *iter;
  		iter++;
  		double b = *iter;
  		iter--;
  		if (a < b - 0.00000001) {
//			cout << "{" << a << "," << b << "} ";
			int divn = 0;
			double div = 0.0;
			for( set<int>::const_iterator jter = nodes.begin(); jter != --nodes.end(); jter++ ) {
				for( set<int>::const_iterator kter = jter; kter != nodes.end(); kter++ ) {
					if (*jter < *kter) { 
						set<int> s;
						s.insert(*jter);
						s.insert(*kter);
						tree<int> st = extractSubtree(s);
						divn++;
						div += getTreeLength(st);
					}
				}
			}
			div /= divn;
//			cout << div << endl;
			while (step > a && step < b) {
				skylineindex.push_back(step);
				skylinevalue.push_back(div);
				step += stepsize;
			}
			skylineindex.push_back(step);
			skylinevalue.push_back(div);
			step += stepsize;
		}
		loc--;
    }
   
}

/* Skyline for segregating site (S/a1) */
/* At every event, take concurrent lineages and calculate total tree length */
/* a1 = sum 1 through n-1 of 1/i */
void CoalescentTree::tajimaSkyline() { 

	/* have a vector[time points] of sets[node labels] */
	tree<int>::iterator it, jt, end;
	vector< set<int> > lineages (tlist.size());
	it = paddedTree.end();
	end = paddedTree.end(); 		
	int rootdepth=paddedTree.depth(it);
	
	/* go through tree and add to set */
	while(it!=paddedTree.begin()) {
		(lineages.at(paddedTree.depth(it)-rootdepth)).insert(*it);
		it--;
	}
	
	int loc = paddedTree.max_depth();
	double rate;
	int count;
	double step = stepsize;
  	for( set<double>::const_iterator iter = tlist.begin(); iter != --tlist.end(); iter++ ) {
  		set<int> nodes = lineages.at(loc);
  		double a = *iter;
  		iter++;
  		double b = *iter;
  		iter--;
  		if (a < b - 0.00000001) {
	//		cout << "{" << a << "," << b << "} ";
			int divn = 0;
			double div = 0.0;
			for( set<int>::const_iterator jter = nodes.begin(); jter != --nodes.end(); jter++ ) {
				for( set<int>::const_iterator kter = jter; kter != nodes.end(); kter++ ) {
					if (*jter < *kter) { 
						set<int> s;
						s.insert(*jter);
						s.insert(*kter);
						tree<int> st = extractSubtree(s);
						divn++;
						div += getTreeLength(st);
					}
				}
			}
			div /= divn;
			
			double treeS;
			double a1 = 0.0;
			double a2 = 0.0;
			int n = nodes.size();
			tree<int> st = extractSubtree(nodes);
			treeS = getTreeLength(st);
			for (int i = 1; i < n; i++) {
				a1 += 1.0/i;
				a2 += 1.0/(i*i);
			}
	//		cout << " " << treeS << " " << a1 << " " << treeS/a1 << endl;
			
			double e1 = (1.0/a1) * ((double)(n+1) / (3*(n-1)) - (1.0/a1));
			double e2 = (1.0 / (a1*a1 + a2) ) * ( (double)(2*(n*n+n+3)) / (9*n*(n-1)) - (double)(n+2) / (n*a1) + a2/(a1*a1) );
			double denom = sqrt(e1*treeS + e2*treeS*(treeS-1));
			double tajima = (div - treeS / a1) / denom;
			
	//		double tajima = div - treeS / a1;
			
	//		cout << "a1: " << a1 << " a2: " << a2 << " e1: " << e1 << " e2: " << e2 << " denom: " << denom << " D: " << tajima << endl;
			
			while (step > a && step < b) {
				skylineindex.push_back(step);
				skylinevalue.push_back(tajima);
				step += stepsize;
			}
			skylineindex.push_back(step);
			skylinevalue.push_back(tajima);
			step += stepsize;
		}
		loc--;
    }
   
}

/* Skyline with time for each sampled lineage to coalesce with phylogeny trunk */
/* only looks at samples, not full tree */
void CoalescentTree::tcSkyline() { 

	tree<int>::iterator it, jt, end;
	it = ctree.begin();
	end = ctree.end(); 

	int steps = ceil(tmap[0] / stepsize);
	vector<int> counts (steps);
	vector<double> tc (steps);
	
	for (double i = 0.0; i <= tmap[0]; i+=stepsize) {
		skylineindex.push_back(i);
	}
	
	while(it!=ctree.end()) {
		if (leafSet.end() != leafSet.find(*it)) {
			
			jt = it;
			while (trunkSet.end() == trunkSet.find(*jt)) {
				jt = ctree.parent(jt);
			}
	//		cout << tmap[*it] << " " << tmap[*jt] - tmap[*it] << endl;
			counts.at( floor(tmap[*it] / stepsize + 0.5) )++;
			tc.at( floor(tmap[*it] / stepsize + 0.5) ) += tmap[*jt] - tmap[*it];
			
		}
		it++;
	}
    
    for (int i = 0; i < steps; i++) {
    	skylinevalue.push_back( tc[i] / counts[i] );
	}
		
	
}

/* Skyline for proportion with a particular label, take mean proportion for concurrent lineages */
void CoalescentTree::labelSkyline(int label) { 

	tree<int>::iterator it, jt, end;
	vector<int> lineages (tlist.size());
	vector<double> proportions (tlist.size());	
	it = paddedTree.end();
	end = paddedTree.end(); 		
	int rootdepth=paddedTree.depth(it);
	
	/* lineages contains a running tally of mean rate at this depth */
	while(it!=paddedTree.begin()) {
		lineages.at(paddedTree.depth(it)-rootdepth)++;
		if (paddedlmap[*it] == label)
			proportions.at(paddedTree.depth(it)-rootdepth) += 1;
//		cout << *it << " " << paddedrmap[*it] << endl;
		it--;
	}
	
	int loc = paddedTree.max_depth();
	double proportion;
	int count;
	double step = stepsize;
	
  	for( set<double>::const_iterator iter = tlist.begin(); iter != --tlist.end(); iter++ ) {
  		double a = *iter;
  		iter++;
  		double b = *iter;
  		iter--;
  		if (a < b - 0.00000001) {
			proportion = proportions.at(loc) / lineages.at(loc);
//			cout << "loc = " << loc << endl;
//			cout << "{" << a << "," << b << "} " << lineages.at(loc) << " " << proportion << endl;
			while (step > a && step < b) {
//				cout << step << "\t" << proportion << endl;
				skylineindex.push_back(step);
				skylinevalue.push_back(proportion);
				step += stepsize;
			}
//			cout << step << "\t" << proportion << endl;
			skylineindex.push_back(step);
			skylinevalue.push_back(proportion);
			step += stepsize;
		}
		if (loc == 0)
			break;
		else
			loc--;
    }
    
}

vector<double> CoalescentTree::getSkylineIndex() {   
    return skylineindex;
}

vector<double> CoalescentTree::getSkylineValue() {   
    return skylinevalue;
}

void CoalescentTree::setStepSize(double ss) {
	stepsize = ss;
}

tree<int> CoalescentTree::extractSubtree(set<int> subset) {

	tree<int> stree;
	tree<int>::pre_order_iterator it, jt, end;

	/* fill top of tree with subset */
//	it = stree.set_head(*subset.begin());
	it = stree.set_head(0);
	for( set<int>::const_iterator jter = subset.begin(); jter != subset.end(); jter++ ) {
//		it = stree.insert_after(it, *jter);
		stree.append_child(it,*jter);
	}
	
	map<int,int> anc;
	
	while (stree.number_of_children(stree.begin()) > 1) {
			
		/* update stree based upon this mapping */
		for (tree<int>::sibling_iterator st = stree.begin(stree.begin()); st != stree.end(stree.begin()); st++) {
			st = stree.wrap(st,paddedAncMap[*st]);
		}	
		
		/* compare all pairs of nodes at top level of stree, if match, merge children of these nodes */	
		tree<int>::sibling_iterator st, tt;
		st = stree.begin(stree.begin());
		while (st != --stree.end(stree.begin())) {
			tt = st;
			tt++;
			while (tt != stree.end(stree.begin())) {
				if (*st == *tt) {
					stree.reparent(st,tt);		// moves children of tt to be children of st
				}
				tt++;
			}
			st++;
		}
		
		/* remove nodes in stree without any children */
		for (tree<int>::sibling_iterator st = stree.begin(stree.begin()); st != stree.end(stree.begin()); st++) {
			if (stree.number_of_children(st) == 0) {
				st = stree.erase(st);
			}
		}	
		
		/* new subset is the top level of tree */
		subset.clear();
		for (tree<int>::sibling_iterator st = stree.begin(stree.begin()); st != stree.end(stree.begin()); st++) {
			subset.insert(*st);
		}	
	
	}
	
	/* remove 0 from head of tree */
	stree.move_after(stree.begin(),++stree.begin());
	stree.erase(stree.begin());
			  
	return stree;	

}

double CoalescentTree::getTreeLength(tree<int> &tr) {

	double length = 0.0;

	tree<int>::pre_order_iterator it, end;	
	it = tr.begin();
	end = tr.end();
	while(it!=end) {
		if (tr.depth(it) > 0 )
			length += paddedtmap[*tr.parent(it)] - paddedtmap[*it];
		it++;
	}
	
	return length;

}

