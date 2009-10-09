/* node.h
Copyright 2009 Trevor Bedford <bedfordt@umich.edu>
Node class definition
This object represents a node in the CoalescentTree object, stores relevant attributes
*/

#ifndef NODE_H
#define NODE_H

#include <string>
using std::string;

class Node {

public:
	Node();							// defaults to -1
	Node(int);						// set number on creation
	
	bool operator==(Node &other);
	
	// GET FUNCTIONS
	int getNumber();
	string getName();
	double getLength();
	double getTime();
	int getLabel();
	double getCoord();
	bool getLeaf();
	bool getTrunk();
	bool getInclude();

	// SET FUNCTIONS
	void setNumber(int);
	void setName(string);
	void setLength(double);
	void setTime(double);
	void setLabel(int);
	void setCoord(double);
	void setLeaf(bool);
	void setTrunk(bool);
	void setInclude(bool);
																			
private:
	int number;						// number of node, must be unique
	string name;					// name of node, doesn't have to exist
	double length;					// length of the branch leading into the node
	double time;					// date of the node	
	int label;						// numeric label associated with node
	double coord;					// y-axis coordinate, used for tree drawing
	bool leaf;						// is this node a leaf?
	bool trunk;						// is this node part of the trunk?
	bool include;					// include this node when performing calculations?


};

#endif