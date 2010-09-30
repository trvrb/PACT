/* node.h
Copyright 2009-2010 Trevor Bedford <bedfordt@umich.edu>
Node class definition
This object represents a node in the CoalescentTree object, stores relevant attributes
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
	string getLabel();
	double getX();
	double getY();	
	double getCoord();
	bool getLeaf();
	bool getTrunk();
	bool getInclude();

	// SET FUNCTIONS
	void setNumber(int);
	void setName(string);
	void setLength(double);
	void setTime(double);
	void setLabel(string);
	void setX(double);
	void setY(double);	
	void setCoord(double);
	void setLeaf(bool);
	void setTrunk(bool);
	void setInclude(bool);
																			
private:
	int number;						// number of node, must be unique
	string name;					// name of node, doesn't have to exist
	double length;					// length of the branch leading into the node
	double time;					// date of the node	
	string label;					// arbitrary label associated with node
	double xloc;					// x-axis location of the node
	double yloc;					// y-axis location of the node	
	double coord;					// y-axis coordinate, used for tree drawing
	bool leaf;						// is this node a leaf?
	bool trunk;						// is this node part of the trunk?
	bool include;					// include this node when performing calculations?


};

#endif