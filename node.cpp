/* node.cpp
Copyright 2009 Trevor Bedford <bedfordt@umich.edu>
Member function definitions for Node class
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

#include <string>
using std::string;

#include "node.h"

Node::Node() {
	
	number = -1;
	name = "";
	length = 0.0;
	time = 0.0;
	label = "1";
	coord = 0.0;
	leaf = false;
	trunk = false;
	include = true;
		
}

Node::Node(int n) {
	
	number = n;
	name = "";
	length = 0.0;
	time = 0.0;
	label = "1";
	coord = 0.0;
	leaf = false;
	trunk = false;
	include = true;
		
}

// overloaded comparison operator
bool Node::operator==(Node &other) {

	int i = (*this).getNumber();
	int j = other.getNumber();
	return i == j;

}

/* Get functions */
int Node::getNumber() { return number; }
string Node::getName() { return name; }
double Node::getLength() { return length; }
double Node::getTime() { return time; }
string Node::getLabel() { return label; }
double Node::getCoord() { return coord; }
bool Node::getLeaf() { return leaf; }
bool Node::getTrunk() { return trunk; }
bool Node::getInclude() { return include; }

/* Set functions */
void Node::setNumber(int n) { number = n; }
void Node::setName(string n) { name = n; }
void Node::setLength(double n) { length = n; }
void Node::setTime(double n) { time = n; }
void Node::setLabel(string n) { label = n; }
void Node::setCoord(double n) { coord = n; }
void Node::setLeaf(bool n) { leaf = n; }
void Node::setTrunk(bool n) { trunk = n; }
void Node::setInclude(bool n) { include = n; }
