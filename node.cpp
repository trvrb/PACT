/* node.cpp
Member function definitions for Node class
*/

#include "node.h"

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

Node::Node(int n) {
	
	number = n;
	name = "";
	length = 0.0;
	time = 0.0;
	label = 1;
	leaf = false;
	trunk = false;
	include = true;
		
}

/* Get functions */
int Node::getNumber() { return number; }
string Node::getName() { return name; }
double Node::getLength() { return length; }
double Node::getTime() { return time; }
int Node::getLabel() { return label; }
bool Node::getLeaf() { return leaf; }
bool Node::getTrunk() { return trunk; }
bool Node::getInclude() { return include; }

/* Set functions */
void Node::setNumber(int n) { number = n; }
void Node::setName(string n) { name = n; }
void Node::setLength(double n) { length = n; }
void Node::setTime(double n) { time = n; }
void Node::setLabel(int n) { label = n; }
void Node::setLeaf(bool n) { leaf = n; }
void Node::setTrunk(bool n) { trunk = n; }
void Node::setInclude(bool n) { include = n; }
