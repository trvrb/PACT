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

Node::Node(int inName) {
	
	name = inName;
	length = 0.0;
	time = 0.0;
	label = 1;
	leaf = false;
	trunk = false;
		
}

/* Simple get functions */
int Node::getName() { return name; }
double Node::getLength() { return length; }
double Node::getTime() { return time; }
int Node::getLabel() { return label; }
bool Node::getLeaf() { return leaf; }
bool Node::getTrunk() { return trunk; }

