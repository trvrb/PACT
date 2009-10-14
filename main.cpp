/* 	
	(P)osterior (A)nalysis of (C)oalescent (T)rees 
	Copyright 2009 Trevor Bedford <bedfordt@umich.edu>

	This program is designed to interpret and manipulate labeled evolutionary trees.  Statistics 
	regarding the structured coalescent may be calculated.
*/

/*	tree.hh: Copyright 2001-2006 Kasper Peeters <kasper.peeters@aei.mpg.de> */

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

/*	Agenda:
	Fix padTree
	Include migration events in PrintParen
	Is there a way to split coal_tree?
		Seems like a good split would be printing / information functions vs. manipulation functions.
		Also, there are a number of low level functions that would be nice to add.  Things like move_branch.
		What sort of primitives would be necessary to make this work?
	Build command-line options.
	Tests?
*/

// Template class for standard library style tree
#include "tree.hh"

// Class for coalescent nodes within a tree object
#include "node.h"

// Extension of the tree class to deal specifically with coalescent trees
#include "coaltree.h"

// Collects a series of measurements, usually from multiple trees
#include "series.h"

// Input Migrate and Beast tree files and output Mathematica trees and tables of statistics
#include "io.h"

// Set default parameters and modify via parameter file
#include "param.h"

#include <iostream>
using std::cout;
using std::endl;

int main() {
			
	cout << "PACT Copyright 2009 Trevor Bedford" << endl << endl;			
	IO trees;
	trees.treeManip();
	trees.printTree();
	trees.printStatistics();
	trees.printTips();
	trees.printSkylines();

	return 0;
}


