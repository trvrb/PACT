/* 	
	(P)osterior (A)nalysis of (C)oalescent (T)rees 
	Copyright 2009 Trevor Bedford <bedfordt@umich.edu>

	This program is designed to interpret and manipulate labeled phylogenetic trees.  Statistics 
	regarding the structured coalescent may be calculated.
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
	Finish with parameters (skylines need implementing).
*/

/*
	tree.hh: Copyright 2001-2006 Kasper Peeters <kasper.peeters@aei.mpg.de>
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

int main() {
			
	IO trees;
	trees.treeManip();
	trees.printHPTree();
	trees.printStatistics();
	trees.printTips();
	trees.printSkylines();

	return 0;
}


