The program PACT computes a variety of statistics from a sample of genealogical trees. It is meant 
to extend the functionality of already existing coalescent inference programs such as Migrate, 
BEAST, IM  and LAMARC. PACT reads genealogies in NEWICK format and performs various operations on 
these genealogies. PACT is designed to work with both structured genealogies and also with 
genealogies assembled from temporally spaced sequence data.

The functionality of PACT is highly modular, relying on combinations of tree manipulation operations 
and summary statistics to produce useful results. For example, the operation tmrca returns the TMRCA 
of the entire tree, but when combined with the prune to label, it returns the TMRCA of samples with 
a specified label. In another example, diversity at a specified point in time may be calculated with 
a combination of time slice and diversity.

On UNIX systems you can compile with:

    make

PACT requires two files to run:

    in.param (parameter listing)
    in.trees (NEWICK trees)
  
Available parameters can be found in parameters.txt.  Full documentation can be found in 
pact_manual.pdf.

Copyright 2011 Trevor Bedford.  Distributed under the GPL v3.
