
....................................................................

This code was written by Long T. Le at Rutgers University.  
The Eliassi Lab owns the copyright to it.
This code is associated with this paper:

  "MET: A Fast Algorithm for Minimizing Propagation 
   in Large Graphs with Small Eigen-Gaps" 
   by Long T. Le, Tina Eliassi-Rad, and Hanghang Tong, 
   appeared in the Proceedings of the 2015 SIAM Conference 
   on Data Mining (SDM'15), Vancouver, Canada, April 2015.  
   URL: http://eliassi.org/papers/le-sdm15.pdf

....................................................................

This code was written and tested in Matlab 2014a.

To start the code, run Main_MET.m, which loads a graph and calls 
the function “IE_DeltaLam_k_MET.m”, which has the main algorithm.

There are two input parameters that you can change:

  1. The input graph 
     - We provide a sample input in sample-graphs/sample.csv.

  2. The edge-deletion budget, k

The code outputs the percentage decrease in the leading eigenvalue 
of the adjacency matrix: 

  100*(lambda1_before_edge_deletion - lambda1_after_k_edge_deletions)
....................................................................

