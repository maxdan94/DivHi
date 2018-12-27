
## Info:

Feel free to use these lines as you wish.
Computes a hierarchical clustering

## To compile:

"make"

## To execute:

"./hierarchy edgelist.txt tree.txt [bisection]".

"edgelist.txt" should contain the undirected unweighted graph: one edge on each line (two unsigned long (nodes' ID)) separated by a space.

"tree.txt": will contain the hierarchy. One set of nodes on each line: "n : node1 node2 ... nodek".
- "n" is the id of the tree's node: 1 for the root and the id of the children of node "n" are "2n" and "2n+1".
- "node1", "node2", ..., "nodek" are the id of nodes of the original graph  ontained in the set of nodes associated with tree's node "n".
- nodes that have degree zero at some point are automaticaly removed
- only tree's nodes associated to a set of nodes containing at least 2 nodes are represented.

option "bisection": the algorithm to use, default is 1.
- 0: random bisection
- 1: greedy spasest cut
- 2: greedy densest cut
- 3: greedy Maxcut


## Initial contributor:
Maximilien Danisch  
Technical consultants: Ayan Bhowmick, Bivas Mitra and Jean-Loup Guillaume  
December 2018  
http://bit.ly/danisch  
maximilien.danisch@gmail.com
