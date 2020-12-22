
## Info:

Feel free to use these lines as you wish.
Computes a hierarchical clustering.
DivHi pronounce "Dive high".  
Maybe should be changed into "High Dive": https://en.wikipedia.org/wiki/High_diving

## To compile:

"make"

## To execute:

"./divhi edgelist.txt tree.txt [bisection]".

"edgelist.txt" should contain an undirected unweighted graph: one edge on each line (two unsigned long (nodes' ID)) separated by a space.

"tree.txt": will contain the hierarchy. One set of nodes on each line: "n : node1 node2 ... nodek".
- "n" is the ID of the tree's node: 1 for the root and the ID of the children of node "n" are "2n" and "2n+1".
- "node1", "node2", ..., "nodek" are the IDs of nodes of the original graph in the set of nodes associated with tree's node "n".
- nodes that have degree zero are automatically removed
- only tree's nodes associated to a set of nodes containing at least 2 nodes are represented.

option "bisection": the algorithm to use, default is 1.
- 0: random bisection
- 1: greedy spasest cut
- 2: greedy densest cut
- 3: greedy max cut
- 4: greedy min cut
- 5: greedy min cut in the complement
- 6 minimum of 4 and 5

## 

## Initial contributor:
Maximilien Danisch  
December 2018  
http://bit.ly/danisch  
maximilien.danisch@gmail.com
