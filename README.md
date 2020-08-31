# mm_kernel
This is the kernelization program used for the bachelor's thesis 
"The Effect of Data Reduction Rules for the Maximum Matching Problem".

The goal is to minimize the input graph while partly solving the 
maximum matching problem on it.

There are 4 Data Reductiton Rules implemented:
* deleting vertices with degree 0 and 1 (external vertices),
* contracting vertices with degree 2,
* deleting crowns,
* contracting relaxed crowns.

The weightify program is a tool that helps reducing the maximum matching 
problem to the minimum weight perfect matching, to prepare the graphs for
using the Blossom V solver (https://link.springer.com/article/10.1007/s12532-009-0002-8).

# Expected Graph Format
The input data should be of the following format:
start, n, m  
from.1 to.1
from.2 to.2
...
from.m to.m

Where everything is a number:
* start is the minimum vertex label, 
* (n+start-1) is the maximum vertex label, 
* m is the number of edges, and
* from.x to.x is an edge with from.x and to.x being vertices with labels between start and (n+start-1).

# Usage
There are 5 modes implemented, namely:
* 1 exhaustively applying  together rules regarding external vertices and redexes,
* 2 same as 1 but afterwards removing all crowns,
* 3 same as 2 but afterwards contracting relaxed crowns until there are not found enough,
* 4 removing all crowns only,
* 5 removing all crowns and afterwards contracting relaxed crowns until there are not found enough.

The program can be compiled using the favorite C++ compiler with C++17 support, ex.g.
``` g++ -g main.cpp Graph.cpp -o kernelizer ```

To run a specific mode $i on some input graph "graph.txt", one should execute
(assuming the code was compiled to the program ```kernelizer```):
``` ./kernelizer $i graph.txt solution.txt outgraph.txt ```

The results will be:
* a file "solution.txt" in which is displayed how much each of the rules decreased the solution size
    (in the same order as explained, with the last number being the sum),
* a file "outgraph.txt" with the remaining graph (kernel), of same format as the input graph.
* The time needed for the kernelization in the standard output.

## weightify 
The weightify tool can also be compiled using ``` g++ -g weightify.cpp Graph.cpp -o weightify ```.
It takes the graph with the mentioned format as the standard input, and gives back the weighted 
instance as the standard output (without start, as it is not needed by BlossomV), so use it like this:
``` ./weightify < graph_in.txt > graph_out.txt ```