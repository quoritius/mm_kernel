#include <iostream>
#include <fstream>
#include <cstdio>
#include <set>
#include <utility>
#include "Graph.hpp"

typedef vertex_interface my_vertex;
typedef edge_interface my_edge;

/*
 * Reducing a Maximum Cardinality Matching instance G to a 
 * Minimum Weight Perfect Matching instance G', to be able to use BLOSSOM V.
 * V(G') is V(G) and one more copy of V(G)( call it V'(G)). 
 * E(G') is E(G) with cost -1 and copies of them on V'(G) with cost 0, and 
 * edges between vertices and their copies with cost 0.
 */
void weightify(Graph& G) {
    // G.clean_vertices();
    int n = G.get_vertex_size();
    int m = G.get_actual_edge_size();
    std::cout << 2*n << " " << n+2*m << std::endl;

    std::vector<my_vertex> vertices = G.get_vertex_interfaces();
    for(my_vertex v : vertices) {
        std::vector<my_edge> edges = G.get_v_edge_interfaces(v);
        for(my_edge e : edges) {
            my_vertex w = G.get_neighbor(v, e);
            // E(G)
            std::cout << v << " " << w << " " << -1 << std::endl;
            // E'(G)
            std::cout << v+n << " " << w+n << " " << 0 << std::endl;
            G.delete_edge(e);
        }
    }
    for(my_vertex v : vertices) {
        // Edges between v from V(G) and its copy v' from V'(G)
        std::cout << v << " " << v+n << " " << 0 << std::endl;
    }
}

int main() {
    int n, m, start;
    // freopen("input/big/roadNet-PA.txt", "r", stdin);
    std::cin >> start >> n >> m;
    Graph G = Graph(n, m, start);
    // fclose(stdin);
    weightify(G);
    return 0;
}