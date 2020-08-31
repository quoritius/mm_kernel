#pragma once
#include <queue>
#include <map>
#include <utility>
#include <set>
#include "Graph.hpp"
#include "DirectedGraph.hpp"
#include "helping_algos.hpp"

typedef vertex_interface my_vertex;
typedef edge_interface my_edge;

int solve_isolated(Graph& G) {
    std::vector<my_vertex> vertices = G.get_vertex_interfaces();
    for(my_vertex v : vertices) {
        if(G.get_deg(v) == 0) {
            G.delete_vertex(v);
        }
    }
    return 0;
}

int solve_externals(Graph& G) {
    int solution = 0;
    std::vector<my_vertex> vertices = G.get_vertex_interfaces();
    std::queue<my_vertex> externals = {};
    // find all deg0 and deg1 vertices
    for(my_vertex v : vertices) {
        int deg = G.get_deg(v);
        if(deg == 0) {
            G.delete_vertex(v);
        } else if(deg == 1) {
            externals.push(v);
        }
    }

    // exhaustively remove external vertices
    while(!externals.empty()) {
        my_vertex v = externals.front();
        externals.pop();

        int deg = G.get_deg(v);
        if(deg == 0) {
            G.delete_vertex(v);
            continue;
        } else if(deg != 1) {
            std::cerr << "\n\n\nDEG01 ERROR!\n\n\n";
            return(-1);
        }

        // there is only one neighbor for v, call it u
        my_vertex u = G.get_neighbors(v).front();
        
        G.delete_vertex(v);

        // neighbors of u might become external after deleting u
        std::vector<my_vertex> u_neighbors = G.get_neighbors(u);
        for(my_vertex u_neighbor : u_neighbors) {
            int deg = G.get_deg(u_neighbor);
            if(deg == 1) {
                G.delete_vertex(u_neighbor);
            } else if(deg == 2) {
                externals.push(u_neighbor);
            }
        }
        G.delete_vertex(u);
        ++solution;
    }
    return(solution);
}

int solve_redexes(Graph& G) {
    auto other_neighbor = [&G](my_vertex v, my_vertex u) -> my_vertex {
        std::vector<my_vertex> neighbors = G.get_neighbors(v);
        return(neighbors[0] == u ? neighbors[1] : neighbors[0]);
    };

    int solution = 0;
    int start = G.get_start();
    std::vector<bool> added(G.get_vertex_size(), false);
    std::queue<my_vertex> redex_centers;
    std::vector<my_vertex> vertices = G.get_vertex_interfaces();
    
    // go through the whole graph, add deg2 to a list/queue
    for(my_vertex v : vertices) {
        if(G.get_deg(v) == 2 && !added[v]) {
            redex_centers.push(v);
            added[v] = true;
        }
    }

    // go through the redex_centers
    while(!redex_centers.empty()) {
        // c is the redex center used now, u and w are the focal vertices
        my_vertex c = redex_centers.front();
        redex_centers.pop();

        // if vertex deleted or not deg2 anymore: next element of list/queue
        if(G.check_deleted(c) || G.get_deg(c) != 2) {
            continue;
        }
        
        // handle the redex_grove now defined by c
        std::queue<vertex_interface> inner_redex_centers;
        inner_redex_centers.push(c);
        // choose any focal vertex, it is the sink now
        my_vertex sink = G.get_neighbors(c)[0];
        //std::cout << "New sink: " << sink << std::endl;
        // also add all other redexes that sink is focus of
        std::vector<my_vertex> sink_neighbors = G.get_neighbors(sink);
        for(my_vertex neighbor : sink_neighbors) {
            if(G.get_deg(neighbor) == 2 && neighbor != c) {
                inner_redex_centers.push(neighbor);
            }
        }
        while(!inner_redex_centers.empty()) {
            c = inner_redex_centers.front();
            inner_redex_centers.pop();

            // if vertex deleted or not deg2 anymore: next element of list/queue
            if(G.check_deleted(c) || G.get_deg(c) != 2) {
                continue;
            }

            my_vertex focus = other_neighbor(c, sink);
            std::vector<my_edge> f_edges = G.get_v_edge_interfaces(focus);
            for(auto e : f_edges) {
                my_vertex neighbor = G.get_neighbor(focus, e);
                // we don't want to copy c's edges, we know both anyways.
                // also don't touch sink
                if(neighbor == c || neighbor == sink) {
                    continue;
                }
                G.move_edge(e, focus, sink);
                // for each moved edge: 
                // check whether the other vertex is a new redex,
                // if yes: add to inner queue
                if(G.get_deg(neighbor) == 2) {
                    inner_redex_centers.push(neighbor);
                    added[neighbor] = true;
                }
            }
            //std::cout << c+start << " matched with " << focus+start << std::endl;
            G.delete_vertex(c);
            G.delete_vertex(focus);
            ++solution;
        }
        if(G.get_deg(sink) == 2 && !added[sink]) {
            redex_centers.push(sink);
            added[sink] = true;
        }
    }
    return(solution);
}

int solve_crowns(Graph& G) {
    Graph B = build_bipartite_graph(G);
    crown_result result = solve_mm_bipartite(B);

    std::vector<my_vertex> vertices = G.get_vertex_interfaces();
    int start = G.get_start();
    int solution = 0;
    // remove the first bunch of integer-valued vertices
    for(my_vertex v : vertices) {
        int res = result.first_LP_solution[v];
        // could actually ask result[v] != 1
        if(res == 0 || res == 2) {
            G.delete_vertex(v);
            solution += (res==2 ? 1 : 0);
            // std::cout << "Deleted " << v+G.get_start() << 
            //     " (" << res << ")" << std::endl;
        }
    }

    // exhaustive search of integer-valued vertices
    DirectedGraph D = build_flow_residual(G, result.R_matched);
    component_partition sccs = kosaraju(D);

    // go thorugh the dag in topological reverse order
    // if scc has no outgoing edges and S_L \cap S_R = \emptyset, 
    // then remove all vertices in it, increase solution by size of S_R
    // remove all ingoing edges for scc.
    int size = sccs.size();
    for(int scc = size-1; scc >= 0; --scc) {
        if(!sccs.has_both(scc) and sccs.dag.get_neighbors(scc).empty()) {
            sccs.dag.delete_vertex(scc);
            solution += G.delete_vertices(sccs.right(scc));
            G.delete_vertices(sccs.left(scc));
        }
    }

    return solution;
}

int solve_relaxed_crowns(Graph& G) {
    int solution = 0;
    Graph B = build_bipartite_graph(G);
    crown_result result = solve_mm_bipartite(B);
    DirectedGraph D = build_flow_residual(G, result.R_matched);
    component_partition sccs = kosaraju(D);

    // ** Use the nearly tail strongly connected components **
    // go through the dag in topological reverse order
    // if scc has exactly one outgoing edge and S_L \cap S_R = \emptyset,
    // then handle it right away.
    int size = sccs.size();
    for(int scc = size-1; scc >= 0; --scc) {
        // get scc neighbors
        std::vector<vertex_interface> scc_out = D.get_neighbors(sccs.get_component(scc));
        auto last = std::unique(scc_out.begin(), scc_out.end());
        scc_out.erase(last, scc_out.end()); 
        if(!sccs.has_both(scc) and scc_out.size() == 1) {
            vertex_interface x = scc_out[0] - sccs.n_;
            std::vector<vertex_interface> H = sccs.right(scc);
            std::vector<vertex_interface> I = sccs.left(scc);
            G.delete_vertices(I);
            G.add_edges_pl(x, G.get_neighbors(std::set<vertex_interface>(H.begin(), H.end())));
            solution += G.delete_vertices(H);
        }
    }
    if(solution != 0) {
        return solution;
    }

    // ** Use the strong articulation points **
    DirectedGraph D_contracted = contract(D, result.R_matched, G);
    sccs = kosaraju(D_contracted);
    size = sccs.size();
    for(int scc = size-1; scc >= 0; --scc) {
        std::set<vertex_interface> s = sccs.get_component(scc);
        std::vector<vertex_interface> scc_vertices(s.size());
        std::copy(s.begin(), s.end(), scc_vertices.begin());
        std::set<vertex_interface> A = get_strong_artpoints(D_contracted, scc_vertices);
        for(vertex_interface x : A) {
            if(G.is_deleted(x)) {
                continue;
            }
            DirectedGraph D_x = cut_graph(D, x);
            component_partition sccs_x = kosaraju(D_x);
            int size_x = sccs_x.size();
            for(int scc_x = size_x-1; scc_x >= 0; --scc_x) {
                if(sccs_x.has_both(scc_x)) {
                    continue;
                }
                std::vector<my_vertex> scc_x_out = D.get_neighbors(sccs_x.get_component(scc_x));
                auto last = std::unique(scc_x_out.begin(), scc_x_out.end());
                scc_x_out.erase(last, scc_x_out.end()); 
                // check S_L cap S_R = emptyset
                if(scc_x_out.size() == 1) {
                    std::vector<vertex_interface> H = sccs_x.right(scc_x);
                    std::vector<vertex_interface> I = sccs_x.left(scc_x);
                    G.delete_vertices(I);
                    G.add_edges_pl(x, G.get_neighbors(std::set<vertex_interface>(H.begin(), H.end())));
                    solution += G.delete_vertices(H);
                    return(solution);
                }
            }
        }
    }
    return solution;
}

int find_concurrent_sets(Graph& G) {
    std::multimap<std::set<my_vertex>, my_vertex> neighborhoods;
    std::vector<my_vertex> vertices = G.get_vertex_interfaces();
    for(my_vertex v : vertices) {
        std::vector<my_vertex> neighbors = G.get_neighbors(v);
        neighborhoods.emplace(std::set<my_vertex>(neighbors.begin(), neighbors.end()), v);
    }
    for(auto& [neighborhood, v] : neighborhoods) {
        int n = neighborhoods.count(neighborhood);
        int size = neighborhood.size();
        if(n >= size-1) {
            auto [begin, end] = neighborhoods.equal_range(neighborhood);
            std::cout << "NEW CONCURRENT SET(" << n << "):" << std::endl;
            for(auto i = begin; i != end; ++i) {
                std::cout << i->second << ", ";
            }
            std::cout << std::endl;
            std::cout << "with neighborhood(" << size << "): " << std:: endl;
            for(my_vertex neighbor : neighborhood) {
                std::cout << neighbor << ", ";
            }
            std::cout << std::endl << std::endl;
        }
    }
    return 0;
}

int solve_concurrent_sets(Graph& G, bool print = false) {
    int result = 0;

    std::multimap<std::set<my_vertex>, my_vertex> neighborhoods;
    std::vector<my_vertex> vertices = G.get_vertex_interfaces();
    for(my_vertex v : vertices) {
        std::vector<my_vertex> neighbors = G.get_neighbors(v);
        neighborhoods.emplace(std::set<my_vertex>(neighbors.begin(), neighbors.end()), v);
    }
    for(auto& [neighborhood_const, v] : neighborhoods) {
        // check deleted neighborhood
        std::list<my_vertex> neighborhood(neighborhood_const.begin(), neighborhood_const.end());
        neighborhood.remove_if([&](my_vertex v)->bool{return G.check_deleted(v);});

        // check deleted twins
        auto [begin, end] = neighborhoods.equal_range(neighborhood_const);
        std::list<my_vertex> concurrents;
        auto it = begin;
        while(it != end) {
            concurrents.push_back(it->second);
            ++it;
        }
        concurrents.remove_if([&](my_vertex v){return G.check_deleted(v);});

        // check n == size-1 or n > size-1
        // handle accordingly
        // n is number of concurrents, l is the size of their neighborhood
        int n = concurrents.size();
        int l = neighborhood.size();

        if(print) {
            if(n >= l-1) {
                std::cout << "DELETED" << std::endl;
                for(my_vertex v : concurrents) {
                    std::cout << v << " ";
                }
                std::cout << std::endl << "with" << std::endl;
                for(my_vertex w : neighborhood) {
                    std::cout << w << " ";
                }
                std::cout << std::endl << std::endl;
            }
        }

        // |H| = |I| + 1
        if(n == l-1) {
            for(my_vertex v : concurrents) {
                G.delete_vertex(v);
            }
            // choose first neighbor as sink, for others w: move w's edges to sink and remove w
            my_vertex sink = neighborhood.front();
            neighborhood.pop_front();
            for(my_vertex w : neighborhood) {
                for(my_edge& e : G.get_v_edge_interfaces(w)) {
                    G.move_edge(e, w, sink);
                }
                G.delete_vertex(w);
            }
            result += l-1;
        // |H| < |I|+1
        } else if(n > l-1) {
            for(my_vertex v : concurrents) {
                G.delete_vertex(v);
            }
            for(my_vertex w : neighborhood) {
                G.delete_vertex(w);
            }
            result += l;
        }
    }
    return result;
}