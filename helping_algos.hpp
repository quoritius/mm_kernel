#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <stack>
#include <set>
#include <cassert>


#include "Graph.hpp"
#include "DirectedGraph.hpp"
#include "components.hpp"
#include "dominators.hpp"

typedef vertex_interface my_vertex;
typedef edge_interface my_edge;

struct crown_result {
    std::vector<my_vertex> R_matched;
    std::vector<int> first_LP_solution;
};

// builds the bipartite graph B to get an BIP2 solution of Vertex Cover
// The vertices of B are divided in V_L={0..n-1} and V_R{n..2*n-1}
// If G had an edge (v, w), B has two edges (v, w+n), (v+n, w)
// Could return B by copy, complexity would not change?? Performance matters.
Graph build_bipartite_graph(Graph& G) {
    int n = G.get_vertex_size();
    int m = G.get_edge_size();
    int start = G.get_start();
    Graph B(n*2, m*2, start, false);
    std::vector<my_vertex> vertices = G.get_vertex_interfaces();
    for(my_vertex v : vertices) {
        std::vector<my_edge> edges = G.get_v_edge_interfaces(v);
        for(my_edge e : edges) {
            B.add_edge(e.first, e.second + n);
            B.add_edge(e.first + n, e.second);
        }
    }
    return B;
}

// solve bipartite maximum matching with hopcroft karp.
// Input graph B should be bipartite and have exactly 2*n vertices. 
// 0..n-1 should be in one partition, n..2*n-1 in the other partition.
//
// result[i] = 0 if i not used for VC in l and r
// result[i] = 1 if i used for VC only in l or only in r
// result[i] = 2 if i used for VC in both l and r
crown_result solve_mm_bipartite(Graph& B){
    int n = B.get_vertex_size() / 2;
    int start = B.get_start();
    int solution = 0;
    std::set<my_vertex> L_free;
    std::set<my_vertex> R_free;
    std::vector<my_vertex> L_matched(n, -1);
    std::vector<my_vertex> R_matched(n, -1);
    std::list<my_vertex> F;
    std::vector<bool> visited(n*2, false);
    // needed_depth is actually not needed
    int needed_depth = 0;
    // optimization: used_depth can substitute visited
    std::vector<int> used_depth(n*2, 0);

    std::vector<std::vector<my_vertex>> adj_list(n*2, std::vector<my_vertex>());
    for(my_vertex v = 0; v < n*2; v++) {
        adj_list[v] = B.get_neighbors(v);
    }

    auto is_r = [&](my_vertex v)->bool{ return v>=n; };
    auto is_l = [&](my_vertex v)->bool{ return v<n; };
    auto is_l_and_free = [&](my_vertex v) -> bool {
        return is_l(v) ? (L_matched[v]==-1) : false;
    };

    // Could do closure of F to not copy it, using the reference instead.
    // Layers: (1 layer) L->R (2 layer) R->L (3 layer) ... 
    //       .(odd layer) L->R (even layer) R->L (odd layer) L->R (even layer).
    // use matching when R->L
    auto bfs = [&]() {
        bool found = false;
        visited = std::vector<bool>(n*2, false);
        used_depth = std::vector<int>(n*2, 0);
        F.clear();
        needed_depth = 0;
        std::queue<my_vertex> q;  
        std::queue<int> levels;
        for(my_vertex l : L_free) {
            visited[l] = true;
            q.push(l);
            levels.push(1);
        }
        my_vertex v;
        int layer;
        while(!q.empty()) {
            v = q.front();
            q.pop();
            layer = levels.front();
            levels.pop();
            used_depth[v] = layer;
            bool even = is_r(v);
            if(even) {
                if(found) {
                    if(R_matched[v-n]==-1) {
                        F.push_back(v);
                        needed_depth = layer;
                    }
                } else if(!visited[R_matched[v-n]]) {
                    visited[R_matched[v-n]] = true;
                    q.push(R_matched[v-n]);
                    levels.push(layer+1);
                }
            } else {
                std::vector<my_vertex> neighbors = adj_list[v];
                for(my_vertex u : neighbors) {
                    if(!visited[u] && L_matched[v] != u) {
                        // not even, v is in L and u is in R
                        if(R_matched[u-n] == -1) {
                            found = true;
                        } 
                        visited[u] = true;
                        q.push(u);
                        levels.push(layer+1);
                    }
                }
            }
        }
        return;
    };

    // v comes from F and as such from R
    // trying to use y combinator (?) here
    // maybe just use explicit type with std::function (might get overhead)
    // Go the layers in reverse order from F.
    // .(even layer) R->L (odd layer) L->R (even layer) R->L (odd layer).
    // use matching when L->R
    auto dfs = [&](my_vertex v) -> bool {
        return ([&](auto dfs_util, my_vertex v) {
            std::stack<my_vertex> s;
            s.push(v);
            my_vertex first = v;
            return dfs_util(dfs_util, v, needed_depth);
        })(
            [&](auto dfs_util,  my_vertex v, int depth) -> bool {
                bool even = is_r(v);
                if(is_l_and_free(v)) {
                    L_free.erase(v);
                    return true;
                }
                // if is_r(v) (even), we need to check all neighbors,
                // that were visited by dft and are not "devisited"
                if(even) {
                    for(my_vertex w : adj_list[v]) {
                        if(visited[w] && R_matched[v-n] != w) {
                            if(depth == 2 && !is_l_and_free(w)) {
                                continue;
                            } else if (used_depth[w] != used_depth[v]-1) {
                                continue;
                            }
                            // mark w to not visit it later in the path again.
                            visited[w] = false;
                            if(dfs_util(dfs_util, w, depth-1)) {
                                L_matched[w] = v;
                                R_matched[v-n] = w;
                                visited[v] = visited[w] = false;
                                return true;
                            }
                        }
                    }
                    visited[v] = false;
                    return false;
                } else {
                    // if is_l(v), we can use L_matched[v] (odd layer)
                    // odd layer can actually be handled from the even layer.
                    my_vertex w = L_matched[v];
                    // marking not needed, v uses only w, and only v can use w.
                    if(visited[w]) {
                        return dfs_util(dfs_util, w, depth-1);
                    } else {
                        return false;
                    }
                }
            },
            v // v should come from F.
        );
    };
    std::vector<my_vertex> vertices = B.get_vertex_interfaces();
    for(my_vertex v : vertices) {
        if(B.get_deg(v) == 0) {
            continue;
        }
        if(is_r(v)) {
            R_free.insert(v);
        } else {
            L_free.insert(v);
        }
    }
    int step = 1;
    while(!R_free.empty() && !L_free.empty()) {
        // std::cout << "start step " << step << std::endl;
        bfs();
        // std::cout << "depth needed: " << needed_depth << std::endl;
        if(F.empty()) {
            break;
        }
        for(my_vertex v : F) {
            if(dfs(v)) {
                R_free.erase(v);
                ++solution;
            }
        }
        ++step;
    }
    
    // std::cout << "Solution size is: " << solution << std::endl;
    // for(my_vertex v = 0; v < n; v++) {
    //     std::cout << v+start << " matched with " << L_matched[v]+start-n << std::endl;
    // }

    std::vector<int> result(n, 0);
    for(int i = 0; i < n; i++) {
        if(L_matched[i] != -1) {
            ++result[i];
        }
    }

    auto konig_bft = [&]() {
        visited = std::vector<bool>(n*2, false);
        std::queue<my_vertex> q;
        for(my_vertex l : L_free) {
            visited[l] = true;
            q.push(l);
        }
        my_vertex v;
        while(!q.empty()) {
            v = q.front();
            q.pop();
            bool even = is_r(v);
            if(even) {
                ++result[v-n];
                if(R_matched[v-n]==-1) {
                    // SHOULD NOT HAPPEN!
                } else if(!visited[R_matched[v-n]]) {
                    visited[R_matched[v-n]] = true;
                    q.push(R_matched[v-n]);
                }
            } else {
                result[v] -= L_matched[v]!=-1;
                std::vector<my_vertex> neighbors = adj_list[v];
                for(my_vertex u : neighbors) {
                    if(!visited[u] && L_matched[v] != u) {
                        // not even, v is in L and u is in R
                        if(R_matched[u-n] == -1) {
                            int a = 0;
                            // SHOULD NOT HAPPEN! This would be part of an augmenting path.
                        } 
                        visited[u] = true;
                        q.push(u);
                    }
                }
            }
        }
        return;
    };

    // for(int i = 0; i < n; i++) {
    //     // u is from R, v is from L
    //     my_vertex u = L_matched[i];
    //     my_vertex v = R_matched[i];
    //     result[i] = (u==-1 ? 0 : 1) + (v==-1 ? 0 : 1);
    // }
    konig_bft();
    return {R_matched, result};
}

// TODO: move the konig to solve_crown, return R_matched from solve_mm_bipartite
std::vector<int> solve_mm_bipartite_(Graph& B){
    int n = B.get_vertex_size() / 2;
    int start = B.get_start();
    int solution = 0;
    std::set<my_vertex> L_free;
    std::set<my_vertex> R_free;
    std::vector<my_vertex> L_matched(n, -1);
    std::vector<my_vertex> R_matched(n, -1);
    std::list<my_vertex> F;
    std::vector<bool> visited(n*2, false);
    // needed_depth is actually not needed
    int needed_depth = 0;
    // optimization: used_depth can substitute visited
    std::vector<int> used_depth(n*2, 0);

    std::vector<std::vector<my_vertex>> adj_list(n*2, std::vector<my_vertex>());
    for(my_vertex v = 0; v < n*2; v++) {
        adj_list[v] = B.get_neighbors(v);
    }

    auto is_r = [&](my_vertex v)->bool{ return v>=n; };
    auto is_l = [&](my_vertex v)->bool{ return v<n; };
    auto is_l_and_free = [&](my_vertex v) -> bool {
        return is_l(v) ? (L_matched[v]==-1) : false;
    };

    // Could do closure of F to not copy it, using the reference instead.
    // Layers: (1 layer) L->R (2 layer) R->L (3 layer) ... 
    //       .(odd layer) L->R (even layer) R->L (odd layer) L->R (even layer).
    // use matching when R->L
    auto bfs = [&]() {
        bool found = false;
        visited = std::vector<bool>(n*2, false);
        used_depth = std::vector<int>(n*2, 0);
        F.clear();
        needed_depth = 0;
        std::queue<my_vertex> q;  
        std::queue<int> levels;
        for(my_vertex l : L_free) {
            visited[l] = true;
            q.push(l);
            levels.push(1);
        }
        my_vertex v;
        int layer;
        while(!q.empty()) {
            v = q.front();
            q.pop();
            layer = levels.front();
            levels.pop();
            used_depth[v] = layer;
            bool even = is_r(v);
            if(even) {
                if(found) {
                    if(R_matched[v-n]==-1) {
                        F.push_back(v);
                        needed_depth = layer;
                    }
                } else if(!visited[R_matched[v-n]]) {
                    visited[R_matched[v-n]] = true;
                    q.push(R_matched[v-n]);
                    levels.push(layer+1);
                }
            } else {
                std::vector<my_vertex> neighbors = adj_list[v];
                for(my_vertex u : neighbors) {
                    if(!visited[u] && L_matched[v] != u) {
                        // not even, v is in L and u is in R
                        if(R_matched[u-n] == -1) {
                            found = true;
                        } 
                        visited[u] = true;
                        q.push(u);
                        levels.push(layer+1);
                    }
                }
            }
        }
        return;
    };

    // v comes from F and as such from R
    // trying to use y combinator (?) here
    // maybe just use explicit type with std::function (might get overhead)
    // Go the layers in reverse order from F.
    // .(even layer) R->L (odd layer) L->R (even layer) R->L (odd layer).
    // use matching when L->R
    auto dfs = [&](my_vertex v) -> bool {
        return ([&](auto dfs_util, my_vertex v) {
            std::stack<my_vertex> s;
            s.push(v);
            my_vertex first = v;
            return dfs_util(dfs_util, v, needed_depth);
        })(
            [&](auto dfs_util,  my_vertex v, int depth) -> bool {
                bool even = is_r(v);
                if(is_l_and_free(v)) {
                    L_free.erase(v);
                    return true;
                }
                // if is_r(v) (even), we need to check all neighbors,
                // that were visited by dft and are not "devisited"
                if(even) {
                    for(my_vertex w : adj_list[v]) {
                        if(visited[w] && R_matched[v-n] != w) {
                            if(depth == 2 && !is_l_and_free(w)) {
                                continue;
                            } else if (used_depth[w] != used_depth[v]-1) {
                                continue;
                            }
                            // mark w to not visit it later in the path again.
                            visited[w] = false;
                            if(dfs_util(dfs_util, w, depth-1)) {
                                L_matched[w] = v;
                                R_matched[v-n] = w;
                                visited[v] = visited[w] = false;
                                return true;
                            }
                        }
                    }
                    visited[v] = false;
                    return false;
                } else {
                    // if is_l(v), we can use L_matched[v] (odd layer)
                    // odd layer can actually be handled from the even layer.
                    my_vertex w = L_matched[v];
                    // marking not needed, v uses only w, and only v can use w.
                    if(visited[w]) {
                        return dfs_util(dfs_util, w, depth-1);
                    } else {
                        return false;
                    }
                }
            },
            v // v should come from F.
        );
    };
    std::vector<my_vertex> vertices = B.get_vertex_interfaces();
    for(my_vertex v : vertices) {
        if(B.get_deg(v) == 0) {
            continue;
        }
        if(is_r(v)) {
            R_free.insert(v);
        } else {
            L_free.insert(v);
        }
    }
    int step = 1;
    while(!R_free.empty() && !L_free.empty()) {
        // std::cout << "start step " << step << std::endl;
        bfs();
        // std::cout << "depth needed: " << needed_depth << std::endl;
        if(F.empty()) {
            break;
        }
        for(my_vertex v : F) {
            if(dfs(v)) {
                R_free.erase(v);
                ++solution;
            }
        }
        ++step;
    }

    return R_matched;
}

DirectedGraph build_flow_residual(Graph& G, std::vector<my_vertex>& R_matched) {
    int n = G.get_vertex_size();
    int m = G.get_edge_size();
    int start = G.get_start();
    DirectedGraph D = DirectedGraph(2*n, 2*m, start);
    std::vector<my_vertex> vertices = G.get_vertex_interfaces();
    for(my_vertex v : vertices) {
        std::vector<my_vertex> neighbors = G.get_neighbors(v);
        for(my_vertex w : neighbors) {
            D.add_edge(v, w+n);
        }
        my_vertex v_matched = R_matched[v];
        if(v_matched != -1) {
            D.add_edge(v+n, v_matched);
        }
    }
    return D;
}

component_partition find_sccs(DirectedGraph& D) {
    int n = D.get_vertex_size();
    DirectedGraph D_transposed = D.get_transpose();
    std::stack<my_vertex> s;
    component_partition sccs(n);
    std::vector<bool> visited(n, false);

    std::vector<my_vertex> vertices = D.get_vertex_interfaces();
    for(my_vertex v : vertices) {
        visit(v, sccs, visited, s, D);
    }
    while(!s.empty()) {
        my_vertex v = s.top();
        s.pop();
        if(!sccs.contains(v)) {
            sccs.add_component();
            assign(v, sccs, D_transposed);
        }
    }
    sccs.build_dag(D);
    return sccs;
}

// Iterative version of Kosaraju's algorithm to find strongly connected components.
// Thanks to sof @stackexchange for his post
// https://codereview.stackexchange.com/questions/242161/recursive-and-iterative-implementation-on-kosaraju-algorithm
component_partition kosaraju(DirectedGraph& D) {
    int n = D.get_vertex_size();
    DirectedGraph D_transposed = D.get_transpose();
    std::stack<my_vertex> s;
    component_partition sccs(n);
    std::vector<bool> visited(n, false);
    std::vector<int> skip(n, 0);
    std::stack<my_vertex> l;
    bool done = false;
    // postorder DFS on D to push root vertices to L
    std::vector<my_vertex> vertices = D.get_vertex_interfaces();
    for(my_vertex v : vertices) {
        if(!visited[v]) {
            visited[v] = true;
            s.push(v);
            while(!s.empty()) {
                v = s.top();
                done = true;
                std::vector<my_vertex> neighbors = D.get_neighbors(v, skip[v]);
                for(my_vertex w : neighbors) {
                    skip[v]++;
                    if(!visited[w]) {
                        visited[w] = true;
                        done = false;
                        s.push(w);
                        break;
                    }
                }
                if(done) {
                    s.pop();
                    l.push(v);
                }
            }
        }
    }
    skip = std::vector<int>(n, 0);
    // postorder DFS on D_transposed to pop root vertices from L and mark SCCs
    while(!l.empty()) {
        my_vertex r = l.top(); l.pop();
        s.push(r);
        if(visited[r]) {
            visited[r] = false;
            sccs.add_component();
            sccs.insert(r);
        }
        while(!s.empty()) {
            my_vertex v = s.top();
            done = true;
            std::vector<my_vertex> neighbors = D_transposed.get_neighbors(v, skip[v]);
            for(my_vertex w : neighbors) {
                skip[v]++;
                if(visited[w]) {
                    visited[w] = false;
                    done = false;
                    s.push(w);
                    sccs.insert(w);
                    break;
                }
            }
            if(done) {
                s.pop();
            }
        }
    }
    sccs.build_dag(D);
    return sccs;
}


std::set<vertex_interface> get_strong_artpoints(
        DirectedGraph &D, std::vector<vertex_interface> &scc
) {
    std::set<vertex_interface> artpoints;
    // get scc-subgraph (just the edges inside the scc are enough) of D
    // thus D_ := D[scc]
    // Keep in mind we have other vertex-labels in D[scc] than in D
    // store the new labels to translate edges from D to D_
    int i = 0;
    std::map<int, int> labels;
    for(auto& v : scc) {
        if(!D.is_deleted(v)) {
            labels[v] = i;
            ++i;
        }
    }
    std::vector<edge_interface> E_ = D.get_multiple_edge_interfaces(scc);
    DirectedGraph D_(scc.size(), E_.size(), 0);
    for(edge_interface e : E_) {
        D_.add_edge(labels[e.first], labels[e.second]);
    }
    DirectedGraph D_t = D_.get_transpose();
    // choose root
    int root = 0;
    // check if root is strong artpoint
    DirectedGraph D_r = D_.get_copy();
    D_r.delete_vertex(root);
    component_partition sccs_r = find_sccs(D_r);
    if(sccs_r.size() != 1) {
        artpoints.insert(scc[root]);
    }
    // get dominators in D[scc]
    Dominator_Solver DomSol(root, scc.size(), D_.get_sorted_neighbor_set());
    std::vector<vertex_interface> dom = DomSol.get_doms();
    // get dominators in D_T[scc]
    Dominator_Solver DomSol_t(root, scc.size(), D_t.get_sorted_neighbor_set());
    std::vector<vertex_interface> dom_t = DomSol_t.get_doms();
    // return union (+?root)
    for(auto v : dom) {
        if(v != root and v >= 0) {
            artpoints.insert(scc[v]);
        }
    }
    for(auto v : dom_t) {
        if(v != root and v >= 0) {
            artpoints.insert(scc[v]);
        }
    }
    return artpoints;
}

std::set<vertex_interface> get_strong_artpoints(
        DirectedGraph &D, std::vector<vertex_interface> &&scc
) {
    return get_strong_artpoints(D, scc);
}

DirectedGraph cut_graph(DirectedGraph& D, vertex_interface x) {
    int n = D.get_vertex_size()/2;
    DirectedGraph D_x = D.get_copy();
    x += n;
    std::vector<edge_interface> Ex_in = D.get_v_in_edge_interfaces(x);
    std::vector<vertex_interface> x_out = D.get_neighbors(x);
    // If the assert fails, then either D or x are built wrong.
    // A is a vertex from the right side of D, and has thus exactly one outgoing edge.
    if(x_out.size() != 1) {
        int a = 0;
    }
    assert(x_out.size() == 1);
    for(edge_interface e : Ex_in) { 
        if(e.first != x_out[0]) {
            D_x.delete_edge(e);
        } else continue;
    }
    return D_x;
}

DirectedGraph contract(DirectedGraph& D, std::vector<my_vertex>& R_matched, Graph& G) {
    int n = D.get_vertex_size()/2;
    int m = D.get_edge_size() - 2*n;
    DirectedGraph C(n, m, D.get_start());
    for(vertex_interface v = 0; v < n; ++v) {
        if(G.is_deleted(v)) {
            C.delete_vertex(v);
            continue;
        }
        std::vector<vertex_interface> W_out = D.get_neighbors(R_matched[v]);
        for(vertex_interface w : W_out) {
            if(v != w-n) {
                C.add_edge(v, w-n);
            }
        }
    }

    return C;
}