#pragma once

#include "DirectedGraph.hpp"
#include <vector>
#include <list>
#include <stack>
#include <set>
#include <functional>

typedef vertex_interface my_vertex;
typedef edge_interface my_edge;

struct component_partition {
    std::vector<std::set<my_vertex>> forest;
    DirectedGraph dag = DirectedGraph(0,0,0);
    std::vector<my_vertex> component;
    std::vector<bool> both;
    int n_;

    component_partition(int n) : n_(n/2) {
        forest = {};
        both = {};
        component = std::vector<my_vertex>(n, -1);
    };

    void build_dag(DirectedGraph& D) {
        int size = forest.size();
        dag = DirectedGraph(size, size, D.get_start());
        // this is still O(n+m), we visit each vertex v deg-(v)+1 times.
        for(auto& scc : forest) {
            for(my_vertex v : scc) {
                for(my_vertex w : D.get_neighbors(v)) {
                    int component_v = get_component_id(v);
                    int component_w = get_component_id(w);
                    if(component_v != component_w) {
                        dag.add_edge(component_v, component_w);
                    }
                }
            }
        }
    };

    my_vertex get_component_id(my_vertex v) { return component[v]; };
    std::set<my_vertex> get_component(int id) { return forest[id]; };
    int size() { return forest.size(); };
    void add_component() { 
        forest.push_back({}); 
        both.push_back(false);
    };
    void insert(my_vertex v) {
        forest.back().insert(v);
        component[v] = forest.size()-1;
        if(v < n_) {
            if(component[v] == component[v+n_]) {
                both[component[v]] = true;
            }
        } else {
            if(component[v] == component[v-n_]) {
                both[component[v]] = true;
            }
        }
    };
    bool contains(my_vertex v) { return component[v] != -1; };

    // int reduce_reverse(std::function<int ()> f) {
    //     int acc = 0;
    //     int size = forest.size();
    //     for(int i = size-1; i > -1; --i) {

    //     }
    // };
    bool has_both(int scc) {
        return both[scc];
    };
    std::vector<my_vertex> right(int scc) {
        std::vector<my_vertex> res;
        res.reserve(forest[scc].size());
        for(my_vertex v : forest[scc]) {
            if(v >= n_) {
                res.push_back(v-n_);
            }
        }
        return res;
    };

    std::vector<my_vertex> left(int scc) {
        std::vector<my_vertex> res;
        res.reserve(forest[scc].size());
        for(my_vertex v : forest[scc]) {
            if(v < n_) {
                res.push_back(v);
            }
        }
        return res;
    };
};

void visit(my_vertex v, component_partition& sccs, std::vector<bool>& visited, std::stack<my_vertex>& s, DirectedGraph& D) {
    if(visited[v]) {
        return;
    }
    visited[v] = true;
    std::vector<my_vertex> neighbors = D.get_neighbors(v);
    for(my_vertex w : neighbors) {
        visit(w, sccs, visited, s, D);
    }
    s.push(v);
}
void assign(my_vertex v, component_partition& sccs, DirectedGraph& D_transposed) {
    if(sccs.contains(v)) {
        return;
    }
    sccs.insert(v);
    std::vector<my_vertex> neighbors = D_transposed.get_neighbors(v);
    for(my_vertex w : neighbors) {
        assign(w, sccs, D_transposed);
    }
}