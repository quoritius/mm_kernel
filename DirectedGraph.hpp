#pragma once
#include <set>
#include "Graph.hpp"


class DirectedGraph : Graph {
public:
    DirectedGraph(int n, int m, int start) : Graph::Graph(n, m, start, false) {
        deg_in = std::vector<int>(n, 0);

        vertices.reserve(vertex_size);
        adj_in_list.reserve(vertex_size);
        for(int i = 0; i < vertex_size; ++i) {
            adj_in_list.push_back(std::list<int>());
            neighbor_in_set.push_back({});
        }
    }

    void init() {
        int first, second;
        for(int j = 0; j < edge_size; ++j) {
            std::cin >> first >> second;
            add_edge(first-start, second-start);
        }
    }

    void add_edge(int first, int second) {
        if(neighbor_set[first].find(second) != neighbor_set[first].end()) {
            return;
        }

        int new_id = edges.size();
        edge e = edge(new_id, first, second);
        edges.push_back(e);
        edge_interfaces.push_back({e.id, e.first, e.second});

        // update outgoing
        adj_list[first].push_back(new_id);
        neighbor_set[first].insert(second);
        ++deg[first];

        // update ingoiong
        adj_in_list[second].push_back(new_id);
        neighbor_in_set[second].insert(first);
        ++deg_in[second];
    }

    void delete_edge(edge& e) {
        if(e.deleted) {
            return;
        }
        e.deleted = true;
        --actual_edge_size;

        // update outgoing
        neighbor_set[e.first].erase(e.second);
        --deg[e.first];

        // update ingoing
        neighbor_in_set[e.second].erase(e.first);
        --deg_in[e.second];
    }

    void delete_edge(edge_interface& e) { 
        delete_edge(edges[e.id]); 
    }

    void delete_vertex(Graph::vertex& v) {
        if(v.deleted) {
            return;
        }
        v.deleted = true;
        --actual_vertex_size;

        // update outgoing
        for(int e_id : adj_list[v.id]) {
            edge& e = edges[e_id];
            if(!(e.deleted || (e.first != v.id && e.second != v.id))) {
                delete_edge(edges[e_id]);
            }
        }
        adj_list[v.id].clear();

        // update ingoing
        for(int e_id : adj_in_list[v.id]) {
            edge& e = edges[e_id];
            if(!e.deleted && e.second == v.id) {
                delete_edge(edges[e_id]);
            }
        }
        adj_in_list[v.id].clear();

        // the neighbor sets and degrees are updated inside delete_edge function
    }

    void delete_vertex(vertex_interface& v) { 
        delete_vertex(vertices[v]); 
    }

    bool is_deleted(const vertex_interface& v) { return vertices[v].deleted; } 

    std::set<vertex_interface> get_sorted_neighbor_set(const vertex_interface& v) { 
        return std::set<vertex_interface>(neighbor_set[v].begin(), neighbor_set[v].end()); 
    }

    std::vector<std::set<vertex_interface>> get_sorted_neighbor_set() { 
        std::vector<std::set<vertex_interface>> result(neighbor_set.size());
        // transform the edges into edge_interfaces
        std::transform(
            neighbor_set.begin(), neighbor_set.end(), result.begin(), 
            [this](std::unordered_set<vertex_interface>& s)->std::set<vertex_interface>{
                return std::set<vertex_interface>(s.begin(), s.end());
            }
        );
        return result;
    }

    DirectedGraph get_transpose() {
        DirectedGraph D_ = DirectedGraph(get_vertex_size(), get_edge_size(), get_start());
        std::vector<vertex_interface> vertices =  get_vertex_interfaces();
        for(vertex_interface v : vertices) {
            std::vector<vertex_interface> neighbors = get_neighbors(v);
            for(vertex_interface w: neighbors) {
                D_.add_edge(w, v);
            }
        }
        return D_;
    }

    DirectedGraph get_copy() {
        DirectedGraph D_(vertex_size, edge_size, start);
        for(edge_interface e : edge_interfaces) {
            if(!edges[e.id].deleted) {
                D_.add_edge(e.first, e.second);
            }
        }
        return D_;
    }

    DirectedGraph get_subgraph(std::vector<vertex_interface>& V) {
        std::vector<edge_interface> E_ = get_multiple_edge_interfaces(V);
        DirectedGraph D_(vertex_size, E_.size(), start);
        for(edge_interface e : E_) {
            D_.add_edge(e.first, e.second);
        }
        return D_;
    }

    std::vector<vertex_interface> get_neighbors(const vertex_interface& v, int skip = 0) { 
        std::vector<vertex_interface> ns = Graph::get_neighbors(v);
        if(skip >= ns.size()) {
            return {};
        }
        std::vector<vertex_interface> neighbors(ns.size());
        std::copy_if(
            ns.begin()+skip, ns.end(), 
            neighbors.begin(), 
            [&](vertex_interface v)->bool{return !is_deleted(v);}
        );
        return neighbors;
    }

    std::vector<vertex_interface> get_neighbors(const std::set<vertex_interface>& V) {
        std::vector<vertex_interface> result;
        result.reserve(V.size());
        for(vertex_interface v : V) {
            std::vector<vertex_interface> to_add = get_neighbors(v);
            std::copy_if(
                to_add.begin(), to_add.end(),
                std::inserter(result, result.end()),
                [&](vertex_interface x) { return V.find(x) == V.end(); }
            );
        }
        return result;
    }

    std::vector<vertex_interface> get_neighbors(const std::set<vertex_interface>&& V) {
        std::vector<vertex_interface> result;
        result.reserve(V.size());
        for(vertex_interface v : V) {
            std::vector<vertex_interface> to_add = get_neighbors(v);
            std::copy_if(
                to_add.begin(), to_add.end(),
                std::inserter(result, result.end()),
                [&](vertex_interface x) { return V.find(x) == V.end(); }
            );
        }
        return result;
    }

    std::vector<vertex_interface> get_neighbors(const std::vector<vertex_interface>& V) {
        std::set<vertex_interface> V_set(V.begin(), V.end());
        std::vector<vertex_interface> result;
        result.reserve(V.size());
        for(vertex_interface v : V) {
            std::vector<vertex_interface> to_add = get_neighbors(v);
            std::copy_if(
                to_add.begin(), to_add.end(),
                std::inserter(result, result.end()),
                [&](vertex_interface x) { return V_set.find(x) == V_set.end(); }
            );
        }
        return result;
    }

    std::vector<vertex_interface> get_neighbors(const std::vector<vertex_interface>&& V) {
        std::set<vertex_interface> V_set(V.begin(), V.end());
        std::vector<vertex_interface> result;
        result.reserve(V.size());
        for(vertex_interface v : V) {
            std::vector<vertex_interface> to_add = get_neighbors(v);
            std::copy_if(
                to_add.begin(), to_add.end(),
                std::inserter(result, result.end()),
                [&](vertex_interface x) { return V_set.find(x) == V_set.end(); }
            );
        }
        return result;
    }

    // just take over from Graph **********************************************
    int get_vertex_size() { return(vertex_size); }
    int get_actual_vertex_size() { return actual_vertex_size; }
    int get_edge_size() { return(edge_size); }

    std::vector<vertex_interface> get_vertex_interfaces() { return Graph::get_vertex_interfaces(); }

    std::vector<edge_interface> get_v_edge_interfaces(const vertex_interface& v) { return Graph::get_v_edge_interfaces(v); } 

    std::vector<Graph::edge> get_v_in_edges(const Graph::vertex& v) {
        // update adj_list: remove deleted and moved edges
        adj_in_list[v.id].remove_if(
            [this, v](int id) {
                edge& e = edges[id];
                return(e.deleted || (e.first != v.id && e.second != v.id));
            }
        );

        // transform ids into the edges
        std::vector<edge> v_in_edges(adj_in_list[v.id].size());
        std::transform(
            adj_in_list[v.id].begin(), adj_in_list[v.id].end(), v_in_edges.begin(), 
            [this](int e_id)->edge{return(edges[e_id]);});
        return(v_in_edges);
    }

    std::vector<edge_interface> get_v_in_edge_interfaces(const vertex_interface& v) {
        std::vector<Graph::edge> v_in_edges = get_v_in_edges(vertices[v]);
        std::vector<edge_interface> v_in_edge_interfaces(v_in_edges.size());
        // transform the edges into edge_interfaces
        std::transform(
            v_in_edges.begin(), v_in_edges.end(), v_in_edge_interfaces.begin(), 
            [this](const edge& e)->edge_interface&{return(edge_interfaces[e.id]);});
        return(v_in_edge_interfaces);
    }

    std::vector<edge_interface> get_multiple_edge_interfaces(const std::vector<vertex_interface>& V) {
        std::set<vertex_interface> V_set(V.begin(), V.end());
        std::vector<edge_interface> result;
        for(vertex_interface v : V) {
            std::vector<edge_interface> to_add = get_v_edge_interfaces(v);
            std::copy_if(
                to_add.begin(), to_add.end(),
                std::inserter(result, result.end()),
                [&](edge_interface e) { 
                    vertex_interface w = get_neighbor(v, e);
                    return V_set.find(w) != V_set.end(); 
                }
            );
        }
        return result;
    }

    vertex_interface get_neighbor(const vertex_interface& v, const edge_interface& e) { return Graph::get_neighbor(v, e); }

    // vertex decorator methods start
    int get_deg(const vertex_interface& id) { return(deg[id]); }
    bool check_deleted(vertex_interface& v) { return(vertices[v].deleted); }
    // vertex decorator methods end

    // print methods
    void print_actual_size() { Graph::print_actual_size(); }
    void print_adj() { Graph::print_adj(); }
    void print_adj_modec() { Graph::print_adj_modec(); }
    void print_edges() { Graph::print_edges(); }
    void print_me() { Graph::print_me(); }
    int get_start() { return(start); }

protected:
    // the regular ones coming from Graph are the "out" ones.
    std::vector<int> deg_in;
    std::vector<std::list<int>> adj_in_list;
    std::vector<std::set<int>> neighbor_in_set;
};