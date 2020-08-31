#pragma once
#include <list>
#include <vector>
#include <set>
#include <utility>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <fstream>

/*
 * The user of this class gets only int (for vertex id) and 
 * {int, vertex_interface, vertex_interface} (for edge id between two vertices) 
 * as the interface.
 * The user should never get or use vertices or edges.
 */

typedef int vertex_interface;
struct edge_interface {
    edge_interface(int id, vertex_interface first, vertex_interface second) : 
        id(id), first(first), second(second) {};
    edge_interface() : id(-1), first(-1), second(-1) {};
    int id;
    vertex_interface first;
    vertex_interface second;
};

class Graph {
protected: 
    struct vertex {
        vertex() : id(-1) { deleted = true; }
        vertex(int id) : id(id) {}

        int id;
        bool deleted = false;
    };

    struct edge {
        edge() : id(-1), first(-1), second(-1) { deleted = true; }
        edge(int id, int v, int u) : id(id), first(v), second(u) {}

        int id;
        int first;
        int second;
        bool deleted = false;
    };
public:
    Graph() {};
    Graph(int n, int m, int start = 1, bool init = true);
    Graph(int n, std::initializer_list<std::pair<int, int>> l, int start);

    int get_vertex_size() { return(vertex_size); }
    int get_actual_vertex_size() { return actual_vertex_size; }
    int get_edge_size() { return(edge_size); }
    int get_actual_edge_size() { return actual_edge_size; }

    std::vector<vertex_interface> get_vertex_interfaces();

    void add_edge(int first, int second);
    void add_edge_pl(int first, int second);
    void add_edges(vertex_interface v, std::vector<vertex_interface>& U);
    void add_edges(vertex_interface v, std::vector<vertex_interface>&& U);
    void add_edges(std::vector<vertex_interface>& U, vertex_interface v);
    void add_edges(std::vector<vertex_interface>&& U, vertex_interface v);
    void add_edges_pl(vertex_interface v, std::vector<vertex_interface>&& U);
    void delete_edge(edge_interface& e);
    void move_edge(edge_interface& e, const vertex_interface& d, const vertex_interface& a);

    void delete_vertex(vertex_interface& v);
    int delete_vertices(std::vector<vertex_interface>& vertices);
    int delete_vertices(std::vector<vertex_interface>&& vertices);
    bool is_deleted(const vertex_interface& v);
    bool is_deleted(const vertex_interface&& v);
    /* time complexity: O(deg[v]) but does not do what it is intended to do.*/
    void unificate(const vertex_interface& v);
    void clean_vertices();

    std::vector<edge_interface> get_v_edge_interfaces(const vertex_interface& v);
    std::vector<vertex_interface> get_neighbors(const std::set<vertex_interface>& V);
    std::vector<vertex_interface> get_neighbors(const std::set<vertex_interface>&& V);

    vertex_interface get_neighbor(const vertex_interface&, const edge_interface& e);
    // ids instead of vectors as we don't want to give (use) any pointers.
    std::vector<vertex_interface> get_neighbors(const vertex_interface& v);

    // vertex decorator methods start
    int get_deg(const vertex_interface& id) { return(deg[id]); }
    bool check_deleted(vertex_interface& v) { return(vertices[v].deleted); }
    // vertex decorator methods end

    // print methods
    void print_actual_size();
    void print_actual_size(std::ofstream& out);
    void print_adj();
    void print_adj_modec();
    void print_edges();
    void print_edges_relabeled();
    void print_edges_relabeled(std::ofstream& out);
    void print_me();
    int get_start() { return(start); }
protected:
    int start;
    int vertex_size;
    int actual_vertex_size;
    int edge_size;
    int actual_edge_size;
    std::vector<vertex> vertices;
    std::vector<edge> edges;
    std::vector<edge_interface> edge_interfaces;

    // vertex decorators start
    std::vector<int> deg;
    std::vector<std::list<int>> adj_list;
    std::vector<std::unordered_set<int>> neighbor_set;
    // vertex decorators end

    vertex& get_neighbor(const vertex& v, const edge& e);
    std::vector<int> get_neighbors_ids(const vertex& v);
    void delete_vertex(vertex& v);

    void delete_edge(edge& e);
    std::vector<edge> get_v_edges(const vertex& v);
};