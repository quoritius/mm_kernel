#include "Graph.hpp"

Graph::Graph(int n, int m, int start, bool init) : 
    vertex_size(n), edge_size(m), start(start), 
    actual_vertex_size(n), actual_edge_size(0),
    vertices{}, edges{}, edge_interfaces{},
    deg(std::vector<int>(n, 0))
{
    vertices.reserve(vertex_size);
    adj_list.reserve(vertex_size);
    for(int i = 0; i < vertex_size; ++i) {
        vertices.push_back(vertex(i));
        adj_list.push_back(std::list<int>());
        neighbor_set.push_back({});
    }
    edges.reserve(edge_size);
    if(init) {
        int first, second;
        for(int j = 0; j < edge_size; ++j) {
            std::cin >> first >> second;
            add_edge(first-start, second-start);
        }
    }
}

Graph::Graph(int n, std::initializer_list<std::pair<int, int>> l, int start) : 
    vertex_size(n), edge_size(l.size()), start(start),
    actual_vertex_size(vertex_size), actual_edge_size(edge_size),
    vertices{}, edges{}, edge_interfaces{},
    deg(std::vector<int>(vertex_size, 0))
{
    vertices.reserve(vertex_size);
    adj_list.reserve(vertex_size); 
    for(int i = 0; i < vertex_size; ++i) {
        vertices.push_back(vertex(i));
        adj_list.push_back(std::list<int>());
        neighbor_set.push_back({});
    }
    edges.reserve(edge_size);
    int first, second;
    for(auto& edge : l) {
        add_edge(edge.first-start, edge.second-start);
    }
}

std::vector<vertex_interface> Graph::get_vertex_interfaces() {
    std::vector<vertex_interface> vertex_interfaces;
    vertex_interfaces.reserve(actual_vertex_size);
    for(vertex_interface v = 0; v < vertex_size; ++v) {
        if(!vertices[v].deleted) {
            vertex_interfaces.push_back(v);
        }
    }
    return(vertex_interfaces);
}

void Graph::clean_vertices() {
    for(vertex_interface v = 0; v < vertex_size; ++v) {
        if(deg[v] == 0) {
            delete_vertex(v);
        }
    }
}

void Graph::add_edge(int first, int second) {
    if(first == second) {
        return;
    }
    if(neighbor_set[first].find(second) != neighbor_set[first].end()) {
        return;
    }

    ++actual_edge_size;
    int new_id = edges.size();
    edge e = edge(new_id, first, second);
    edges.push_back(e);
    edge_interfaces.push_back({e.id, e.first, e.second});
    adj_list[first].push_back(new_id);
    neighbor_set[first].insert(second);
    ++deg[first];
    adj_list[second].push_back(new_id);
    neighbor_set[second].insert(first);
    ++deg[second];
}

void Graph::add_edges(vertex_interface v, std::vector<vertex_interface>& U) {
    for(auto u : U) {
        add_edge(v, u);
    }
}
void Graph::add_edges(vertex_interface v, std::vector<vertex_interface>&& U) {
    add_edges(v, U);
}
void Graph::add_edges(std::vector<vertex_interface>& U, vertex_interface v) {
    add_edges(v, U);
}
void Graph::add_edges(std::vector<vertex_interface>&& U, vertex_interface v) {
    add_edges(v, U);
}

void Graph::add_edge_pl(int first, int second) {
    if(neighbor_set[first].find(second) != neighbor_set[first].end()) {
        return;
    }
    if(first == second) {
        return;
    }

    add_edge(first, second);
    ++actual_edge_size;
}

void Graph::add_edges_pl(vertex_interface v, std::vector<vertex_interface>&& U) {
    for(auto u : U) {
        add_edge_pl(v, u);
    }
}

void Graph::delete_edge(edge& e) {
    if(e.deleted) {
        return;
    }
    e.deleted = true;
    neighbor_set[e.first].erase(e.second);
    --deg[e.first];
    neighbor_set[e.second].erase(e.first);
    --deg[e.second];
    --actual_edge_size;
}

void Graph::delete_edge(edge_interface& e) { 
    delete_edge(edges[e.id]); 
}

void Graph::move_edge(edge_interface& e, const vertex_interface& d, const vertex_interface& a) {
    edge& e_e = edges[e.id];
    edge_interface& e_i = edge_interfaces[e.id];

    if(e_e.deleted) {
        return;
    }

    vertex_interface v = get_neighbor(d, e);
    // do not create self edges
    if(v == a) {
        delete_edge(e_e);
        return;
    }
    // do not create double edge
    if(neighbor_set[v].find(a) != neighbor_set[v].end()) {
        delete_edge(e_e);
        return;
    }

    if(d == e_e.first) {
        e_e.first = a;
        e_i.first = a;
    } else {
        e_e.second = a;
        e_i.second = a;
    }
    neighbor_set[d].erase(v);
    neighbor_set[v].erase(d);
    --deg[d];
    adj_list[a].push_back(e_e.id);
    neighbor_set[a].insert(v);
    neighbor_set[v].insert(a);
    ++deg[a];
    
}

void Graph::delete_vertex(Graph::vertex& v) {
    if(v.deleted) {
        return;
    }
    v.deleted = true;
    for(int e_id : adj_list[v.id]) {
        edge& e = edges[e_id];
        if(!(e.deleted || (e.first != v.id && e.second != v.id))) {
            delete_edge(edges[e_id]);
        }
    }
    adj_list[v.id].clear();
    --actual_vertex_size;
}

void Graph::delete_vertex(vertex_interface& v) { 
    delete_vertex(vertices[v]); 
}

bool Graph::is_deleted(const vertex_interface& v) {
    return vertices[v].deleted;
}

bool Graph::is_deleted(const vertex_interface&& v) {
    return vertices[v].deleted;
}

int Graph::delete_vertices(std::vector<vertex_interface>&& vertices) {
    int num_deleted = 0;
    for(vertex_interface v : vertices) {
        if(!is_deleted(v)) {
            delete_vertex(v);
            ++num_deleted;
        }
    }
    return num_deleted;
}
int Graph::delete_vertices(std::vector<vertex_interface>& vertices) {
    int num_deleted = 0;
    for(vertex_interface v : vertices) {
        if(!is_deleted(v)) {
            delete_vertex(v);
            ++num_deleted;
        }
    }
    return num_deleted;
}

void Graph::unificate(const vertex_interface& v) {
    adj_list[v].remove_if(
        [this, v](int id)->bool {
            edge& e = edges[id];
            return(e.deleted || (e.first != v && e.second != v));
        }
    );
    // get unique edge, delete the extra ones
    adj_list[v].unique(
        [v, this](int e_id, int f_id)->bool {
            edge_interface& e = edge_interfaces[e_id];
            edge_interface& f = edge_interfaces[f_id];
            vertex_interface u = get_neighbor(v, e);
            vertex_interface w = get_neighbor(v, f);
            if(u == w) {
                delete_edge(edges[e_id]);
            } 
            return(u == w);
        }
    );
}

std::vector<Graph::edge> Graph::get_v_edges(const Graph::vertex& v) {
    // update adj_list: remove deleted and moved edges
    adj_list[v.id].remove_if(
        [this, v](int id) {
            edge& e = edges[id];
            return(e.deleted || (e.first != v.id && e.second != v.id));
        }
    );

    // transform ids into the edges
    std::vector<edge> v_edges;
    v_edges.resize(adj_list[v.id].size());
    std::transform(
        adj_list[v.id].begin(), adj_list[v.id].end(), v_edges.begin(), 
        [this](int e_id)->edge{return(edges[e_id]);});
    return(v_edges);
}

std::vector<edge_interface> Graph::get_v_edge_interfaces(const vertex_interface& v) {
    std::vector<edge> v_edges = get_v_edges(vertices[v]);
    std::vector<edge_interface> v_edge_interfaces(v_edges.size());
    // transform the edges into edge_interfaces
    std::transform(
        v_edges.begin(), v_edges.end(), v_edge_interfaces.begin(), 
        [this](const edge& e)->edge_interface&{return(edge_interfaces[e.id]);});
    return(v_edge_interfaces);
}

Graph::vertex& Graph::get_neighbor(const vertex& v, const edge& e) {
    int neighbor_id = v.id == e.first ? e.second : e.first;
    return(vertices[neighbor_id]);
}

int Graph::get_neighbor(const vertex_interface& v, const edge_interface& e) {
    return(get_neighbor(vertices[v], edges[e.id]).id);
}

std::vector<int> Graph::get_neighbors_ids(const Graph::vertex& v) {
    std::vector<edge> v_edges = get_v_edges(v);
    std::vector<int> neighbors(v_edges.size());
    std::transform(
        v_edges.begin(), v_edges.end(), neighbors.begin(),
        [&v, this](edge& e){return(get_neighbor(v, e).id);});
    return(neighbors);
}

std::vector<vertex_interface> Graph::get_neighbors(const vertex_interface& v) {
    return(get_neighbors_ids(vertices[v]));
}

std::vector<vertex_interface> Graph::get_neighbors(const std::set<vertex_interface>& V) {
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

std::vector<vertex_interface> Graph::get_neighbors(const std::set<vertex_interface>&& V) {
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

void Graph::print_actual_size() {
    std::cout << actual_vertex_size << " " << actual_edge_size << std::endl;
}

void Graph::print_actual_size(std::ofstream& out) {
    out << actual_vertex_size << " " << actual_edge_size << std::endl;
}

void Graph::print_adj() {
    print_actual_size();
    for(auto& v : vertices) {
        if(v.deleted) {
            continue;
        }
        std::cout << v.id + start << " -> ";
        for(int neighbor_id : get_neighbors_ids(v)) {
            std::cout << neighbor_id + start << " ";
        }
        std::cout << std::endl;
    }
}

void Graph::print_adj_modec() {
    for(auto& v : vertices) {
        if(v.deleted) {
            continue;
        }
        std::cout << v.id + start << "->";
        for(int neighbor_id : get_neighbors_ids(v)) {
            std::cout << neighbor_id + start << ",";
        }
        std::cout << std::endl;
    }
}

void Graph::print_edges() {
    print_actual_size();
    for(auto& e : edges) {
        if(e.deleted) {
            continue;
        }
        std::cout << e.first + start << " " << e.second + start << std::endl;
    }
}

void Graph::print_edges_relabeled() {
    int i = start;
    std::vector<int> labels(vertex_size+start, -1);
    for(auto& v : vertices) {
        if(!v.deleted) {
            labels[v.id] = i;
            ++i;
        }
    }
    print_actual_size();
    for(auto& e : edges) {
        if(!e.deleted) {
            std::cout << labels[e.first] << " " << labels[e.second] << std::endl;
        }
    }
}
void Graph::print_edges_relabeled(std::ofstream& out) {
    int i = start;
    std::vector<int> labels(vertex_size+start, -1);
    for(auto& v : vertices) {
        if(!v.deleted) {
            labels[v.id] = i;
            ++i;
        }
    }
    print_actual_size(out);
    for(auto& e : edges) {
        if(!e.deleted) {
            out << labels[e.first] << " " << labels[e.second] << std::endl;
        }
    }
}

void Graph::print_me() {
    print_actual_size();
    std::cout << "VERTICES" << std::endl;
    int cnt = 0;
    for(auto& v : vertices) {
        cnt = (cnt+1)%5;
        std::cout << v.id + start;
        if(cnt == 0) {
            std::cout << std::endl;
        } else {
            std::cout << " ";
        }
    }
    if(cnt != 0) {
        std::cout << std::endl;
    }
    std::cout << "EDGES" << std::endl;
    for(auto& e : edges) {
        if(e.deleted) {
            continue;
        }
        std::cout << e.first + start << " " << e.second + start << std::endl;
    }
}