#pragma once
#include <vector>
#include <set>
#include <map>


struct Dominator_Solver {
    // class params
    int r, n;
    std::vector<std::set<int>> succ;

    // n-size vectors
    std::vector<int> dom;
    std::vector<int> parent;
    std::vector<int> ancestor;
    std::vector<int> vertex;

    // n+1-size vectors
    std::vector<int> label;
    std::vector<int> semi;

    // sets of (actual) size n
    std::vector<std::set<int>> pred;
    std::vector<std::set<int>> bucket;

    int u, v, x;

    Dominator_Solver(int r, int n, std::vector<std::set<int>> &succ) : r(r), n(n), succ(succ) {
        dom.resize(n, -1);
        parent.resize(n, -1);
        ancestor.resize(n, -1);
        vertex.resize(n, -1);

        label.resize(n+1);
        semi.resize(n+1, -1);

        pred.resize(n);
        bucket.resize(n);
    }
    Dominator_Solver(int r, int n, std::vector<std::set<int>> &&succ) : r(r), n(n), succ(succ) {
        dom.resize(n, -1);
        parent.resize(n, -1);
        ancestor.resize(n, -1);
        vertex.resize(n, -1);

        label.resize(n+1);
        semi.resize(n+1, -1);

        pred.resize(n);
        bucket.resize(n);
    }
    void dfs(int v) {
        ++n;
        semi[v] = n;
        vertex[n] = label[v] = v;
        ancestor[v] = -1;

        for(auto w : succ[v]) {
            if(semi[w] == -1) {
                parent[w] = v;
                dfs(w);
            }
            pred[w].insert(v);
        }
    }

    void compress(int v) {
        if(ancestor[ancestor[v]] != -1) {
            compress(ancestor[v]);
            if(semi[label[ancestor[v]]] < semi[label[v]]) {
                label[v] = label[ancestor[v]];
            }
            ancestor[v] = ancestor[ancestor[v]];
        }
    }

    int eval(int v) {
        if(ancestor[v] == -1) {
            return v;
        } else {
            compress(v);
            return label[v];
        }
    }

    void link(int v, int w) {
        ancestor[w] = v;
    }

    std::vector<int> get_doms() {
        int w = -1;
        // step 1
        for(int v = 0; v < n; ++v) {
            pred[v].clear(); 
            bucket[v].clear();
            semi[v] = -1;
        }
        n = -1;
        dfs(r);

        for(int i = n; i > 0; --i) {
            w = vertex[i];
            // step 2
            for(auto v : pred[w]) {
                u = eval(v);
                if(semi[u] < semi[w]) {
                    semi[w] = semi[u];
                }
            }
            bucket[vertex[semi[w]]].insert(w);
            link(parent[w], w);
            
            // step 3
            for(auto& v : bucket[parent[w]]) {
                bucket[parent[w]].erase(v);
                u = eval(v);
                dom[v] = (semi[u] < semi[v])? u : parent[w];
            }
        }

        // step 4
        for(int i = 1; i < n; ++i) {
            u = vertex[i];
            if(dom[w] != vertex[semi[w]]) {
                dom[w] = dom[dom[w]];
            }
        }
        dom[r] = -1;

        return dom;
    }
};