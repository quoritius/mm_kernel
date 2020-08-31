#include <iostream>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <chrono>
#include "Graph.hpp"
#include "DirectedGraph.hpp"
#include "Kernelizer.hpp"
#include "helping_algos.hpp"

std::vector<int> full_test_exhaustive(Graph& G, bool deg012, bool crown, bool relcrown) {
    std::vector<int> solutions;
    int tmp = 0;
    solve_isolated(G);

    // exhaustively remove deg0+deg1 and deg2
    if(deg012) {
        tmp = 0;
        solutions.push_back(0);
        do {
            solutions.back() += tmp;
            tmp = solve_externals(G) + solve_redexes(G);
        } while(solutions.back() + tmp != solutions.back());
    }

    // exhaustively remove crowns
    if(crown) {
        solutions.push_back(0);
        solutions.back() += solve_crowns(G); 
    }

    // exhaustively remove crowns and deg0 and relcrowns
    if(relcrown) {
        tmp = 0;
        solutions.push_back(0);
        do {
            solutions.back() += tmp;
            tmp = solve_crowns(G) + solve_isolated(G) + solve_relaxed_crowns(G);
        } while(tmp > 10 && tmp >= G.get_actual_edge_size()/50);
        solutions.back() += tmp;
    }

    return solutions;
}


int main(int argc, char* argv[]) {
    // check usage validity
    if(argc != 5) {
        std::cerr << "Use mode and 3 filepaths (input, solution decrease, output) as arguments please." << std::endl;
        return -1;
    }
    int mode = atoi(argv[1]);
    if(mode == 0 && mode > 5) {
        std::cerr << "Wrong mode given as argument." << std::endl;
        return -1;
    }

    // read graph
    freopen(argv[2], "r", stdin);
    std::ofstream sol_file, output_graph;
    int n, m, start;
    std::cin >> start >> n >> m;
    Graph G = Graph(n, m, start);
    fclose(stdin);

    // get mode
    bool deg012, crown, relcrown;
    deg012 = mode < 4;
    crown = mode > 1;
    relcrown = mode == 3 || mode == 5;

    // kernelize
    auto start_time = std::chrono::steady_clock::now();
    std::vector<int> solutions = full_test_exhaustive(G, deg012, crown, relcrown);
    auto end_time = std::chrono::steady_clock::now();
    std::cout << "Kernelization needed time: " << 
        std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << 
        " ms" << std::endl;

    // save solution decrease
    sol_file.open(argv[3]);
    for(auto sol : solutions) {
        sol_file << sol << "; ";
    }
    sol_file << std::accumulate(solutions.begin(), solutions.end(), 0);
    sol_file.close();

    // save the kernel
    output_graph.open(argv[4]);
    output_graph << start << " ";
    G.print_edges_relabeled(output_graph);
    output_graph.close();
    return 0;
}