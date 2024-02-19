#include "header/graph_visualizator.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

std::vector<std::vector<int>>
GraphVisualizator::ReadGraph(const std::string& filename) {
  std::ifstream input(filename);

  if (!input) {
    std::cerr << "File wasn't opened";
    return;
  }
  
  input >> vertex_num_ >> edge_num_;
  
  std::vector<std::vector<int>> graph(vertex_num_);
  for (size_t i = 0; i < edge_num_; i++) {
    size_t node_1, node_2;
    input >> node_1 >> node_2;

    graph[node_1 - 1].push_back(node_2);
    graph[node_2 - 1].push_back(node_1);
  }

  input.close();

  return graph;
}