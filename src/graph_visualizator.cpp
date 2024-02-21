#include "header/graph_visualizator.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stack>

void GraphVisualizator::ReadGraph(const std::string& filename) {
  std::ifstream input(filename);

  if (!input.is_open()) {
    throw std::runtime_error("Failed to open file");
  }
  
  input >> vertex_num_ >> edge_num_;
  
  std::vector<std::vector<size_t>> graph(vertex_num_);
  for (size_t i = 0; i < edge_num_; i++) {
    size_t node_1, node_2;
    input >> node_1 >> node_2;

    graph[node_1 - 1].push_back(node_2 - 1);
    graph[node_2 - 1].push_back(node_1 - 1);
  }

  input.close();

  graph_ = graph;
}

std::vector<std::vector<size_t>>
GraphVisualizator::SearchForConnectivityComponents() {
  std::vector<std::vector<size_t>> connectivity_components;
  std::vector<bool> visited_vertices(vertex_num_, false);

  for (size_t i = 0; i < vertex_num_; ++i) {
    if (!visited_vertices[i]) {
      FindConnectivityComponent(i, &visited_vertices, 
                                &connectivity_components);
    }
  }
  
  return connectivity_components;
}

void GraphVisualizator::FindConnectivityComponent(
    const size_t& vertex_ind, 
    std::vector<bool> *visited,
    std::vector<std::vector<size_t>> *connectivity_components) {
  std::vector<size_t> connectivity_comp;
  std::stack<size_t> stack;
  stack.push(vertex_ind);

  while (!stack.empty()) {
    size_t current = stack.top();
    stack.pop();

    if ((*visited)[current]) {
      continue;
    }
    
    (*visited)[current] = true;
    connectivity_comp.push_back(current);

    for (size_t adjacent : graph_[current]) {
      if (!(*visited)[adjacent]) {
        stack.push(adjacent);
      }
    }
  }

  connectivity_components->push_back(connectivity_comp);
}