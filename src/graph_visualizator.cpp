#include "header/graph_visualizator.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stack>
#include <utility>

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

std::vector<std::pair<size_t, size_t>> GraphVisualizator::SearchForBridges() {
  std::vector<bool> visited_vertices(vertex_num_, false);
  std::vector<size_t> time_in(vertex_num_);
  std::vector<size_t> min_time_in_reachable(vertex_num_);
  std::vector<std::pair<size_t, size_t>> bridges;

  for (size_t i = 0; i < vertex_num_; ++i) {
    if (!visited_vertices[i]) {
      FindBridge(i, -1, &visited_vertices, &time_in, 
                 &min_time_in_reachable, &bridges);
    }
  }

  return bridges;
}

void GraphVisualizator::FindBridge(
    const size_t& vertex_ind,
    const size_t& parent_ind,
    std::vector<bool> *visited,
    std::vector<size_t> *time_in,
    std::vector<size_t> *min_time_in_reachable,
    std::vector<std::pair<size_t, size_t>> *bridges) {
  static size_t timer = 0;

  (*visited)[vertex_ind] = true;
  (*time_in)[vertex_ind] = (*min_time_in_reachable)[vertex_ind] = timer++;
  for (size_t i = 0; i < graph_[vertex_ind].size(); ++i) {
    size_t to = graph_[vertex_ind][i];
    if (to == parent_ind) {
      continue;
    }

    if ((*visited)[to]) {
      (*min_time_in_reachable)[vertex_ind] = std::min(
          (*min_time_in_reachable)[vertex_ind],
          (*time_in)[to]);
    } else {
      FindBridge(to, vertex_ind, visited, time_in, 
                 min_time_in_reachable, bridges);
      
      (*min_time_in_reachable)[vertex_ind] = std::min(
        (*min_time_in_reachable)[vertex_ind],
        (*min_time_in_reachable)[to]);
      
      if ((*min_time_in_reachable)[to] > (*time_in)[vertex_ind]) {
        bridges->push_back(std::make_pair(vertex_ind, to));
      }
    }
  }
}