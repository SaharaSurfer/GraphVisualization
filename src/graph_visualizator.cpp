#include "header/graph_visualizator.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <memory>

void GraphVisualizator::ReadGraph(const std::string& filename) {
  std::ifstream input(filename);

  if (!input.is_open()) {
    throw std::runtime_error("Failed to open file");
  }
  
  input >> vertex_num_ >> edge_num_;
  
  std::vector<std::shared_ptr<Vertex>> graph(vertex_num_, nullptr);
  for (size_t i = 0; i < vertex_num_; ++i) {
    graph[i] = std::make_shared<Vertex>(i);
  }

  for (size_t i = 0; i < edge_num_; ++i) {
    size_t ind_1, ind_2;
    input >> ind_1 >> ind_2;

    graph[ind_1 - 1]->neighbours.push_back(graph[ind_2 - 1]);
    graph[ind_2 - 1]->neighbours.push_back(graph[ind_1 - 1]);
  }

  input.close();

  graph_ = graph;
}

std::vector<size_t> GraphVisualizator::BFS(std::shared_ptr<Vertex> root) {
  std::vector<size_t> dist(vertex_num_, vertex_num_);
  dist[root->number] = 0;
  
  std::queue<std::shared_ptr<Vertex>> queue;
  queue.push(root);

  while (!queue.empty()) {
    std::shared_ptr<Vertex> vertex = queue.front();
    queue.pop();

    for (auto neighbour : vertex->neighbours) {
      if (dist[neighbour->number] > dist[vertex->number] + 1) {
        dist[neighbour->number] = dist[vertex->number] + 1;
        queue.push(neighbour);
      }
    }
  }

  return dist;
}