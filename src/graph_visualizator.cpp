#include "header/graph_visualizator.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <deque>
#include <memory>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>

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

std::vector<size_t> GraphVisualizator::CreateVertexFiltration() {
  std::deque<std::shared_ptr<Vertex>> filtration(graph_.begin(), graph_.end());
  std::vector<size_t> borders{vertex_num_};

  int debug = 1000;

  std::srand(std::time(NULL));
  while (borders.back() != 3) {
    std::deque<std::shared_ptr<Vertex>> predecessor = filtration;
    predecessor.resize(borders.back());
    std::deque<std::shared_ptr<Vertex>> successor;

    size_t successor_size = 0;
    while (!predecessor.empty()) {
      auto vertex = predecessor[std::rand() % predecessor.size()];
      successor.push_front(vertex);
      ++successor_size;

      std::vector<size_t> dist = BFS(vertex);
      auto it = predecessor.begin();
      while (it != predecessor.end()) {
        auto v = *it;
  
        if (dist[v->number] <= std::pow(2, (borders.size() - 1))) {
          successor.push_back(v);
          it = predecessor.erase(it);
        } else {
          ++it;
        }
      }
    }
    
    if (successor_size < 3) {
      debug -= 1;

      if (debug == 0) {
        debug = 1000;
        borders.pop_back();
      }

      continue;
    }

    for (size_t i = 0; i < borders.back(); ++i) {
      filtration[i] = successor[i];
    }
    borders.push_back(successor_size);
  }

  graph_ = std::vector<std::shared_ptr<Vertex>>(filtration.begin(),
                                                filtration.end());
  std::reverse(borders.begin(), borders.end());
  return borders;
}