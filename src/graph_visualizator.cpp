#include "header/graph_visualizator.h"

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <string>
#include <vector>

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

    for (const auto& neighbour : vertex->neighbours) {
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

  size_t repeat_before_reset = kRepeatBeforeReset_;

  std::mt19937 gen(std::time(nullptr));
  while (borders.back() != kNeighbourNumber_) {
    std::deque<std::shared_ptr<Vertex>> predecessor = filtration;
    predecessor.resize(borders.back());

    std::deque<std::shared_ptr<Vertex>> successor;
    size_t successor_size = 0;

    while (!predecessor.empty()) {
      std::uniform_int_distribution<> distribution(0, predecessor.size() - 1);
      size_t rand_ind = distribution(gen);
      auto vertex = predecessor[rand_ind];
      
      successor.push_front(vertex);
      predecessor.erase(predecessor.begin() + rand_ind);
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
    
    if (successor_size < kNeighbourNumber_) {
      if (--repeat_before_reset == kResetThreshhold_) {
        repeat_before_reset = kRepeatBeforeReset_;
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

void GraphVisualizator::PlaceCoreVertices() {
  std::vector<size_t> triangle_dist(kNeighbourNumber_, 0); // 1-2, 2-3, 3-1

  std::vector<size_t> dist = BFS(graph_[0]);
  triangle_dist[0] = dist[graph_[1]->number];
  triangle_dist[2] = dist[graph_[2]->number];

  dist = BFS(graph_[1]);
  triangle_dist[1] = dist[graph_[2]->number];

  graph_[1]->x = dist[0];

  const double kHalfPerimeter = std::reduce(triangle_dist.begin(), 
                                            triangle_dist.end(), 0.0) / 2;
  const double kArea = std::sqrt(kHalfPerimeter * 
                                 (kHalfPerimeter - triangle_dist[0]) *
                                 (kHalfPerimeter - triangle_dist[1]) * 
                                 (kHalfPerimeter - triangle_dist[2]));
  const double kHeight = (kArea / triangle_dist[0]) * 2;
  const double kXCord = std::sqrt(std::pow(triangle_dist[2], 2) - 
                                  std::pow(kHeight, 2));

  graph_[2]->x = kXCord;
  graph_[2]->y = kHeight;
}