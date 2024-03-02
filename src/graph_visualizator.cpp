#include "header/graph_visualizator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
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
  
  // Create and initialize vertices for the graph
  std::vector<std::shared_ptr<Vertex>> graph(vertex_num_, nullptr);
  for (size_t i = 0; i < vertex_num_; ++i) {
    graph[i] = std::make_shared<Vertex>(i);
  }

  // Read the edges from the file and update the adjacency lists of the vertices
  for (size_t i = 0; i < edge_num_; ++i) {
    size_t ind_1, ind_2;
    input >> ind_1 >> ind_2;

    // Check if the indices are valid
    if (ind_1 < 1 || ind_1 > vertex_num_ || ind_2 < 1 || ind_2 > vertex_num_) {
      throw std::runtime_error("Invalid vertex index in edge definition");
    }

    graph[ind_1 - 1]->neighbours.push_back(graph[ind_2 - 1]);
    graph[ind_2 - 1]->neighbours.push_back(graph[ind_1 - 1]);
  }

  input.close();

  graph_ = graph;
}

std::vector<size_t> GraphVisualizator::BFS(const std::shared_ptr<Vertex>& root) {
  // Since the distance in the graph between 2 vertices cannot be greater
  // than vertex_num_ - 1, the value of vertex_num_ can be used
  // to initialise unreachability.
  std::vector<size_t> dist(vertex_num_, vertex_num_);
  dist[root->number] = 0;

  std::queue<std::shared_ptr<Vertex>> queue;
  queue.push(root);

  while (!queue.empty()) {
    std::shared_ptr<Vertex> vertex = queue.front();
    queue.pop();

    // Traverse the neighbors of the current vertex
    for (const auto& neighbour_weak : vertex->neighbours) {
      auto neighbour = neighbour_weak.lock();

      if (!neighbour) {
        throw std::logic_error("Neighbour is a null pointer");
      }

      // Update distance if shorter path found
      if (dist[neighbour->number] > dist[vertex->number] + 1) {
        dist[neighbour->number] = dist[vertex->number] + 1;
        queue.push(neighbour);
      }
    }
  }

  return dist;
}

std::vector<std::vector<size_t>> GraphVisualizator::FindAllPairsShortestPath() {
  std::vector<std::vector<size_t>> APSP(vertex_num_);

  // Compute shortest paths from each vertex to all other vertices using BFS
  for (size_t i = 0; i < vertex_num_; ++i) {
    APSP[i] = BFS(graph_[i]);
  }

  return APSP;
}

void GraphVisualizator::SetUpRandomLayout() {
  std::mt19937 gen(std::time(nullptr));
  std::uniform_int_distribution<> distribution(0, vertex_num_ * kEdgeLen_);

  for (auto vertex : graph_) {
    vertex->x = distribution(gen);
    vertex->y = distribution(gen);
  }
}

std::vector<size_t>
GraphVisualizator::KCenters(const size_t& k,
                            const std::vector<std::vector<size_t>>& distances) {
  std::vector<bool> in_centers(vertex_num_, false);
  std::vector<size_t> centers;

  // Randomly choose the first center
  std::mt19937 gen(std::time(nullptr));
  std::uniform_int_distribution<> distribution(0, vertex_num_ - 1);

  size_t rand_ind = distribution(gen);
  centers.push_back(rand_ind);
  in_centers[rand_ind] = true;

  // Iteratively select the remaining centers
  for (size_t i = 1; i < k; ++i) {
    // Find the farthest vertex from the already selected centers
    size_t vertex_ind = FindFarthestVertex(distances, in_centers);
    centers.push_back(vertex_ind);
    in_centers[vertex_ind] = true;
  }

  return centers;
}

size_t GraphVisualizator::FindFarthestVertex(
    const std::vector<std::vector<size_t>>& distances,
    const std::vector<bool>& in_centers) {
  size_t farthest_vert_ind = 0;
  size_t max_dist = 0;

  for (size_t i = 0; i < vertex_num_; ++i) {
    // Skip vertices that are already selected as centers
    if (in_centers[i]) { continue; }
    
    // Iterate over all selected centers to find the minimum distance to centers
    size_t min_dist_to_centers = SIZE_MAX;
    for (size_t j = 0; j < vertex_num_; ++j) {
      if (!in_centers[j]) { continue; }

      min_dist_to_centers = std::min(min_dist_to_centers, distances[i][j]);
    }

    // Update the farthest vertex if it has a greater minimum distance to centers
    if (min_dist_to_centers >= max_dist) {
      max_dist = min_dist_to_centers;
      farthest_vert_ind = i;
    }
  }

  return farthest_vert_ind;
}

// Computes the derivatives of the energy function and uses them to determine
// the change in position (delta) along the x and y axes.
std::pair<double, double> GraphVisualizator::FindKSmallDelta(
    const std::vector<std::vector<size_t>>& distances,
    const size_t& vertex_ind,
    const size_t& k) {
  std::vector<double> partial_derivatives = FindKPartialDerivatives(
                                                distances, vertex_ind, k);

  std::pair<double, double> first_derivatives = FindKEnergyDerivative(
                                                    distances, vertex_ind, k);
  
  double delta_x = (-first_derivatives.first * partial_derivatives[2] + 
                    partial_derivatives[1] * -first_derivatives.second) /
                   (std::pow(partial_derivatives[1], 2) + 
                   partial_derivatives[2] * partial_derivatives[0]);
  
  double delta_y = (partial_derivatives[1] * -first_derivatives.first +
                    first_derivatives.second * partial_derivatives[0]) /
                   (std::pow(partial_derivatives[1], 2) + 
                   partial_derivatives[2] * partial_derivatives[0]);
  
  return {delta_x, delta_y};
}

std::vector<double> GraphVisualizator::FindKPartialDerivatives(
    const std::vector<std::vector<size_t>>& distances,
    const size_t& vertex_ind,
    const size_t& k) {
  std::vector<size_t> k_neighbourhood = FindKNeighbourhood(distances, 
                                                           vertex_ind, k);
  double derivative_x_x = 0.0;
  double derivative_x_y = 0.0;
  double derivative_y_y = 0.0;

  for (size_t i : k_neighbourhood) {
    if (i == vertex_ind) { continue; }

    double k_ij = 1 / distances[vertex_ind][i];
    double dist = FindEuclideanDistance(vertex_ind, i);

    Vertex v_ind = *graph_[vertex_ind];
    Vertex v_i = *graph_[i];

    derivative_x_x += 2 * k_ij * (1 - kEdgeLen_ * 
                                      distances[vertex_ind][i] * 
                                      std::pow(v_ind.y - v_i.y, 2) /
                                      std::pow(dist, 3));
    
    derivative_x_y += 2 * k_ij * kEdgeLen_ * 
                      distances[vertex_ind][i] * 
                      (v_ind.y - v_i.y) * (v_ind.x - v_i.x) /
                      std::pow(dist, 3);

    derivative_y_y += 2 * k_ij * (1 - kEdgeLen_ * 
                                      distances[vertex_ind][i] * 
                                      std::pow(v_ind.x - v_i.x, 2) /
                                      std::pow(dist, 3));
  }

  return {derivative_x_x, derivative_x_y, derivative_y_y};
}

double GraphVisualizator::FindKDelta(
    const std::vector<std::vector<size_t>>& distances,
    const size_t& vertex_ind, 
    const size_t& k) {
  auto derivatives = FindKEnergyDerivative(distances, vertex_ind, k);

  return std::sqrt(std::pow(derivatives.first, 2) +
                   std::pow(derivatives.second, 2));
}

std::pair<double, double> GraphVisualizator::FindKEnergyDerivative(
      const std::vector<std::vector<size_t>>& distances,
      const size_t& vertex_ind,
      const size_t& k) {
  // Find the k-neighbourhood of the vertex
  std::vector<size_t> k_neighbourhood = FindKNeighbourhood(distances, 
                                                           vertex_ind, k);
  double x_k_energy_derivative = 0.0;
  double y_k_energy_derivative = 0.0;

  // Iterate over the vertices in the k-neighbourhood to compute derivatives
  for (size_t i : k_neighbourhood) {
    if (i == vertex_ind) { continue; }

    // weighting constant that can be either:
    // (1 / distance[u][v]) or (1 / distance[u][v]^2 )
    double k_ij = 1 / distances[vertex_ind][i];

    double dist = FindEuclideanDistance(vertex_ind, i);
    
    Vertex v_ind = *graph_[vertex_ind];
    Vertex v_i = *graph_[i];

    x_k_energy_derivative += 2 * k_ij * (v_ind.x - v_i.x) * 
                             (1 - kEdgeLen_ * distances[vertex_ind][i] / dist);
    
    y_k_energy_derivative += 2 * k_ij * (v_ind.y - v_i.y) * 
                             (1 - kEdgeLen_ * distances[vertex_ind][i] / dist);
  }

  return {x_k_energy_derivative, y_k_energy_derivative};
}

std::vector<size_t> GraphVisualizator::FindKNeighbourhood(
    const std::vector<std::vector<size_t>>& distances,
    const size_t& vertex_ind,
    const size_t& k) {
  std::vector<size_t> k_neighbourhood;

  for (size_t i = 0; i < vertex_num_; ++i) {
    if (distances[vertex_ind][i] < k) {
      k_neighbourhood.push_back(i);
    }
  }

  return k_neighbourhood;
}

double GraphVisualizator::FindEuclideanDistance(const size_t& a, 
                                                const size_t& b) {
  Vertex v_1 = *graph_[a];
  Vertex v_2 = *graph_[b];
  
  return std::sqrt(std::pow(v_1.x - v_2.x, 2) + std::pow(v_1.y - v_2.y, 2));
}
