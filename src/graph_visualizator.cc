#include "header/graph_visualizator.h"

#include <algorithm>  // For std::min, std::max
#include <cmath>  // For std::sqrt, std::pow, std::round
#include <cstdint>  // For size_t, uint8_t
#include <ctime>  // For std::time
#include <fstream>  // For std::ifstream
#include <limits>  // For std::numeric_limits
#include <queue>  // For std::queue
#include <random>  // For random number generation
#include <string>  // For std::string
#include <vector>  // For std::vector

// Reads the graph data from a file and initializes 
// the graph vertices and edges.
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

  // Read the edges from the file and update the adjacency 
  // lists of the vertices
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

// Computes the global layout adjustments for vertices 
// based on the K-Centers algorithm.
std::pair<int, int> GraphVisualizator::ComputeGlobalLayout() {
  // Compute the All Pairs Shortest Path (APSP) matrix for the graph
  auto APSP = FindAllPairsShortestPath();

  // Set up a random layout for the initial positions of vertices
  SetUpRandomLayout();

  size_t k = kMinSize_;
  while (k <= vertex_num_) {
    // Find the centers using the K-Centers algorithm 
    // based on the current value of k
    auto centers = KCenters(k, APSP);
    size_t max_radius = 0;
    
    // Calculate the maximum radius among the centers
    for (size_t i : centers) {
      size_t min_dist = std::numeric_limits<size_t>::max();

      // Find the minimum distance from the current center to any other vertex
      for (size_t j = 0; j < vertex_num_; ++j) {
        if (j == i) { continue; }

        min_dist = std::min(min_dist, APSP[i][j]);
      }

      // Update max_radius if necessary
      if (min_dist > max_radius) {
        max_radius = min_dist;
      }
    }

    // Scale the max_radius by kRadius_
    max_radius *= kRadius_;

    // Perform local layout adjustments for vertices within
    // the calculated maximum radius
    ComputeLocalLayout(APSP, centers, max_radius);

    // Perturb the positions of vertices slightly
    // to avoid clustering around centers
    std::mt19937 gen(std::time(nullptr));
    std::uniform_real_distribution<> distribution(0, 1);
    for (auto vertex : graph_) {
      double min_dist = std::numeric_limits<double>::max();
      size_t min_ind = 0;

      // Find the closest center to each vertex and perturb
      // its position around that center
      for (size_t i : centers) {
        if (APSP[vertex->number][i] < min_dist) {
          min_dist = APSP[vertex->number][i];
          min_ind = i;
        }
      }

      vertex->x = graph_[min_ind]->x + distribution(gen);
      vertex->y = graph_[min_ind]->y + distribution(gen);
    }

    k *= kRatio_;  // Increase k for the next iteration
  }

  // Round vertex coordinates to integers and correct for negative values
  RoundVertexCoordinates();

  // Return the maximum x and y coordinates for drawing bmp
  std::pair<int, int> borders = CorrectCoordinates();
  return borders;
}

std::vector<std::vector<uint8_t>> GraphVisualizator::GetData() {
  // Get borders
  std::pair<int, int> borders = ComputeGlobalLayout();

  // Add padding
  borders.first += (4 - borders.first % 4) % 4;

  // Initialize data matrix with 255 (white colour)
  std::vector<std::vector<uint8_t>>
  data(borders.second, std::vector<uint8_t>(borders.first, 255));  // Trouble with type

  // Add circle around vertex position
  for (const auto& vertex : graph_) {
    // ((vertex->x + x) - vertex->x)^2 + (y - vertex->y)^2 = kCircleRadius_^2
    for (int x = -kCircleRadius_; x < kCircleRadius_; ++x) {
      double delta_y_sq = std::pow(kCircleRadius_, 2) -
                          std::pow(x, 2);
      
      double y_1 = std::round(std::sqrt(delta_y_sq) + vertex->y);
      double y_2 = std::round(-std::sqrt(delta_y_sq) + vertex->y);

      // Shade a pixel at the desired position
      data[y_1][vertex->x + x] = 0;
      data[y_2][vertex->x + x] = 0;
    }

    for (int y = -kCircleRadius_; y < kCircleRadius_; ++y) {
      double delta_x_sq = std::pow(kCircleRadius_, 2) -
                          std::pow(y, 2);
      
      double x_1 = std::round(std::sqrt(delta_x_sq) + vertex->x);
      double x_2 = std::round(-std::sqrt(delta_x_sq) + vertex->x);

      // Shade a pixel at the desired position
      data[vertex->y + y][x_1] = 0;
      data[vertex->y + y][x_2] = 0;
    }
  }

  // Add edges between vertices
  for (const auto& vertex : graph_) {
    for (const auto& neighbour_weak : vertex->neighbours) {
      auto neighbour = neighbour_weak.lock();

      if (!neighbour) {
        throw std::logic_error("Neighbour is a null pointer");
      }
      
      // Draw an edge
      for (size_t x = vertex->x; x < neighbour->x; ++x) {
        double y = std::round((x - vertex->x) * (neighbour->y - vertex->y) / 
                   (neighbour->x - vertex->x) + vertex->y);
        
        // Shade a pixel at the desired position
        if (std::pow(y - neighbour->y, 2) + 
            std::pow(x - neighbour->x, 2) < std::pow(kRadius_, 2)) {
          continue;
        }

        if (std::pow(y - vertex->y, 2) + 
            std::pow(x - vertex->x, 2) < std::pow(kRadius_, 2)) {
          continue;
        }

        data[y][x] = 0;
      }

      for (size_t y = vertex->y; y < neighbour->y; ++y) {
        double x = std::round((y - vertex->y) * (neighbour->x - vertex->x) / 
                   (neighbour->y - vertex->y) + vertex->x);
        
        // Shade a pixel at the desired position
        if (std::pow(y - neighbour->y, 2) + 
            std::pow(x - neighbour->x, 2) < std::pow(kRadius_, 2)) {
          continue;
        }

        if (std::pow(y - vertex->y, 2) + 
            std::pow(x - vertex->x, 2) < std::pow(kRadius_, 2)) {
          continue;
        }

        data[y][x] = 0;
      }
    }
  }
  return data;
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
  std::uniform_real_distribution<> distribution(0, vertex_num_ * kEdgeLen_);

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

void GraphVisualizator::ComputeLocalLayout(
    const std::vector<std::vector<size_t>>& distances,
    const std::vector<size_t>& centers,
    const size_t& k) {
  // Perform local layout adjustments for a specified number of iterations
  for (size_t i = 1; i < kIterations_ * vertex_num_; ++i) {
    size_t ind = ChooseVertex(distances, centers, k);

    // Find displacement
    std::pair<double, double> displacement = FindKSmallDelta(distances, 
                                                             ind, k);
    
    // Update vertex
    graph_[ind]->x += displacement.first;
    graph_[ind]->y += displacement.second;
  }
}

size_t GraphVisualizator::ChooseVertex(
    const std::vector<std::vector<size_t>>& distances,
    const std::vector<size_t>& centers,
    const size_t& k) {
  double max_big_delta = 0.0;
  size_t max_ind = 0;

  for (size_t j : centers) {
    double big_delta = FindKDelta(distances, j, k);

    if (big_delta >= max_big_delta) {
      max_big_delta = big_delta;
      max_ind = j;
    }
  }

  return max_ind;
}

// Computes the derivatives of the energy function and uses them to determine
// the change in position (delta) along the x and y axes.
std::pair<double, double> GraphVisualizator::FindKSmallDelta(
    const std::vector<std::vector<size_t>>& distances,
    const size_t& vertex_ind,
    const size_t& k) {
  // Compute partial derivatives of the energy function
  std::vector<double> partial_derivatives = FindKPartialDerivatives(
                                                distances, vertex_ind, k);

  // Compute first derivatives of the energy function
  std::pair<double, double> first_derivatives = FindKEnergyDerivative(
                                                    distances, vertex_ind, k);
  
  // Extract coefficients for the linear system
  double a_1 = partial_derivatives[0];
  double b_1 = partial_derivatives[1];
  double c_1 = -first_derivatives.first;

  double a_2 = partial_derivatives[1];
  double b_2 = partial_derivatives[2];
  double c_2 = -first_derivatives.second;

  // Solve the linear system to find the delta values
  // a_1 * delta_x + b_1 * delta_y = c_1
  // a_2 * delta_x + b_2 * delta_y = c_2
  double delta_x = 0, delta_y = 0;
  if (a_1 != 0) { // Normalize first coeff
    b_1 /= a_1;
    c_1 /= a_1;
    a_1 = 1.0;

    if (a_2 != 0) { // Normalize first coeff
      b_2 /= a_2;
      c_2 /= a_2;
      a_2 = 1.0;

      // Subtract the 2nd equation from the 1st equation
      a_1 = 0.0;
      b_1 -= b_2;
      c_1 -= c_2;

      if (b_1 != 0) {
        delta_y = c_1 / b_1;
        delta_x = c_2 - delta_y * b_2;
      }

    } else if (b_2 != 0) {
      delta_y = c_2 / b_2;
      delta_x = (c_1 - delta_y * b_1) / a_1;
    }
  
  } else if (b_1 != 0) {
    delta_y = c_1 / b_1;

    if (a_2 != 0) {
      delta_x = (c_2 - delta_y * b_2) / a_2;
    }
  }
  
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

    double k_ij = 1.0 / distances[vertex_ind][i];
    double dist = FindEuclideanDistance(vertex_ind, i);

    Vertex v_ind = *graph_[vertex_ind];
    Vertex v_i = *graph_[i];

    derivative_x_x += k_ij * (1.0 - kEdgeLen_ * 
                              distances[vertex_ind][i] * 
                              std::pow(v_ind.y - v_i.y, 2) /
                              std::pow(dist, 3));
    
    derivative_x_y += k_ij * kEdgeLen_ * 
                      distances[vertex_ind][i] * 
                      (v_ind.y - v_i.y) * (v_ind.x - v_i.x) /
                      std::pow(dist, 3);

    derivative_y_y += k_ij * (1.0 - kEdgeLen_ * 
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
    double k_ij = 1.0 / distances[vertex_ind][i];

    double dist = FindEuclideanDistance(vertex_ind, i);
    
    Vertex v_ind = *graph_[vertex_ind];
    Vertex v_i = *graph_[i];

    x_k_energy_derivative += k_ij * (v_ind.x - v_i.x) * 
                             (1.0 - kEdgeLen_ * distances[vertex_ind][i] / dist);
    
    y_k_energy_derivative += k_ij * (v_ind.y - v_i.y) * 
                             (1.0 - kEdgeLen_ * distances[vertex_ind][i] / dist);
  }

  return {x_k_energy_derivative, y_k_energy_derivative};
}

std::vector<size_t> GraphVisualizator::FindKNeighbourhood(
    const std::vector<std::vector<size_t>>& distances,
    const size_t& vertex_ind,
    const size_t& k) {
  std::vector<size_t> k_neighbourhood;

  for (size_t i = 0; i < vertex_num_; ++i) {
    // Check if the distance from the current vertex to the 
    // target vertex is within the threshold
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

void GraphVisualizator::RoundVertexCoordinates() {
  for (size_t i = 0; i < vertex_num_; ++i) {
    graph_[i]->x = std::round(graph_[i]->x);
    graph_[i]->y = std::round(graph_[i]->y);
  }
}

std::pair<int, int> GraphVisualizator::CorrectCoordinates() {
  int min_x = std::numeric_limits<int>::max();
  int min_y = std::numeric_limits<int>::max();
  
  int max_x = std::numeric_limits<int>::min();
  int max_y = std::numeric_limits<int>::min();
  
  // Iterate over vertices to find minimum and maximum coordinates
  for (const auto& vertex : graph_) {
    min_x = std::min(min_x, static_cast<int>(vertex->x));
    min_y = std::min(min_y, static_cast<int>(vertex->y));

    max_x = std::max(max_x, static_cast<int>(vertex->x));
    max_y = std::max(max_y, static_cast<int>(vertex->y));
  }

  // Calculate displacements to ensure all coordinates are positive
  int disp_x = kEdgeLen_ - min_x;
  int disp_y = kEdgeLen_ - min_y;

  // Update vertex coordinates to ensure positivity
  for (size_t i = 0; i < vertex_num_; ++i) {
    graph_[i]->x += disp_x;
    graph_[i]->y += disp_y;
  }

  // Return corrected coordinates with additional margin
  return {max_x + disp_x + kEdgeLen_, max_y + disp_y + kEdgeLen_};
}