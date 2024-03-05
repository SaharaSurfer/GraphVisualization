#include "header/graph_visualizator.h"

#include <algorithm>  // For std::min, std::max, std::transform, std::max_element
#include <cmath>  // For std::sqrt, std::pow, std::round
#include <cstdint>  // For size_t, uint8_t
#include <fstream>  // For std::ifstream
#include <iterator> // For std::distance
#include <limits>  // For std::numeric_limits
#include <queue>  // For std::queue, std::priority_queue
#include <random>  // For random number generation
#include <stdexcept> // For std::runtime_error, std::logic_error
#include <string>  // For std::string
#include <vector>  // For std::vector
#include <memory> // For std::shared_ptr

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
      for (size_t j : centers) {
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
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> distribution(0.1, 1);
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
  data(borders.second, std::vector<uint8_t>(borders.first, 255));

  // Add circle around vertex position
  for (const auto& vertex : graph_) {
    for (int y = -kCircleRadius_; y <= kCircleRadius_; ++y) {
      for (int x = -kCircleRadius_; x <= kCircleRadius_; ++x) {
        if (x * x + y * y <= kCircleRadius_ * kCircleRadius_) {
          int posX = static_cast<int>(std::round(vertex->x)) + x;
          int posY = static_cast<int>(std::round(vertex->y)) + y;

          if (posX >= 0 && posY >= 0 && posX < borders.first && posY < borders.second) {
            data[posY][posX] = 0; // Set pixel to black
          }
        }
      }
    }
  }

  // Add edges between vertices
  for (const auto& vertex : graph_) {
    for (const auto& neighbour_weak : vertex->neighbours) {
      auto neighbour = neighbour_weak.lock();

      if (!neighbour) {
        throw std::logic_error("Neighbour is a null pointer");
      }

      // Bresenham's line algorithm to draw a line between two points
      int x0 = static_cast<int>(std::round(vertex->x));
      int y0 = static_cast<int>(std::round(vertex->y));
      int x1 = static_cast<int>(std::round(neighbour->x));
      int y1 = static_cast<int>(std::round(neighbour->y));

      int dx = std::abs(x1 - x0);
      int dy = std::abs(y1 - y0);
      int sx = x0 < x1 ? 1 : -1;
      int sy = y0 < y1 ? 1 : -1;
      int err = (dx > dy ? dx : -dy) / 2;
      int e2;

      while (true) {
        if (x0 >= 0 && y0 >= 0 && x0 < borders.first && y0 < borders.second) {
          data[y0][x0] = 0; // Set pixel to black
        }

        if (x0 == x1 && y0 == y1) break;
        e2 = err;
        if (e2 > -dx) { err -= dy; x0 += sx; }
        if (e2 <  dy) { err += dx; y0 += sy; }
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
  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<> distribution(0, vertex_num_ * kEdgeLen_);

  for (auto vertex : graph_) {
    vertex->x = distribution(gen);
    vertex->y = distribution(gen);
  }
}

std::vector<size_t>
GraphVisualizator::KCenters(const size_t& k,
                            const std::vector<std::vector<size_t>>& distances) {
  std::vector<size_t> centers;
  centers.reserve(k);

  // Randomly choose the first center
  std::mt19937 gen(std::random_device{}());
  std::uniform_int_distribution<> distribution(0, vertex_num_ - 1);

  size_t rand_ind = distribution(gen);
  centers.emplace_back(rand_ind);

  // Memorize the current distance of every vertex from the centers
  std::vector<size_t> min_dist_to_centers(vertex_num_);
  std::transform(distances.begin(), distances.end(), min_dist_to_centers.begin(),
                   [rand_ind](const std::vector<size_t>& dist) { return dist[rand_ind]; });

  // Iteratively select the remaining centers
  for (size_t i = 1; i < k; ++i) {
    // Find the farthest vertex from the already selected centers
    auto farthest_iter = std::max_element(min_dist_to_centers.begin(), min_dist_to_centers.end());
    size_t farthest_vertex = std::distance(min_dist_to_centers.begin(), farthest_iter);
    centers.emplace_back(farthest_vertex);

    // Update min_dist_to_centers with distances to the new center
    std::transform(distances[farthest_vertex].begin(), distances[farthest_vertex].end(),
                   min_dist_to_centers.begin(), min_dist_to_centers.begin(),
                   [](size_t dist, size_t min_dist) { return std::min(dist, min_dist); });
  }

  return centers;
}

void GraphVisualizator::ComputeLocalLayout(
    const std::vector<std::vector<size_t>>& distances,
    const std::vector<size_t>& centers,
    const size_t& k) {
  // Perform local layout adjustments for a specified number of iterations
  for (size_t i = 0; i < kIterations_ * vertex_num_; ++i) {
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
  // Define a priority queue to store vertices based on their big_delta value
  std::priority_queue<std::pair<double, size_t>> pq;

  // Calculate big_delta for each center and push it into the priority queue
  for (size_t j : centers) {
    double big_delta = FindKDelta(distances, j, k);
    pq.push({big_delta, j});
  }

  // Retrieve the vertex with the maximum big_delta from the priority queue
  return pq.top().second;
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

    // Weighting constant that can be either:
    // (1 / distance[u][v]) or (1 / distance[u][v]^2 )
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
  std::pair<double, double> derivatives = FindKEnergyDerivative(distances, 
                                                                vertex_ind, k);

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

    // Weighting constant that can be either:
    // (1 / distance[u][v]) or (1 / distance[u][v]^2 )
    double k_ij = 1.0 / distances[vertex_ind][i];

    double dist = FindEuclideanDistance(vertex_ind, i);
    
    Vertex v_ind = *graph_[vertex_ind];
    Vertex v_i = *graph_[i];

    x_k_energy_derivative += k_ij * ((v_ind.x - v_i.x) - (v_ind.x - v_i.x) * 
                             kEdgeLen_ * distances[vertex_ind][i] / dist);
    
    y_k_energy_derivative += k_ij * ((v_ind.y - v_i.y) - (v_ind.y - v_i.y) *
                             kEdgeLen_ * distances[vertex_ind][i] / dist);
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