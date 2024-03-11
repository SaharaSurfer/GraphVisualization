#include "header/graph_layout.h"

#include <cmath>
#include <queue>
#include <random>

GraphLayout::GraphLayout(Graph graph) {
  graph_ = graph;
  coordinates_ = std::vector(graph.GetVertexNum(), Point(0.0, 0.0));
}

void GraphLayout::ComputeGlobalLayout() {
  // Compute the All Pairs Shortest Path (APSP) matrix for the graph
  auto APSP = graph_.GetAPSP();

  // Set up a random layout for the initial positions of vertices
  SetUpRandomLayout();

  uint32_t centers_num = kMinSize_;
  while (centers_num <= coordinates_.size()) {
    // Find the centers using the K-Centers algorithm 
    // based on the current value of k
    auto centers = graph_.FindKCenters(APSP, centers_num);
    uint32_t max_radius = 0;
    
    // Calculate the maximum radius among the centers
    for (size_t i : centers) {
      uint32_t min_dist = std::numeric_limits<uint32_t>::max();

      // Find the minimum distance from the current center to any other center
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
    // the calculated maximum radius for each center
    for (size_t center : centers) {
      ComputeLocalLayout(APSP, max_radius, center);
    }

    // Perturb the positions of vertices slightly
    // to avoid clustering around centers
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> distribution(0.1, 1);
    for (auto& [x, y] : coordinates_) {
      x += distribution(gen);
      y += distribution(gen);
    }

    // Increase number of the centers for the next iteration
    centers_num *= kRatio_;
  }
}

void GraphLayout::SetUpRandomLayout() {
  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<> distribution(0, coordinates_.size() * kEdgeLen_);

  for (auto& [x, y] : coordinates_) {
    x = distribution(gen);
    y = distribution(gen);
  }
}

std::vector<Point> GraphLayout::GetCoordinates() {
  return coordinates_;
}

void GraphLayout::ComputeLocalLayout(
    const std::vector<std::vector<uint32_t>>& APSP,
    const uint32_t& neighbourhood_radius,
    const size_t& center) {
  std::vector<size_t> k_neighbourhood = graph_.FindNeighbourhood(
                                            APSP, neighbourhood_radius, center);

  // Perform local layout adjustments for a specified number of iterations
  for (size_t i = 0; i < kIterations_ * coordinates_.size(); ++i) {
    size_t ind = ChooseVertex(APSP, k_neighbourhood, neighbourhood_radius);

    // Find displacement
    auto [disp_x, disp_y] = FindSmallDelta(APSP, neighbourhood_radius, ind);
    
    // Update vertex
    coordinates_[ind].x += disp_x;
    coordinates_[ind].y += disp_y;
  }
}

size_t GraphLayout::ChooseVertex(const std::vector<std::vector<uint32_t>>& APSP,
                                 const std::vector<size_t>& k_neighbourhood,
                                 const uint32_t& neighbourhood_radius) {
  // Define a priority queue to store vertices based on their big_delta value
  std::priority_queue<std::pair<double, size_t>> pq;

  // Calculate big_delta for each center and push it into the priority queue
  for (size_t j : k_neighbourhood) {
    double big_delta = FindBigDelta(APSP, neighbourhood_radius, j);
    pq.push({big_delta, j});
  }

  // Retrieve the vertex with the maximum big_delta from the priority queue
  return pq.top().second;
}

double GraphLayout::FindBigDelta(const std::vector<std::vector<uint32_t>>& APSP,
                                 const uint32_t& neighbourhood_radius,
                                 const size_t& vertex_ind) {
  auto [derivative_x, derivative_y] = FindFirstEnergyDerivative(
                                          APSP, 
                                          neighbourhood_radius, 
                                          vertex_ind);

  return std::sqrt(std::pow(derivative_x, 2) +
                   std::pow(derivative_y, 2));
}

std::pair<double, double> GraphLayout::FindFirstEnergyDerivative(
    const std::vector<std::vector<uint32_t>>& APSP,
    const uint32_t& neighbourhood_radius,
    const size_t& vertex_ind) {
  // Find the k-neighbourhood of the vertex
  std::vector<size_t> k_neighbourhood = graph_.FindNeighbourhood(
                                            APSP, neighbourhood_radius, vertex_ind);
  double derivative_x = 0.0;
  double derivative_y = 0.0;

  // Iterate over the vertices in the k-neighbourhood to compute derivatives
  for (size_t i : k_neighbourhood) {
    if (i == vertex_ind) { continue; }

    // Weighting constant that can be either:
    // (1 / distance[u][v]) or (1 / distance[u][v]^2 )
    const double k_ij = 1.0 / APSP[vertex_ind][i];

    const double dist = FindEuclideanDistance(vertex_ind, i);
    
    const auto& [x1, y1] = coordinates_[vertex_ind];
    const auto& [x2, y2] = coordinates_[i];

    derivative_x += k_ij * ((x1 - x2) - (x1 - x2) * 
                             kEdgeLen_ * APSP[vertex_ind][i] / dist);
    
    derivative_y += k_ij * ((y1 - y2) - (y1 - y2) *
                             kEdgeLen_ * APSP[vertex_ind][i] / dist);
  }

  return {derivative_x, derivative_y};
}

double GraphLayout::FindEuclideanDistance(const size_t& a, 
                                          const size_t& b) {
  const auto& [x1, y1] = coordinates_[a];
  const auto& [x2, y2] = coordinates_[b];
  
  return std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
}

std::pair<double, double> GraphLayout::FindSmallDelta(
    const std::vector<std::vector<uint32_t>>& APSP,
    const uint32_t& neighbourhood_radius,
    const size_t& vertex_ind) {
  // Compute second derivatives of the energy function
  auto second_derivatives = FindSecondEnergyDerivatives(APSP, 
                                                        neighbourhood_radius,
                                                        vertex_ind);

  // Compute first derivatives of the energy function
  auto first_derivatives = FindFirstEnergyDerivative(APSP, 
                                                     neighbourhood_radius,
                                                     vertex_ind);
  
  // Extract coefficients for the linear system
  double a_1 = second_derivatives[0];
  double b_1 = second_derivatives[1];
  double c_1 = -first_derivatives.first;

  double a_2 = second_derivatives[1];
  double b_2 = second_derivatives[2];
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

std::vector<double> GraphLayout::FindSecondEnergyDerivatives(
    const std::vector<std::vector<uint32_t>>& APSP,
    const uint32_t& neighbourhood_radius,
    const size_t& vertex_ind) {
  std::vector<size_t> k_neighbourhood = graph_.FindNeighbourhood(
                                            APSP, neighbourhood_radius, vertex_ind);
  double derivative_x_x = 0.0;
  double derivative_x_y = 0.0;
  double derivative_y_y = 0.0;

  for (size_t i : k_neighbourhood) {
    if (i == vertex_ind) { continue; }

    // Weighting constant that can be either:
    // (1 / distance[u][v]) or (1 / distance[u][v]^2 )
    const double k_ij = 1.0 / APSP[vertex_ind][i];

    const double dist = FindEuclideanDistance(vertex_ind, i);

    const auto& [x1, y1] = coordinates_[vertex_ind];
    const auto& [x2, y2] = coordinates_[i];

    derivative_x_x += k_ij * (1.0 - kEdgeLen_ * 
                              APSP[vertex_ind][i] * 
                              std::pow(y1 - y2, 2) /
                              std::pow(dist, 3));
    
    derivative_x_y += k_ij * kEdgeLen_ * 
                      APSP[vertex_ind][i] * 
                      (y1 - y2) * (x1 - x2) /
                      std::pow(dist, 3);

    derivative_y_y += k_ij * (1.0 - kEdgeLen_ * 
                              APSP[vertex_ind][i] * 
                              std::pow(x1 - x2, 2) /
                              std::pow(dist, 3));
  }

  return {derivative_x_x, derivative_x_y, derivative_y_y};
}
