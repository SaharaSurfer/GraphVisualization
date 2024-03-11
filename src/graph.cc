#include "header/graph.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <queue>
#include <random>
#include <string>
#include <vector>
#include <memory>

// Reads the graph data from a file and initializes 
// the graph vertices and edges.
void Graph::ReadGraph(const std::string& filename) {
  std::ifstream input(filename);

  if (!input.is_open()) {
    throw std::runtime_error("Failed to open file");
  }
  
  input >> vertex_num_ >> edge_num_;
  
  // Create and initialize vertices for the graph
  vertices_.resize(vertex_num_);
  for (size_t i = 0; i < vertex_num_; ++i) {
    vertices_[i] = std::make_shared<Vertex>(i);
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

    vertices_[ind_1 - 1]->neighbours.emplace_back(vertices_[ind_2 - 1]);
    vertices_[ind_2 - 1]->neighbours.emplace_back(vertices_[ind_1 - 1]);
  }

  input.close();
}

std::vector<uint32_t> Graph::BFS(const size_t& root_ind) {
  // Since the distance in the graph between 2 vertices cannot be greater
  // than vertex_num_ - 1, the value of vertex_num_ can be used
  // to initialise unreachability.
  std::vector<uint32_t> dist(vertex_num_, vertex_num_);
  dist[root_ind] = 0;

  std::queue<std::shared_ptr<Vertex>> queue;
  queue.emplace(vertices_[root_ind]);

  while (!queue.empty()) {
    auto vertex = queue.front();
    queue.pop();

    // Traverse the neighbours of the current vertex
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

std::vector<std::vector<uint32_t>> Graph::GetAPSP() {
  std::vector<std::vector<uint32_t>> APSP(vertex_num_);

  // Compute shortest paths from each vertex to all other vertices using BFS
  for (size_t i = 0; i < vertex_num_; ++i) {
    APSP[i] = BFS(i);
  }

  return APSP;
}

std::vector<size_t> Graph::FindKCenters(
    const std::vector<std::vector<uint32_t>>& APSP,
    const uint32_t& centers_num) {
  std::vector<size_t> centers;
  centers.reserve(centers_num);

  // Randomly choose the first center
  std::mt19937 gen(std::random_device{}());
  std::uniform_int_distribution<> distribution(0, vertex_num_ - 1);

  size_t rand_ind = distribution(gen);
  centers.emplace_back(rand_ind);

  // Memorize the current distance of every vertex from the centers
  std::vector<uint32_t> min_dist_to_centers = APSP[rand_ind];

  // Iteratively select the remaining centers
  for (size_t i = 1; i < centers_num; ++i) {
    // Find the farthest vertex from the already selected centers
    auto farthest_iter = std::max_element(min_dist_to_centers.begin(), min_dist_to_centers.end());
    size_t farthest_vertex = std::distance(min_dist_to_centers.begin(), farthest_iter);
    centers.emplace_back(farthest_vertex);

    // Update min_dist_to_centers with distances to the new center
    std::transform(APSP[farthest_vertex].begin(), APSP[farthest_vertex].end(),
                   min_dist_to_centers.begin(), min_dist_to_centers.begin(),
                   [](size_t dist, size_t min_dist) { return std::min(dist, min_dist); });
  }

  return centers;
}

std::vector<size_t> Graph::FindNeighbourhood(
    const std::vector<std::vector<uint32_t>>& APSP,
    const uint32_t& neighbourhood_radius, 
    const size_t& vertex_ind) {
  std::vector<size_t> k_neighbourhood;

  // Traverse all vertices in the graph
  for (size_t i = 0; i < vertex_num_; ++i) {
    // Check if the distance from the current vertex to the 
    // target vertex is within the threshold
    if (APSP[vertex_ind][i] < neighbourhood_radius) {
      k_neighbourhood.emplace_back(i);
    }
  }

  return k_neighbourhood;
}

uint32_t Graph::GetVertexNum() {
  return vertex_num_;
}

std::vector<std::shared_ptr<Vertex>> Graph::GetVertices() {
  return vertices_;
}
