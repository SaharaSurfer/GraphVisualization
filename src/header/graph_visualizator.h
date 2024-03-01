#ifndef GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_
#define GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_

#include <cstdint>
#include <deque>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Vertex struct represents a vertex in a graph.
// It stores the vertex number, coordinates (x, y), 
// and a list of neighboring vertices.
struct Vertex {
  Vertex(size_t _number) { number = _number; }

  size_t number;
  double x = 0.0;
  double y = 0.0;
  std::vector<std::shared_ptr<Vertex>> neighbours;
  
  // The highest number of the filtration set that contains the given vertex.
  int depth = 0;
};

// GraphVisualizator class is responsible for visualizing graphs.
// It provides functions to read a graph from a file, perform graph algorithms,
// and place vertices in a visualization space.
class GraphVisualizator {
public:
  GraphVisualizator() = default;

  // ReadGraph function reads a graph from a file 
  // and initializes the graph data structure.
  // Parameters:
  // - filename: the name of the file containing the graph data
  void ReadGraph(const std::string& filename);

private:
  // BFS function performs breadth-first search traversal on the 
  // graph starting from a given root vertex.
  // Parameters:
  // - root: a shared pointer to the root vertex for BFS traversal
  // Returns:
  // - a vector containing the distances from the root to all vertices
  std::vector<size_t> BFS(std::shared_ptr<Vertex> root);

  // CreateVertexFiltration function generates a vertex filtration 
  // (V = V_0 > V_1 > ... > V_k > 0) for the graph using a random process.
  // Each V_i is a maximal subset of V_{i - 1} for which the graph distance
  // between any pair of its elements is at least 2^{i - 1} + 1.
  // The function swaps the elements of the graph so that the first |V_k| 
  // elements belong to V_k, the second |V_{k-1}| elements 
  // belong to V_{k-1}, and so on.
  // Returns:
  // - a vector containing the sizes of each filtration stage
  std::vector<size_t> CreateVertexFiltration();

  // Finds sets of nearest neighbors of a graph vertex for each filter set.
  // Parameters:
  //   - borders: A vector containing the boundaries of each filter set.
  //   - root: The root vertex from which to start the search.
  // Returns:
  // A vector of vectors containing the nearest neighbor sets for each filter set.
  std::vector<std::vector<std::shared_ptr<Vertex>>>
  FindNearestNeighbourSets(const std::vector<size_t>& borders,
                           std::shared_ptr<Vertex> root);

  // PlaceCoreVertices function calculates the positions of 
  // core vertices in the graph.
  void PlaceCoreVertices();

  // Finds the three nearest vertices to a given root vertex in the graph.
  // Parameters:
  //   - borders: A vector containing the boundaries of each filter set.
  //   - root: The root vertex from which to start the search.
  // Returns:
  // A vector of pairs containing the distance to root and the index in graph_.
  std::vector<std::pair<size_t, size_t>>
  FindThreeNearestVertices(const std::vector<size_t>& borders,
                           std::shared_ptr<Vertex> root);

  size_t vertex_num_ = 0;
  size_t edge_num_ = 0;
  std::vector<std::shared_ptr<Vertex>> graph_;

  // Desired length of an edge in the visualization
  const size_t kEdgeLen_ = 30;

  // Number of elements to remain in V_k
  const size_t kCoreVertexNumber_ = 3;

  // Number of repetitions before resetting last set of filtration
  const size_t kRepeatBeforeReset_ = 1000;
  
  // Threshold for resetting filtration
  const size_t kResetThreshhold_ = 0;
};

#endif  // GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_
