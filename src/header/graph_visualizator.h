#ifndef GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_
#define GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_

#include <cstdint>
#include <memory>
#include <string>
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

  size_t vertex_num_ = 0;
  size_t edge_num_ = 0;
  std::vector<std::shared_ptr<Vertex>> graph_;

  // Desired length of an edge in the visualization
  const size_t kEdgeLen_ = 30;
};

#endif  // GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_
