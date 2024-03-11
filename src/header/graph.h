#ifndef GRAPHVISUALIZATION_HEADER_GRAPH_H_
#define GRAPHVISUALIZATION_HEADER_GRAPH_H_

#include <cstdint>
#include <string>
#include <vector>
#include <memory>

// Vertex struct represents a vertex in a graph.
// It stores the vertex number 
// and a list of neighboring vertices.
struct Vertex {
  Vertex(size_t _number) : number{_number} {}

  size_t number;
  std::vector<std::weak_ptr<Vertex>> neighbours;
};

class Graph {
 public:
  Graph() = default;

  // Reads the graph data from a file and initializes 
  // the graph vertices and edges.
  // Parameters:
  // - filename: the name of the file containing the graph data
  void ReadGraph(const std::string& filename);

  // Performs Breadth-First Search (BFS) traversal on the graph
  // starting from the vertex with the given root index.
  // Parameters:
  // - root_index: the index of the root vertex for BFS traversal
  // Returns:
  // - a vector containing the distances from the root to all vertices
  std::vector<uint32_t> BFS(const size_t& root_index);

  // Computes the All Pairs Shortest Path (APSP) matrix
  // using BFS algorithm.
  // Returns:
  // - A 2D vector representing the APSP matrix, where APSP[i][j] is 
  //   the shortest distance from vertex i to vertex j.
  std::vector<std::vector<uint32_t>> GetAPSP();

  // Finds k centers in the graph using the K-Centers algorithm
  // based on the provided APSP matrix and the number of centers.
  // Parameters:
  // - APSP: a 2D vector representing the All Pairs Shortest Path matrix
  // - centers_num: the number of centers to find
  // Returns:
  // - A vector containing the indices of the selected centers.
  std::vector<size_t> FindKCenters(
      const std::vector<std::vector<uint32_t>>& APSP,
      const uint32_t& centers_num);

  // Finds the neighbourhood of a vertex within a given radius
  // based on the provided APSP matrix, neighbourhood radius, 
  // and vertex index.
  // Parameters:
  // - APSP: a 2D vector representing the All Pairs Shortest Path matrix
  // - neigbourhood_radius: the maximum distance for vertices to be 
  //                        considered in the neighbourhood
  // - vertex_ind: the index of the vertex for which 
  //               the k-neighbourhood is to be found
  // Returns:
  // - A vector containing the indices of vertices within
  //   the k-neighbourhood of the given vertex
  std::vector<size_t> FindNeighbourhood(
      const std::vector<std::vector<uint32_t>>& APSP,
      const uint32_t& neigbourhood_radius, 
      const size_t& vertex_ind);

  // Returns:
  // - The number of vertices in the graph
  uint32_t GetVertexNum();

  // Returns:
  // - A vector containing shared pointers to the vertices of the graph
  std::vector<std::shared_ptr<Vertex>> GetVertices();

 private:
  std::vector<std::shared_ptr<Vertex>> vertices_;
  uint32_t vertex_num_ = 0;
  uint32_t edge_num_ = 0;
};

#endif  // GRAPHVISUALIZATION_HEADER_GRAPH_H_
