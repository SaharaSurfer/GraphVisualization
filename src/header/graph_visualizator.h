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
  std::vector<std::weak_ptr<Vertex>> neighbours;
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

  // Computes the global layout adjustments for vertices 
  // based on the K-Centers algorithm.
  void ComputeGlobalLayout();

private:
  // BFS function performs breadth-first search traversal on the 
  // graph starting from a given root vertex.
  // Parameters:
  // - root: a shared pointer to the root vertex for BFS traversal
  // Returns:
  // - a vector containing the distances from the root to all vertices
  std::vector<size_t> BFS(const std::shared_ptr<Vertex>& root);

  // Finds the shortest paths between all pairs of 
  // vertices in the graph using BFS.
  // Returns:
  // - A 2D vector representing the APSP matrix, where APSP[i][j] is 
  //   the shortest distance from vertex i to vertex j.
  std::vector<std::vector<size_t>> FindAllPairsShortestPath();
  
  // Sets up a random layout for the vertices in the graph.
  // Assigns random x and y coordinates to each vertex.
  void SetUpRandomLayout();

  // Performs K-Centers algorithm
  // Parameters:
  // - k: The number of centers to find.
  // - distances: APSP matrix
  // Returns:
  // - A vector containing the indices of the selected centers.
  std::vector<size_t> KCenters(const size_t& k,
                               const std::vector<std::vector<size_t>>& distances);
  
  // Finds the farthest vertex from the selected centers
  // Parameters:
  // - distances: APSP matrix
  // - in_centers: A vector indicating whether each vertex is already
  //   selected as a center.
  // Returns:
  // - The index of the farthest vertex from the selected centers.
  size_t FindFarthestVertex(const std::vector<std::vector<size_t>>& distances,
                            const std::vector<bool>& in_centers);

  // Computes the local layout adjustments for vertices based on the K-Centers algorithm.
  // Iterates through a specified number of iterations, adjusting the positions
  // of vertices locally according to the K-Centers algorithm.
  // Parameters:
  // - distances: A 2D vector representing the distances between vertices in the graph.
  // - k: The parameter used in the K-Centers algorithm.
  void ComputeLocalLayout(const std::vector<std::vector<size_t>>& distances,
                          const std::vector<size_t>& centers,
                          const size_t& k);

  // Iterates through all vertices and selects the one 
  // with the maximum FindKDelta.
  // Parameters:
  // - distances: A 2D vector representing the distances between
  // vertices in the graph.
  // - k: The parameter used in the K-Centers algorithm.
  // Returns:
  // - The index of the selected vertex with the largest displacement.
  size_t ChooseVertex(const std::vector<std::vector<size_t>>& distances,
                      const std::vector<size_t>& centers,
                      const size_t& k);

  // Calculates the small displacement (delta) for the given
  // vertex based on the k-neighbourhood.
  // Parameters:
  // - distances: A 2D vector representing the distances between
  //   vertices in the graph.
  // - vertex_ind: The index of the vertex for which the 
  //   displacement is calculated.
  // - k: The parameter used to define the k-neighbourhood.
  // Returns:
  // - A pair containing the calculated displacements
  // along the x and y axes.
  std::pair<double, double> FindKSmallDelta(
      const std::vector<std::vector<size_t>>& distances,
      const size_t& vertex_ind,
      const size_t& k);

  // Calculates the partial derivatives of the energy function
  // with respect to x and y.
  // Parameters:
  // - distances: A 2D vector representing the distances between
  //   vertices in the graph.
  // - vertex_ind: The index of the vertex for which the partial
  //   derivatives are calculated.
  // - k: The parameter used to define the k-neighbourhood.
  // Returns:
  // - A vector containing the partial derivatives of the energy 
  //   function with respect to x and y.
  std::vector<double> FindKPartialDerivatives(
      const std::vector<std::vector<size_t>>& distances,
      const size_t& vertex_ind,
      const size_t& k);

  // Calculates the change (∆^k_u) needed to choose the vertex to move.
  // Computation is based on the derivatives of the energy function.
  // Parameters:
  // - distances: A 2D vector representing the distances between
  //   vertices in the graph.
  // - vertex_ind: The index of the vertex for which the change is calculated.
  // - k: The parameter used to define the k-neighbourhood.
  // Returns:
  // - The change (∆^k_u) needed to choose the vertex to move.
  double FindKDelta(const std::vector<std::vector<size_t>>& distances,
                    const size_t& vertex_ind,
                    const size_t& k);
  
  // Finds the derivative of the energy function for both x and y coordinates.
  // The energy function represents the force exerted on a 
  // vertex due to its neighbours.
  // Parameters:
  // - distances: A 2D vector representing the distances between
  //              vertices in the graph.
  // - vertex_ind: The index of the vertex for which
  //               the energy derivative is calculated.
  // - k: The parameter used to define the k-neighbourhood.
  // Returns:
  // - A pair containing the derivatives of the energy function
  //   with respect to x and y coordinates.
  std::pair<double, double> FindKEnergyDerivative(
      const std::vector<std::vector<size_t>>& distances,
      const size_t& vertex_ind,
      const size_t& k);

  // Finds neighbourhood of vertex "v" defined
  // to be N^k(v) = {u ∈ V | 0 ≤ d_{uv} < k}
  // Parameters:
  // - distances: A 2D vector representing the distances between 
  //              vertices in the graph.
  // - vertex_ind: The index of the vertex for which 
  //               the k-neighbourhood is to be found.
  // - k: The maximum distance for vertices to be 
  //      considered in the neighbourhood.
  // Returns:
  // - A vector containing the indices of vertices within
  //   the k-neighbourhood of the given vertex.
  std::vector<size_t> FindKNeighbourhood(
      const std::vector<std::vector<size_t>>& distances,
      const size_t& vertex_ind, 
      const size_t& k);

  // Calculates euclidean distance between 2 vertices with given indices
  double FindEuclideanDistance(const size_t& a, const size_t& b);

  size_t vertex_num_ = 0;
  size_t edge_num_ = 0;
  std::vector<std::shared_ptr<Vertex>> graph_;

  // Desired length of an edge in the visualization
  const int kEdgeLen_ = 30;

  // Determines radius of local neighbourhoods
  const int kRadius_ = 7;

  // Determines number of iteration in local beautification;
  const int kIterations_ = 4;

  // Ratio between number of vertices in two consecutive levels
  const int kRatio_ = 3;

  // Size of the coarsest graph
  const int kMinSize_ = 10;
};

#endif  // GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_
