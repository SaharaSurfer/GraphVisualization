#ifndef GRAPHVISUALIZATION_HEADER_GRAPH_LAYOUT_H_
#define GRAPHVISUALIZATION_HEADER_GRAPH_LAYOUT_H_

#include <vector>

#include "graph.h"

// Point struct represents a point in a 2D space.
// It stores the x and y coordinates of the point.
struct Point {
  Point(double _x, double _y) : x(_x), y(_y) {}

  double x = 0.0;
  double y = 0.0;
};

// GraphLayout class is responsible for computing the layout of a graph.
// It provides functions to compute both global and local layouts of the graph.
class GraphLayout {
 public:
  GraphLayout(Graph graph);

  // Computes the global layout of the graph.
  void ComputeGlobalLayout();
  
  // Sets up a random layout for the vertices in the graph.
  // Assigns random x and y coordinates to each vertex.
  void SetUpRandomLayout();
  
  // Returns:
  // - A vector containing the coordinates (x, y) of the vertices.
  std::vector<Point> GetCoordinates();
  
 private:
  // Desired length of an edge in the visualization
  const int kEdgeLen_ = 30;

  // Determines radius of local neighbourhoods
  const int kRadius_ = 7;

  // Determines number of iteration in local beautification;
  const int kIterations_ = 4;

  // Ratio between number of vertices in two consecutive levels
  const int kRatio_ = 3;

  // Size of the coarsest graph
  const int kMinSize_ = 3;

  // Computes the local layout adjustments for vertices based on the K-Centers algorithm.
  // Iterates through a specified number of iterations, adjusting the positions
  // of vertices locally according to the K-Centers algorithm.
  // Parameters:
  // - APSP: A 2D vector representing the All Pairs Shortest Path matrix
  // - neighbourhood_radius: The maximum distance for vertices to be considered
  //                         in the neighbourhood
  // - center: The index of the center vertex for local layout computation
  void ComputeLocalLayout(const std::vector<std::vector<uint32_t>>& APSP,
                          const uint32_t& neighbourhood_radius,
                          const size_t& center);

  // Chooses the vertex with the maximum Big Delta (just a criterion) 
  // for local layout adjustment.
  // Parameters:
  // - APSP: A 2D vector representing the All Pairs Shortest Path matrix
  // - k_neighbourhood: A vector containing the indices of vertices within
  //                    the k-neighbourhood of the center vertex
  // - neighbourhood_radius: The maximum distance for vertices to be considered
  //                         in the neighbourhood
  // Returns:
  // - The index of the selected vertex with the largest displacement.
  size_t ChooseVertex(const std::vector<std::vector<uint32_t>>& APSP,
                      const std::vector<size_t>& k_neighbourhood,
                      const uint32_t& neighbourhood_radius);
  
  // Finds the maximum Big Delta needed for local layout adjustment.
  // Parameters:
  // - APSP: A 2D vector representing the All Pairs Shortest Path matrix
  // - neighbourhood_radius: The maximum distance for vertices to be considered
  //                         in the neighbourhood
  // - vertex_ind: The index of the vertex for which the big delta is calculated
  // Returns:
  // - The maximum displacement needed to adjust the position of the vertex.
  double FindBigDelta(const std::vector<std::vector<uint32_t>>& APSP,
                      const uint32_t& neighbourhood_radius,
                      const size_t& vertex_ind);

  // Finds the small displacement (delta) needed for local layout adjustment.
  // Parameters:
  // - APSP: A 2D vector representing the All Pairs Shortest Path matrix
  // - neighbourhood_radius: The maximum distance for vertices to be considered
  //                         in the neighbourhood
  // - vertex_ind: The index of the vertex for which the small delta is calculated
  // Returns:
  // - A pair containing the calculated displacements along the x and y axes.
  std::pair<double, double> FindSmallDelta(
      const std::vector<std::vector<uint32_t>>& APSP,
      const uint32_t& neighbourhood_radius,
      const size_t& vertex_ind);

  // Finds the derivative of the energy function for local layout adjustment.
  // Parameters:
  // - APSP: A 2D vector representing the All Pairs Shortest Path matrix
  // - neighbourhood_radius: The maximum distance for vertices to be considered
  //                         in the neighbourhood
  // - vertex_ind: The index of the vertex for which the energy derivative is calculated
  // Returns:
  // - A pair containing the derivatives of the energy function with
  //   respect to x and y coordinates.
  std::pair<double, double> FindFirstEnergyDerivative(
      const std::vector<std::vector<uint32_t>>& APSP,
      const uint32_t& neighbourhood_radius,
      const size_t& vertex_ind);
  
  // Finds the second derivatives of the energy function for local layout adjustment.
  // Parameters:
  // - distances: A 2D vector representing the distances between vertices in the graph.
  // - neighbourhood_radius: The maximum distance for vertices to be considered
  //                         in the neighbourhood
  // - vertex_ind: The index of the vertex for which the derivatives are calculated
  // Returns:
  // - A vector containing the second derivatives of the energy function with
  //   respect to x and y (x_x, x_y, y_y).
  std::vector<double> FindSecondEnergyDerivatives(
      const std::vector<std::vector<uint32_t>>& distances,
      const uint32_t& neighbourhood_radius,
      const size_t& vertex_ind);

  // Calculates the euclidean distance between two vertices with given indices.
  // Parameters:
  // - a: The index of the first vertex
  // - b: The index of the second vertex
  // Returns:
  // - The euclidean distance between the two vertices.
  double FindEuclideanDistance(const size_t& a, const size_t& b);

  Graph graph_;
  std::vector<Point> coordinates_;
};

#endif  // GRAPHVISUALIZATION_HEADER_GRAPH_LAYOUT_H_
