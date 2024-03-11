#ifndef GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_
#define GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "graph.h"
#include "graph_layout.h"

// GraphVisualizator class is responsible for visualizing graphs.
// It provides functions to draw graphs and save them as bitmap images.
class GraphVisualizator {
public:
  GraphVisualizator() = default;

  // Draws the graph and saves it as a bitmap image.
  // Parameters:
  // - graph: The graph to be visualized
  // - filename: The name of the file to save the bitmap image
  void DrawBmp(Graph graph, const std::string& filename);

private:
  // Desired radius of a circle presenting vertex in the visualization
  const int kCircleRadius_ = 4;

  const int kPadding_ = 30;

  // Prepares the data for visualization, including vertex coordinates.
  // Parameters:
  // - graph: Pointer to the graph to be visualized
  // - coordinates: Pointer to the vector containing the coordinates of the vertices
  // Returns:
  // - A 2D vector representing the bitmap image data
  std::vector<std::vector<uint8_t>> PrepareData(Graph* graph,
                                                std::vector<Point>* coordinates);

  // Rounds all vertex coordinates to the nearest integer number
  // Parameters:
  // - coordinates: Pointer to the vector containing the coordinates of the vertices
  void RoundVertexCoordinates(std::vector<Point>* coordinates);

  // Moves vertices closer to (0; 0) and ensures all coordinates are positive
  // Parameters:
  // - coordinates: Pointer to the vector containing the coordinates of the vertices
  // Returns:
  // - A pair representing the edge of the image.
  std::pair<int, int> CorrectCoordinates(std::vector<Point>* coordinates);

  //void DrawCircle();

  //void DrawNumber();

  //void DrawEdge();
};

#endif  // GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_
