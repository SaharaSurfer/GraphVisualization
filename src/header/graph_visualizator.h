#ifndef GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H
#define GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H

#include <string>
#include <vector>

class GraphVisualizator {
public:
  GraphVisualizator() = default;

  std::vector<std::vector<int>> ReadGraph(const std::string& filename);

private:
  size_t vertex_num_ = 0;
  size_t edge_num_ = 0;
};

#endif  // GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H