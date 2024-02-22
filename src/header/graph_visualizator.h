#ifndef GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_
#define GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_

#include <string>
#include <vector>
#include <utility>

class GraphVisualizator {
public:
  GraphVisualizator() = default;

  void ReadGraph(const std::string& filename);

private:
  size_t vertex_num_ = 0;
  size_t edge_num_ = 0;
  std::vector<std::vector<size_t>> graph_;
};

#endif  // GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_