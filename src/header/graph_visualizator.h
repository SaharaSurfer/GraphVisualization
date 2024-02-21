#ifndef GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_
#define GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_

#include <string>
#include <vector>

class GraphVisualizator {
public:
  GraphVisualizator() = default;

  void ReadGraph(const std::string& filename);
  std::vector<std::vector<size_t>> SearchForConnectivityComponents();
  void FindConnectivityComponent(
      const size_t &vertex_ind,
      std::vector<bool> *visited,
      std::vector<std::vector<size_t>> *connectivity_components);

private:
  size_t vertex_num_ = 0;
  size_t edge_num_ = 0;
  std::vector<std::vector<size_t>> graph_;
};

#endif  // GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_