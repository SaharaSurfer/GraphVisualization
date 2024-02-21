#ifndef GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_
#define GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_

#include <string>
#include <vector>
#include <utility>

class GraphVisualizator {
public:
  GraphVisualizator() = default;

  void ReadGraph(const std::string& filename);

  std::vector<std::vector<size_t>> SearchForConnectivityComponents();
  void FindConnectivityComponent(
      const size_t& vertex_ind,
      std::vector<bool> *visited,
      std::vector<std::vector<size_t>> *connectivity_components);

  std::vector<std::pair<size_t, size_t>> SearchForBridges();
  void FindBridge(const size_t& vertex_ind,
                  const size_t& parent_ind,
                  std::vector<bool> *visited,
                  std::vector<size_t> *time_in,
                  std::vector<size_t> *min_time_in_reachable,
                  std::vector<std::pair<size_t, size_t>> *bridges);
  
  std::vector<size_t> GetCycle();
  bool DetectCycle(const size_t& vertex_ind,
                   size_t *cycle_start, size_t *cycle_end,
                   std::vector<char> *vertex_color,
                   std::vector<size_t> *parents);

private:
  size_t vertex_num_ = 0;
  size_t edge_num_ = 0;
  std::vector<std::vector<size_t>> graph_;
};

#endif  // GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_