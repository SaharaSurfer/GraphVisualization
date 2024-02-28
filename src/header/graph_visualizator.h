#ifndef GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_
#define GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_

#include <cstdint>
#include <deque>
#include <memory>
#include <string>
#include <utility>
#include <vector>

struct Vertex {
  Vertex(size_t _number) { number = _number; }

  size_t number;
  double x;
  double y;
  std::vector<std::shared_ptr<Vertex>> neighbours;
};

class GraphVisualizator {
public:
  GraphVisualizator() = default;

  void ReadGraph(const std::string& filename);

private:
  std::vector<size_t> BFS(std::shared_ptr<Vertex> root);
  std::vector<size_t> CreateVertexFiltration();
  void PlaceCoreVertices();

  size_t vertex_num_ = 0;
  size_t edge_num_ = 0;
  std::vector<std::shared_ptr<Vertex>> graph_;

  const size_t kEdgeLen_ = 30;
  const size_t kNeighbourNumber_ = 3;
  const size_t kRepeatBeforeReset_ = 1000;
  const size_t kResetThreshhold_ = 0;
};

#endif  // GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_