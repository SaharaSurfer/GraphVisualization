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
  size_t x;
  size_t y;
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

  const size_t kLastSetSize_ = 3;
};

#endif  // GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_