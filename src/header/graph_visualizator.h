#ifndef GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_
#define GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_

#include <string>
#include <vector>
#include <deque>
#include <memory>
#include <cstdint>
#include <utility>

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
  std::vector<size_t> BFS(std::shared_ptr<Vertex> root);
  std::vector<size_t> CreateVertexFiltration();

  std::vector<std::shared_ptr<Vertex>> graph_;

private:
  size_t vertex_num_ = 0;
  size_t edge_num_ = 0;
};

#endif  // GRAPHVISUALIZATION_HEADER_GRAPH_VISUALIZATOR_H_