#include "header/graph_visualizator.h"

#include <algorithm>  // For std::min, std::max, std::transform, std::max_element
#include <cmath>  // For std::sqrt, std::pow, std::round
#include <cstdint>  // For size_t, uint8_t
#include <fstream>  // For std::ifstream
#include <iterator> // For std::distance
#include <limits>  // For std::numeric_limits
#include <queue>  // For std::queue, std::priority_queue
#include <random>  // For random number generation
#include <stdexcept> // For std::runtime_error, std::logic_error
#include <string>  // For std::string
#include <vector>  // For std::vector
#include <memory> // For std::shared_ptr

#include "header/bmp.h"
#include "header/graph.h"
#include "header/graph_layout.h"

void GraphVisualizator::DrawBmp(Graph graph, 
                                const std::string& filename) {
  GraphLayout layout(graph);
  layout.ComputeGlobalLayout();

  auto coordinates = layout.GetCoordinates();
  auto data = PrepareData(&graph, &coordinates);

  Bmp image;
  image.Interpret(data);
  image.Write(filename);
}

std::vector<std::vector<uint8_t>> GraphVisualizator::PrepareData(
    Graph* graph,
    std::vector<Point>* coordinates) {
  RoundVertexCoordinates(coordinates);
  auto borders = CorrectCoordinates(coordinates);

  // Add padding
  borders.first += (4 - borders.first % 4) % 4;

  // Initialize data matrix with 255 (white colour)
  std::vector<std::vector<uint8_t>>
  data(borders.second, std::vector<uint8_t>(borders.first, 255));

  // Add circle around vertex position
  for (size_t i = 0; i < (*coordinates).size(); ++i) {
    auto& [x0, y0] = (*coordinates)[i];
    for (int y = -kCircleRadius_; y <= kCircleRadius_; ++y) {
      for (int x = -kCircleRadius_; x <= kCircleRadius_; ++x) {
        if (x * x + y * y <= kCircleRadius_ * kCircleRadius_) {
          int posX = static_cast<int>(x0) + x;
          int posY = static_cast<int>(y0) + y;

          if (posX >= 0 && posY >= 0 && posX < borders.first && posY < borders.second) {
            data[posY][posX] = 0; // Set pixel to black
          }
        }
      }
    }

    // Add numbers of vertices next to them
    std::string vertex_number_string = std::to_string(i + 1);
    size_t start_x = x0 + 2 * kCircleRadius_;
    size_t start_y = y0 - kCircleRadius_;

    for (auto& digit : vertex_number_string) {
      std::string filename(1, digit);
      filename = "../src/digits/" + filename + ".bmp";

      Bmp digit_bmp;
      digit_bmp.Read(filename);

      size_t width = digit_bmp.info_header_.width;
      size_t height = digit_bmp.info_header_.height;

      for (size_t y = 0; y < height; ++y) {
        for (size_t x = 0; x < width; ++x) {
          data[start_y + y][start_x + x] = digit_bmp.data_[y * width + x];
        }
      }

      start_x += width;
    }
  }

  // Add edges between vertices
  auto vertices = (*graph).GetVertices();
  for (const auto& vertex : vertices) {
    for (const auto& neighbour_weak : vertex->neighbours) {
      auto neighbour = neighbour_weak.lock();

      if (!neighbour) {
        throw std::logic_error("Neighbour is a null pointer");
      }

      // Bresenham's line algorithm to draw a line between two points
      int x0 = static_cast<int>((*coordinates)[vertex->number].x);
      int y0 = static_cast<int>((*coordinates)[vertex->number].y);
      int x1 = static_cast<int>((*coordinates)[neighbour->number].x);
      int y1 = static_cast<int>((*coordinates)[neighbour->number].y);

      int dx = std::abs(x1 - x0);
      int dy = std::abs(y1 - y0);
      int sx = x0 < x1 ? 1 : -1;
      int sy = y0 < y1 ? 1 : -1;
      int err = (dx > dy ? dx : -dy) / 2;
      int e2;

      while (true) {
        if (x0 >= 0 && y0 >= 0 && x0 < borders.first && y0 < borders.second) {
          data[y0][x0] = 0; // Set pixel to black
        }

        if (x0 == x1 && y0 == y1) break;
        e2 = err;
        if (e2 > -dx) { err -= dy; x0 += sx; }
        if (e2 <  dy) { err += dx; y0 += sy; }
      }
    }
  }

  return data;
}

void GraphVisualizator::RoundVertexCoordinates(std::vector<Point>* coordinates) {
  for (auto& [x, y] : *coordinates) {
    x = std::round(x);
    y = std::round(y);
  }
}

std::pair<int, int> GraphVisualizator::CorrectCoordinates(
    std::vector<Point>* coordinates) {
  int min_x = std::numeric_limits<int>::max();
  int min_y = std::numeric_limits<int>::max();
  
  int max_x = std::numeric_limits<int>::min();
  int max_y = std::numeric_limits<int>::min();
  
  // Iterate over vertices to find minimum and maximum coordinates
  for (const auto& [x, y] : *coordinates) {
    min_x = std::min(min_x, static_cast<int>(x));
    min_y = std::min(min_y, static_cast<int>(y));

    max_x = std::max(max_x, static_cast<int>(x));
    max_y = std::max(max_y, static_cast<int>(y));
  }

  // Calculate displacements to ensure all coordinates are positive
  int disp_x = kPadding_ - min_x;
  int disp_y = kPadding_ - min_y;

  // Update vertex coordinates to ensure positivity
  for (auto& [x, y] : *coordinates) {
    x += disp_x;
    y += disp_y;
  }

  // Return corrected coordinates with additional margin
  return {max_x + disp_x + kPadding_, max_y + disp_y + kPadding_};
}
