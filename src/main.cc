#include <cstdint>
#include <vector>
#include <deque>
#include <iostream>
#include <fstream>

#include "header/bmp.h"
#include "header/graph_visualizator.h"

int main() {
  Graph graph;
  graph.ReadGraph("../src/graphs/star.txt");
  
  GraphVisualizator bob_ross;
  bob_ross.DrawBmp(graph, "image.bmp");

  return 0;
}