#include <cstdint>
#include <vector>
#include <deque>
#include <iostream>
#include <fstream>

#include "header/bmp.h"
#include "header/graph_visualizator.h"

int main() {
  Bmp bob_ross;

  GraphVisualizator gr_bob_ross;
  gr_bob_ross.ReadGraph("../src/graphs/slipper.txt");
  DataMatrix data = gr_bob_ross.GetData();

  bob_ross.Interpret(data);

  std::string filename = "image.bmp";
  bob_ross.Write(filename);

  return 0;
}