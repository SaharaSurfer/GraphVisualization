#include <cstdint>
#include <vector>
#include <deque>
#include <iostream>
#include <fstream>

#include "header/bmp.h"
#include "header/graph_visualizator.h"

int main() {
  GraphVisualizator gr_bob_ross;
  gr_bob_ross.ReadGraph("../src/test.txt");
  std::vector<std::vector<uint8_t>> data = gr_bob_ross.GetData();

  BmpPainter bob_ross;

  std::vector<uint8_t> flattened_data;
  for (const auto& row : data) {
    for (const uint8_t& value : row) {
      flattened_data.push_back(value);
    }
  }

  bob_ross.data_ = flattened_data;

  FileHeader file_header;
  InfoHeader info_header;

  file_header.file_size += flattened_data.size();
  info_header.width = data[0].size();
  info_header.height = data.size();
  info_header.size_image = flattened_data.size();

  bob_ross.file_header_ = file_header;
  bob_ross.info_header_ = info_header;

  std::string filename = "image.bmp";
  bob_ross.Write(filename);

  return 0;
}