#include "header/bmp.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

void BmpPainter::Write(const std::string& filename) {
  std::ofstream image(filename, std::ios::binary);
  
  if (!image.is_open()) {
    std::cerr << "File wasn't opened";
    return;
  }

  // Write file header
  image.write(reinterpret_cast<char*>(&file_header_), sizeof(FileHeader));
  image.write(reinterpret_cast<char*>(&info_header_), sizeof(InfoHeader));
  
  // Write color table
  for (const RGB& color : color_table_) {
    image.write(reinterpret_cast<const char*>(&color), sizeof(RGB));
  }

  // Write image data
  image.write(reinterpret_cast<char*>(data_.data()), data_.size());

  image.close();
}

void BmpPainter::CreateColorTable() {
  for (int i = 0; i < 256; ++i) {
    color_table_[i].red = i;
    color_table_[i].green = i;
    color_table_[i].blue = i;
  }
}
