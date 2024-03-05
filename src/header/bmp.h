#ifndef GRAPHVISUALIZATION_HEADER_BMP_H_
#define GRAPHVISUALIZATION_HEADER_BMP_H_

#include <cstdint>
#include <string>
#include <vector>

#pragma pack(push, 1)
struct FileHeader {
  uint16_t file_type{ 0x4D42 };
  uint32_t file_size{ 54 };
  uint16_t reserved_1{ 0 };
  uint16_t reserved_2{ 0 };
  uint32_t offset_data{ 54 + 256 * 4 };
};

struct InfoHeader {
  uint32_t header_size{ 40 };
  int32_t width{ 0 };
  int32_t height{ 0 };
  uint16_t planes{ 1 };
  uint16_t bit_count{ 8 };
  uint32_t compression{ 0 };
  uint32_t size_image{ 0 };
  int32_t x_pixels_per_meter{ 0 };
  int32_t y_pixels_per_meter{ 0 };
  uint32_t colors_used{ 256 };
  uint32_t colors_important{ 256 };
};

struct RGB {
  uint8_t red;
  uint8_t green;
  uint8_t blue;
  uint8_t reserved = 0;
};
#pragma pack(pop)

class BmpPainter {
public:
  BmpPainter() { CreateColorTable(); };

  void Read(const std::string& filename);
  void Write(const std::string& filename);

  FileHeader file_header_;
  InfoHeader info_header_;
  std::vector<RGB> color_table_{ 256 };
  std::vector<uint8_t> data_;
  
private:
  void CreateColorTable();
};

#endif  // GRAPHVISUALIZATION_HEADER_BMP_H_