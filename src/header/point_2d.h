#ifndef GRAPHVISUALIZATION_HEADER_POINT_2D_H_
#define GRAPHVISUALIZATION_HEADER_POINT_2D_H_

#include <utility>

class Vector;

class Point {
 public:
  Point(double x = 0.0, double y = 0.0) : x_{x}, y_{y} {}

  Point operator+(const Vector& direction) const;
  Vector operator-(const Point& end) const;

  std::pair<double, double> GetCoordinates() const;

 private:
  double x_ = 0.0;
  double y_ = 0.0;
};

#endif  // GRAPHVISUALIZATION_HEADER_POINT_2D_H_