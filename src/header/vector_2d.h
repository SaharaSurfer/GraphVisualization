#ifndef GRAPHVISUALIZATION_HEADER_VECTOR_2D_H_
#define GRAPHVISUALIZATION_HEADER_VECTOR_2D_H_

#include <utility>

class Point;

class Vector {
 public:
  Vector(double x = 0.0, double y = 0.0) : x_{x}, y_{y} {}

  Point operator+(const Point& start) const;
  Point operator*(const double& scalar) const;

  // Scalar product
  double operator*(const Vector& other) const;
  double GetLength() const;

  std::pair<double, double> GetCoordinates() const;

 private:
  double x_ = 0.0;
  double y_ = 0.0;
};

#endif  // GRAPHVISUALIZATION_HEADER_VECTOR_2D_H_