#include "header/vector_2d.h"

#include <cmath>
#include <utility>

#include "header/point_2d.h"

Point Vector::operator+(const Point& start) const {
  auto [x, y] = start.GetCoordinates();
  return Point(x + x_, y + y_);
}

Point Vector::operator*(const double& scalar) const {
  return Point(x_ * scalar, y_ * scalar);
}

double Vector::operator*(const Vector& other) const {
  auto [x, y] = other.GetCoordinates();
  return x_ * x + y_ * y;
}

double Vector::GetLength() const {
  return std::sqrt(operator*(*this));
}

std::pair<double, double> Vector::GetCoordinates() const {
  return {x_, y_};
}
