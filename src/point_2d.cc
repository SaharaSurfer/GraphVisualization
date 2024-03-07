#include "header/point_2d.h"

#include <utility>

#include "header/vector_2d.h"

Point Point::operator+(const Vector& direction) const {
  auto [x, y] = direction.GetCoordinates();
  return Point(x + x_, y + y_);
}

Vector Point::operator-(const Point& end) const {
  auto [x, y] = end.GetCoordinates();
  return Vector(x - x_ , y - y_);
}

std::pair<double, double> Point::GetCoordinates() const {
  return {x_, y_};
}
