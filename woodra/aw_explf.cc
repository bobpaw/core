#include <iostream>

#include <pcu_util.h> // PCU asserts

#include <ma.h> // MeshAdapt

#include "aw_explf.h"

ExpLengthFunc::ExpLengthFunc(ma::Mesh *mesh) : m(mesh) {
  // Find mesh bounds.
  ma::Vector mins(HUGE_VAL, HUGE_VAL, HUGE_VAL), maxs(-HUGE_VAL, -HUGE_VAL, -HUGE_VAL);
  ma::Iterator *it = m->begin(0);
  while (ma::Entity *ent = m->iterate(it)) {
    ma::Vector v = ma::getPosition(m, ent);

    if (v.x() < mins.x()) mins.x() = v.x();
    if (v.x() > maxs.x()) maxs.x() = v.x();
    if (v.y() < mins.y()) mins.y() = v.y();
    if (v.y() > maxs.y()) maxs.y() = v.y();
    if (v.z() < mins.z()) mins.z() = v.z();
    if (v.z() > maxs.z()) maxs.z() = v.z();
  }
  m->end(it);

  ma::Vector length = maxs - mins;
  if (length.x() >= length.y() && length.x() >= length.z()) {
    // X is the longest axis.
    axis = Axis::X;
    midpoint = (mins.x() + maxs.x()) / 2.0;
    h0 = length.x() / 20.0;
    std::cout << "Adapting along x-axis." << std::endl;
  } else if (length.y() >= length.x() && length.y() >= length.z()) {
    // Y is the longest axis.
    axis = Axis::Y;
    midpoint = (mins.y() + maxs.x()) / 2.0;
    h0 = length.y() / 20.0;
    std::cout << "Adapting along y-axis." << std::endl;
  } else {
    // Z is the longest axis.
    PCU_DEBUG_ASSERT(length.z() >= length.x() && length.z() >= length.y());
    axis = Axis::Z;
    midpoint = (mins.z() + maxs.z()) / 2.0;
    h0 = length.z() / 20.0;
    std::cout << "Adapting along z-axis." << std::endl;
  }
}

double ExpLengthFunc::ldist (const ma::Vector &v) const {
  switch (axis) {
    case Axis::X:
      return std::fabs(v.x() - midpoint);
    case Axis::Y:
      return std::fabs(v.y() - midpoint);
    case Axis::Z:
      return std::fabs(v.z() - midpoint);
  }
}

double ExpLengthFunc::eval (const ma::Vector &v) const {
  double dist = ldist(v);
  double h_min = h0 / 4.0, h_max = h0 * 2.0;
  double a = 0.1524, b = 2;
  return h_max + (h_min - h_max) * std::exp(-a*std::pow(ldist(v), b)/h0);
}

