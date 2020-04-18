#ifndef CAVC_API_TEST_HELPERS_HPP
#define CAVC_API_TEST_HELPERS_HPP
#include "cavaliercontours.h"
#include "testhelpers.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <iostream>

MATCHER(VertexFuzzyEqual, "") {
  auto const &left = std::get<0>(arg);
  auto const &right = std::get<1>(arg);
  return fuzzyEqual(left.x, right.x) && fuzzyEqual(left.y, right.y) &&
         fuzzyEqual(left.bulge, right.bulge);
}

MATCHER(VertexEqual, "") {
  auto const &left = std::get<0>(arg);
  auto const &right = std::get<1>(arg);
  return left.x == right.x && left.y == right.y && left.bulge == right.bulge;
}

MATCHER(PointFuzzyEqual, "") {
  auto const &left = std::get<0>(arg);
  auto const &right = std::get<1>(arg);
  return fuzzyEqual(left.x, right.x) && fuzzyEqual(left.y, right.y);
}

MATCHER(VertexListsFuzzyEqual, "") {
  auto const &left = std::get<0>(arg);
  auto const &right = std::get<1>(arg);
  if (left.size() != right.size()) {
    *result_listener << "sizes of vertex lists do not match";
    return false;
  }

  for (std::size_t i = 0; i < left.size(); ++i) {
    if (!ExplainMatchResult(VertexFuzzyEqual(), std::make_tuple(left[i], right[i]),
                            result_listener)) {
      *result_listener << "at parent index: " << i;
      return false;
    }
  }
  return true;
}

inline std::ostream &operator<<(std::ostream &os, cavc_vertex const &v) {
  os << '[' << v.x << "," << v.y << "," << v.bulge << ']';
  return os;
}

inline std::ostream &operator<<(std::ostream &os, cavc_point const &p) {
  os << '[' << p.x << "," << p.y << ']';
  return os;
}

#endif // CAVC_API_TEST_HELPERS_HPP
