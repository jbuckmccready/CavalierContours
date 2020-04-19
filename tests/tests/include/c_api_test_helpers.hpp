#ifndef CAVC_API_TEST_HELPERS_HPP
#define CAVC_API_TEST_HELPERS_HPP
#include "cavaliercontours.h"
#include "testhelpers.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <iostream>

inline bool vertexesFuzzyEqual(cavc_vertex const &left, cavc_vertex const &right) {
  return fuzzyEqual(left.x, right.x) && fuzzyEqual(left.y, right.y) &&
         fuzzyEqual(left.bulge, right.bulge);
}

template <typename Container>
inline std::size_t nextWrappingIndex(Container const &container, std::size_t index) {
  if (index == container.size() - 1) {
    return 0;
  }
  return index + 1;
}

MATCHER(VertexFuzzyEqual, "") {
  auto const &left = std::get<0>(arg);
  auto const &right = std::get<1>(arg);
  return vertexesFuzzyEqual(left, right);
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

MATCHER_P(VertexListsFuzzyEqual, isClosed, "") {
  auto const &left = std::get<0>(arg);
  auto const &right = std::get<1>(arg);
  if (left.size() != right.size()) {
    *result_listener << "sizes of vertex lists do not match ";
    return false;
  }

  std::size_t const vertexCount = left.size();

  if (!isClosed) {
    // open polyline indexes much match up
    for (std::size_t i = 0; i < vertexCount; ++i) {
      if (!vertexesFuzzyEqual(left[i], right[i])) {
        *result_listener << "vertexes not equal at index: " << i << " ";
        return false;
      }
    }
    return true;
  }

  // vertexes may not have same indexes in the case of a closed polyline, find first matching and
  // start matching from there
  std::size_t startIndex = 0;
  for (; startIndex < vertexCount; ++startIndex) {
    if (vertexesFuzzyEqual(left[0], right[startIndex])) {
      break;
    }
  }

  if (startIndex == vertexCount) {
    *result_listener << "did not find matching vertex to start with ";
    return false;
  }

  *result_listener << " started matching at index: " << startIndex << " ";
  // first one already matched, start at next index
  std::size_t index = nextWrappingIndex(left, startIndex);
  for (std::size_t i = 1; i < vertexCount; ++i) {
    if (!vertexesFuzzyEqual(left[i], right[index])) {
      *result_listener << "vertexes not equal at index: " << index << " ";
      return false;
    }
    index = nextWrappingIndex(left, index);
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
