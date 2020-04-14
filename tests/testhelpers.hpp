#ifndef CAVC_TESTHELPERS_HPP
#define CAVC_TESTHELPERS_HPP
#include "cavaliercontours.h"
#include <cmath>
#define CAVC_PI 3.14159265358979323846264338327950288
inline bool fuzzyEqual(cavc_real const &left, cavc_real const &right, cavc_real epsilon = 1e-5) {
  return std::abs(left - right) < epsilon;
}

#endif // CAVC_TESTHELPERS_HPP
