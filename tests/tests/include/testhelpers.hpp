#ifndef CAVC_TESTHELPERS_HPP
#define CAVC_TESTHELPERS_HPP
#include <cmath>
#include "cavaliercontours.h"
inline cavc_real PI() { return 3.14159265358979323846264338327950288; }
template<typename Real>
inline bool fuzzyEqual(Real const &left, cavc_real const &right, Real epsilon = 1e-5) {
  return std::abs(left - right) < epsilon;
}

#endif // CAVC_TESTHELPERS_HPP
