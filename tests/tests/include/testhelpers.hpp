#ifndef CAVC_TESTHELPERS_HPP
#define CAVC_TESTHELPERS_HPP
#include <cmath>
#include "cavaliercontours.h"
constexpr inline cavc_real PI() { return 3.14159265358979323846264338327950288; }
constexpr inline cavc_real TEST_EPSILON() { return 1e-5; }
template <typename Real>
inline bool fuzzyEqual(Real const &left, cavc_real const &right) {
return std::abs(left - right) < TEST_EPSILON();
}

#endif // CAVC_TESTHELPERS_HPP
