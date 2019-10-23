#ifndef POLYLINE_H
#define POLYLINE_H
#include "intrcircle2circle2.h"
#include "intrlineseg2circle2.h"
#include "intrlineseg2lineseg2.h"
#include "staticspatialindex.h"
#include "vector2.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cavc {

template <typename Real> class PlineVertex {
public:
  PlineVertex() = default;
  PlineVertex(Real x, Real y, Real bulge) : m_position(x, y), m_bulge(bulge) {}
  PlineVertex(Vector2<Real> position, Real bulge)
      : PlineVertex(position.x(), position.y(), bulge) {}

  Real x() const { return m_position.x(); }
  Real &x() { return m_position.x(); }

  Real y() const { return m_position.y(); }
  Real &y() { return m_position.y(); }

  Real bulge() const { return m_bulge; }
  Real &bulge() { return m_bulge; }

  bool bulgeIsZero(Real epsilon = utils::realPrecision<Real>) const {
    return std::abs(m_bulge) < epsilon;
  }

  Vector2<Real> const &pos() const { return m_position; }
  Vector2<Real> &pos() { return m_position; }

private:
  Vector2<Real> m_position;
  Real m_bulge;
};

template <typename Real> class Polyline {
public:
  Polyline() : m_isClosed(false), m_vertexes() {}

  using PlineVertex = PlineVertex<Real>;

  inline PlineVertex const &operator[](std::size_t i) const { return m_vertexes[i]; }

  inline PlineVertex &operator[](std::size_t i) { return m_vertexes[i]; }

  bool isClosed() const { return m_isClosed; }
  bool &isClosed() { return m_isClosed; }

  void addVertex(Real x, Real y, Real bulge) { m_vertexes.emplace_back(x, y, bulge); }
  void addVertex(PlineVertex vertex) { addVertex(vertex.x(), vertex.y(), vertex.bulge()); }

  std::size_t size() const { return m_vertexes.size(); }

  PlineVertex const &lastVertex() const { return m_vertexes.back(); }
  PlineVertex &lastVertex() { return m_vertexes.back(); }

  std::vector<PlineVertex> &vertexes() { return m_vertexes; }
  std::vector<PlineVertex> const &vertexes() const { return m_vertexes; }

private:
  bool m_isClosed;
  std::vector<PlineVertex> m_vertexes;
};

// Polyline traversal/iteration functions.
template <typename Real, typename F> void iterateSegIndices(Polyline<Real> const &pline, F &&f) {
  std::size_t i;
  std::size_t j;
  if (pline.isClosed()) {
    i = 0;
    j = pline.vertexes().size() - 1;
  } else {
    i = 1;
    j = 0;
  }

  while (i < pline.vertexes().size() && f(j, i)) {
    j = i;
    i = i + 1;
  }
}

template <typename Real, typename F> void iterateSegs(Polyline<Real> const &pline, F &&f) {
  auto visitor = [&](std::size_t i, std::size_t j) {
    return std::forward<F>(f)(pline[i], pline[j]);
  };

  iterateSegIndices(pline, visitor);
}

template <typename Real, typename F> void iterateSegs(Polyline<Real> &pline, F &&f) {
  auto visitor = [&](std::size_t i, std::size_t j) {
    return std::forward<F>(f)(pline[i], pline[j]);
  };

  iterateSegIndices(pline, visitor);
}

template <typename Real> struct ArcRadiusAndCenter {
  Real radius;
  Vector2<Real> center;
};

/// Get the arc radius and arc center of a arc segment defined by v1 to v2.
template <typename Real>
ArcRadiusAndCenter<Real> arcRadiusAndCenter(PlineVertex<Real> const &v1,
                                            PlineVertex<Real> const &v2) {
  assert(!v1.bulgeIsZero() && "v1 to v2 must be an arc");
  assert(!fuzzyEqual(v1.pos(), v2.pos()) && "v1 must not equal v2");

  // compute radius
  Real b = std::abs(v1.bulge());
  Vector2<Real> v = v2.pos() - v1.pos();
  Real d = length(v);
  Real r = d * (b * b + Real(1)) / (Real(4) * b);

  // compute center
  Real s = b * d / Real(2);
  Real m = r - s;
  Real offsX = -m * v.y() / d;
  Real offsY = m * v.x() / d;
  if (v1.bulge() < Real(0)) {
    offsX = -offsX;
    offsY = -offsY;
  }

  Vector2<Real> c{v1.x() + v.x() / Real(2) + offsX, v1.y() + v.y() / Real(2) + offsY};
  return ArcRadiusAndCenter<Real>{r, c};
}

/// Result of splitting a segment v1 to v2.
template <typename Real> struct SplitResult {
  /// Updated starting vertex.
  PlineVertex<Real> updatedStart;
  /// Vertex at the split point.
  PlineVertex<Real> splitVertex;
};

/// Split the segment defined by v1 to v2 at some point defined along it.
template <typename Real>
SplitResult<Real> splitAtPoint(PlineVertex<Real> const &v1, PlineVertex<Real> const &v2,
                               Vector2<Real> const &point) {
  SplitResult<Real> result;
  if (v1.bulgeIsZero()) {
    result.updatedStart = v1;
    result.splitVertex = PlineVertex<Real>(point, Real(0));
  } else if (fuzzyEqual(v1.pos(), v2.pos())) {
    result.updatedStart = PlineVertex<Real>(point, Real(0));
    result.splitVertex = PlineVertex<Real>(point, v1.bulge());
  } else {
    auto radiusAndCenter = arcRadiusAndCenter(v1, v2);
    Vector2<Real> arcCenter = radiusAndCenter.center;
    Real a = angle(arcCenter, point);
    Real arcStartAngle = angle(arcCenter, v1.pos());
    Real theta1 = utils::deltaAngle(arcStartAngle, a);
    Real bulge1 = std::tan(theta1 / Real(4));
    Real arcEndAngle = angle(arcCenter, v2.pos());
    Real theta2 = utils::deltaAngle(a, arcEndAngle);
    Real bulge2 = std::tan(theta2 / Real(4));

    result.updatedStart = PlineVertex<Real>(v1.pos(), bulge1);
    result.splitVertex = PlineVertex<Real>(point, bulge2);
  }

  return result;
}

/// Axis aligned bounding box (AABB).
template <typename Real> struct AABB {
  Real xMin;
  Real yMin;
  Real xMax;
  Real yMax;

  void expand(Real val) {
    xMin -= val;
    yMin -= val;
    xMax += val;
    yMax += val;
  }
};

/// Compute the extents of a polyline.
template <typename Real> AABB<Real> extents(Polyline<Real> const &pline) {
  if (pline.vertexes().size() < 2) {
    return AABB<Real>{pline[0].x(), pline[0].y(), pline[0].x(), pline[0].y()};
  }

  AABB<Real> result{std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max(),
                    std::numeric_limits<Real>::min(), std::numeric_limits<Real>::min()};

  auto visitor = [&](std::size_t i, std::size_t j) {
    PlineVertex<Real> const &v1 = pline[i];
    if (v1.bulgeIsZero()) {
      if (v1.x() < result.xMin)
        result.xMin = v1.x();
      if (v1.y() < result.yMin)
        result.yMin = v1.y();
      if (v1.x() > result.yMax)
        result.xMax = v1.x();
      if (v1.y() > result.yMax)
        result.yMax = v1.y();
    } else {
      PlineVertex<Real> const &v2 = pline[j];
      auto arc = arcRadiusAndCenter(v1, v2);

      Real startAngle = angle(arc.center, v1.pos());
      Real endAngle = angle(arc.center, v2.pos());
      Real sweepAngle = utils::deltaAngle(startAngle, endAngle);

      Real arcXMin, arcYMin, arcXMax, arcYMax;

      // crosses PI/2
      if (utils::angleIsWithinSweep(startAngle, sweepAngle, Real(0.5) * utils::pi<Real>)) {
        arcYMax = arc.center.y() + arc.radius;
      } else {
        arcYMax = std::max(v1.y(), v2.y());
      }

      // crosses PI
      if (utils::angleIsWithinSweep(startAngle, sweepAngle, utils::pi<Real>)) {
        arcXMin = arc.center.x() - arc.radius;
      } else {
        arcXMin = std::min(v1.x(), v2.x());
      }

      // crosses 3PI/2
      if (utils::angleIsWithinSweep(startAngle, sweepAngle, Real(1.5) * utils::pi<Real>)) {
        arcYMin = arc.center.y() - arc.radius;
      } else {
        arcYMin = std::min(v1.y(), v2.y());
      }

      // crosses 2PI
      if (utils::angleIsWithinSweep(startAngle, sweepAngle, Real(2) * utils::pi<Real>)) {
        arcXMax = arc.center.x() + arc.radius;
      } else {
        arcXMax = std::max(v1.x(), v2.x());
      }

      if (arcXMin < result.xMin)
        result.xMin = arcXMin;
      if (arcYMin < result.yMin)
        result.yMin = arcYMin;
      if (arcXMax > result.yMax)
        result.xMax = arcXMax;
      if (arcYMax > result.yMax)
        result.yMax = arcYMax;
    }

    // return true to iterate all segments
    return true;
  };

  iterateSegIndices(pline, visitor);

  return result;
}

template <typename Real>
Vector2<Real> closestPointOnSeg(PlineVertex<Real> const &v1, PlineVertex<Real> const &v2,
                                Vector2<Real> const &point) {
  if (v1.bulgeIsZero()) {
    return closestPointOnLineSeg(v1.pos(), v2.pos(), point);
  }

  auto arc = arcRadiusAndCenter(v1, v2);

  if (fuzzyEqual(point, arc.center)) {
    // avoid normalizing zero length vector (point is at center, just return start point)
    return v1.pos();
  }

  if (pointWithinArcSweepAngle(arc.center, v1.pos(), v2.pos(), v1.bulge(), point)) {
    // closest point is on the arc
    Vector2<Real> vToPoint = point - arc.center;
    normalize(vToPoint);
    return arc.radius * vToPoint + arc.center;
  }

  // else closest point is one of the ends
  Real dist1 = distSquared(v1.pos(), point);
  Real dist2 = distSquared(v2.pos(), point);
  if (dist1 < dist2) {
    return v1.pos();
  }

  return v2.pos();
}

template <typename Real> struct ClosestPointAndIndex {
  std::size_t index;
  Vector2<Real> point;
  Real distance;
};

template <typename Real>
ClosestPointAndIndex<Real> closestPointAndIndex(Polyline<Real> const &pline,
                                                Vector2<Real> const &point) {
  assert(pline.vertexes().size > 0 && "empty polyline has no closest point");
  ClosestPointAndIndex<Real> result;
  if (pline.vertexes().size() == 1) {
    result.index = 0;
    result.point = pline[0];
    result.distance = length(point - pline[0].pos());
    return result;
  }

  auto visitor = [&](std::size_t i, std::size_t j) {
    Vector2<Real> cp = closestPointOnSeg(pline[i], pline[j], point);
    Real dist = length(point - cp);
    if (dist < result.distance) {
      result.index = i;
      result.point = cp;
      result.distance = dist;
    }

    // iterate all segments
    return true;
  };

  iterateSegIndices(pline, visitor);

  return result;
}

/// Computes a fast approximate AABB of a segment described by v1 to v2, bounding box may be larger
/// than the true bounding box for the segment
template <typename Real>
AABB<Real> createFastApproxBoundingBox(PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
  AABB<Real> result;
  if (v1.bulgeIsZero()) {
    if (v1.x() < v2.x()) {
      result.xMin = v1.x();
      result.xMax = v2.x();
    } else {
      result.xMin = v2.x();
      result.xMax = v1.x();
    }

    if (v1.y() < v2.y()) {
      result.yMin = v1.y();
      result.yMax = v2.y();
    } else {
      result.yMin = v2.y();
      result.yMax = v1.y();
    }

    return result;
  }

  // For arcs we don't compute the actual extents which requires slow trig functions, instead we
  // create an approximate bounding box from the rectangle formed by extending the chord by the
  // sagitta, NOTE: this approximate bounding box is always equal to or bigger than the true
  // bounding box
  Real b = v1.bulge();
  Real offsX = b * (v2.y() - v1.y()) / Real(2);
  Real offsY = -b * (v2.x() - v1.x()) / Real(2);

  Real pt1X = v1.x() + offsX;
  Real pt2X = v2.x() + offsX;
  Real pt1Y = v1.y() + offsY;
  Real pt2Y = v2.y() + offsY;

  Real endPointXMin, endPointXMax;
  if (v1.x() < v2.x()) {
    endPointXMin = v1.x();
    endPointXMax = v2.x();
  } else {
    endPointXMin = v2.x();
    endPointXMax = v1.x();
  }

  Real ptXMin, ptXMax;
  if (pt1X < pt2X) {
    ptXMin = pt1X;
    ptXMax = pt2X;
  } else {
    ptXMin = pt2X;
    ptXMax = pt1X;
  }

  Real endPointYMin, endPointYMax;
  if (v1.y() < v2.y()) {
    endPointYMin = v1.y();
    endPointYMax = v2.y();
  } else {
    endPointYMin = v2.y();
    endPointYMax = v1.y();
  }

  Real ptYMin, ptYMax;
  if (pt1Y < pt2Y) {
    ptYMin = pt1Y;
    ptYMax = pt2Y;
  } else {
    ptYMin = pt2Y;
    ptYMax = pt1Y;
  }

  result.xMin = std::min(endPointXMin, ptXMin);
  result.yMin = std::min(endPointYMin, ptYMin);
  result.xMax = std::max(endPointXMax, ptXMax);
  result.yMax = std::max(endPointYMax, ptYMax);
  return result;
}

template <typename Real>
StaticSpatialIndex<Real> createApproxSpatialIndex(Polyline<Real> const &pline) {
  assert(pline.size() > 1 && "need at least 2 vertexes to form segments for spatial index");

  std::size_t segmentCount = pline.isClosed() ? pline.size() : pline.size() - 1;
  StaticSpatialIndex<Real> result(segmentCount);

  for (std::size_t i = 0; i < pline.size() - 1; ++i) {
    AABB<Real> approxBB = createFastApproxBoundingBox(pline[i], pline[i + 1]);
    result.add(approxBB.xMin, approxBB.yMin, approxBB.xMax, approxBB.yMax);
  }

  if (pline.isClosed()) {
    // add final segment from last to first
    AABB<Real> approxBB = createFastApproxBoundingBox(pline.lastVertex(), pline[0]);
    result.add(approxBB.xMin, approxBB.yMin, approxBB.xMax, approxBB.yMax);
  }

  result.finish();

  return result;
}

template <typename Real> void invertDirection(Polyline<Real> &pline) {
  if (pline.size() < 2) {
    return;
  }
  std::reverse(std::begin(pline.vertexes()), std::end(pline.vertexes()));

  // shift and negate bulge (to maintain same geometric path)
  Real firstBulge = pline[0].bulge();

  for (std::size_t i = 1; i < pline.size(); ++i) {
    pline[i - 1].bulge() = -pline[i].bulge();
  }

  pline.lastVertex().bulge() = -firstBulge;
}

template <typename Real>
Polyline<Real> pruneSingularities(Polyline<Real> const &pline,
                                  Real epsilon = utils::realPrecision<Real>) {
  Polyline<Real> result;
  result.isClosed() = pline.isClosed();

  if (pline.size() == 0) {
    return result;
  }

  // allocate up front (most of the time the number of repeated positions 0 or much less than the
  // total number of vertexes so we're not using very much more memory than required)
  result.vertexes().reserve(pline.size());

  result.addVertex(pline[0]);

  for (std::size_t i = 1; i < pline.size(); ++i) {
    if (fuzzyEqual(result.lastVertex().pos(), pline[i].pos(), epsilon)) {
      result.lastVertex().bulge() = pline[i].bulge();
    } else {
      result.addVertex(pline[i]);
    }
  }

  if (result.isClosed() && result.size() > 1) {
    if (fuzzyEqual(result.lastVertex().pos(), result[0].pos(), epsilon)) {
      result.vertexes().pop_back();
    }
  }

  return result;
}

/// Compute the area of a closed polyline, assumes no self intersects, returns positive number if
/// polyline direction is counter clockwise, negative if clockwise, zero if not closed
template <typename Real> Real area(Polyline<Real> const &pline) {
  // Implementation notes:
  // Using the shoelace formula (https://en.wikipedia.org/wiki/Shoelace_formula) modified to support
  // arcs defined by a bulge value. The shoelace formula returns a negative value for clockwise
  // oriented polygons and positive value for counter clockwise oriented polygons. The area of each
  // circular segment defined by arcs is then added if it is a counter clockwise arc or subtracted
  // if it is a clockwise arc. The area of the circular segments are computed by finding the area of
  // the arc sector minus the area of the triangle defined by the chord and center of circle.
  // See https://en.wikipedia.org/wiki/Circular_segment
  if (!pline.isClosed() || pline.size() < 2) {
    return Real(0);
  }

  Real doubleEdgeAreaTotal = Real(0);
  Real doubleArcSegAreaTotal = Real(0);

  auto visitor = [&](std::size_t i, std::size_t j) {
    doubleEdgeAreaTotal += pline[i].x() * pline[j].y() - pline[i].y() * pline[j].x();
    if (!pline[i].bulgeIsZero()) {
      // add segment area
      Real b = std::abs(pline[i].bulge());
      Real sweepAngle = Real(4) * std::atan(b);
      Real triangleBase = length(pline[j].pos() - pline[i].pos());
      Real radius = triangleBase * (b * b + Real(1)) / (Real(4) * b);
      Real sagitta = b * triangleBase / Real(2);
      Real triangleHeight = radius - sagitta;
      Real doubleSectorArea = sweepAngle * radius * radius;
      Real doubleTriangleArea = triangleBase * triangleHeight;
      Real doubleArcSegArea = doubleSectorArea - doubleTriangleArea;
      if (pline[i].bulge() < Real(0)) {
        doubleArcSegArea = -doubleArcSegArea;
      }

      doubleArcSegAreaTotal += doubleArcSegArea;
    }

    // iterate all segments
    return true;
  };

  iterateSegIndices(pline, visitor);

  return (doubleEdgeAreaTotal + doubleArcSegAreaTotal) / Real(2);
}

template <typename Real> struct PlineOffsetSegment {
  PlineVertex<Real> v1;
  PlineVertex<Real> v2;
  Vector2<Real> origV2Pos;
  bool collapsedArc;
};

template <typename Real>
std::vector<PlineOffsetSegment<Real>> createUntrimmedOffsetSegments(Polyline<Real> const &pline,
                                                                    Real offset) {
  std::size_t segmentCount = pline.isClosed() ? pline.size() : pline.size() - 1;

  std::vector<PlineOffsetSegment<Real>> result;
  result.reserve(segmentCount);

  auto lineVisitor = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
    result.emplace_back();
    PlineOffsetSegment<Real> &seg = result.back();
    seg.collapsedArc = false;
    seg.origV2Pos = v2.pos();
    Vector2<Real> edge = v2.pos() - v1.pos();
    Vector2<Real> offsetV = offset * unitPerp(edge);
    seg.v1.pos() = v1.pos() + offsetV;
    seg.v1.bulge() = v1.bulge();
    seg.v2.pos() = v2.pos() + offsetV;
    seg.v2.bulge() = v2.bulge();
  };

  auto arcVisitor = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
    auto arc = arcRadiusAndCenter(v1, v2);
    Real offs = v1.bulge() < Real(0) ? offset : -offset;
    Real radiusAfterOffset = arc.radius + offs;
    Vector2<Real> v1ToCenter = v1.pos() - arc.center;
    normalize(v1ToCenter);
    Vector2<Real> v2ToCenter = v2.pos() - arc.center;
    normalize(v2ToCenter);

    result.emplace_back();
    PlineOffsetSegment<Real> &seg = result.back();
    seg.origV2Pos = v2.pos();
    seg.v1.pos() = offs * v1ToCenter + v1.pos();
    seg.v2.pos() = offs * v2ToCenter + v2.pos();
    seg.v2.bulge() = v2.bulge();

    if (radiusAfterOffset < utils::realThreshold<Real>) {
      // collapsed arc, offset arc start and end points towards arc center and turn into line
      // handles case where offset vertexes are equal and simplifies path for clipping algorithm
      seg.collapsedArc = true;
      seg.v1.bulge() = Real(0);
    } else {
      seg.collapsedArc = false;
      seg.v1.bulge() = v1.bulge();
    }
  };

  auto offsetVisitor = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
    if (v1.bulgeIsZero()) {
      lineVisitor(v1, v2);
    } else {
      arcVisitor(v1, v2);
    }
  };

  for (std::size_t i = 1; i < pline.size(); ++i) {
    offsetVisitor(pline[i - 1], pline[i]);
  }

  if (pline.isClosed()) {
    offsetVisitor(pline.lastVertex(), pline[0]);
  }

  return result;
}

namespace detail {
struct IndexPairHash {
  std::size_t operator()(std::pair<std::size_t, std::size_t> const &pair) const {
    return std::hash<std::size_t>()(pair.first) ^ std::hash<std::size_t>()(pair.second);
  }
};

template <typename Real> bool falseIntersect(Real t) { return t < 0.0 || t > 1.0; }

template <typename Real> Real segX(Vector2<Real> const &v1, Vector2<Real> const &v2, Real t) {
  return v1.x() + t * (v2.x() - v1.x());
}
template <typename Real> Real segY(Vector2<Real> const &v1, Vector2<Real> const &v2, Real t) {
  return v1.y() + t * (v2.y() - v1.y());
}

template <typename Real>
Vector2<Real> segMidpoint(PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
  if (v1.bulgeIsZero()) {
    return midpoint(v1.pos(), v2.pos());
  }

  auto arc = arcRadiusAndCenter(v1, v2);
  Real a1 = angle(arc.center, v1.pos());
  Real a2 = angle(arc.center, v2.pos());
  Real midAngle = a1 + utils::deltaAngle(a1, a2) / Real(2);
  return pointOnCircle(arc.radius, arc.center, midAngle);
}

template <typename Real>
void addOrReplaceIfSamePos(Polyline<Real> &pline, PlineVertex<Real> const &vertex,
                           Real epsilon = utils::realPrecision<Real>) {
  if (pline.size() == 0) {
    pline.addVertex(vertex);
    return;
  }

  if (fuzzyEqual(pline.lastVertex().pos(), vertex.pos(), epsilon)) {
    pline.lastVertex().bulge() = vertex.bulge();
    return;
  }

  pline.addVertex(vertex);
}
// Gets the bulge to describe the arc going from start point to end point with the given arc center
// and curve orientation, if orientation is negative then bulge is negative otherwise it is positive
template <typename Real>
Real bulgeForConnection(Vector2<Real> const &arcCenter, Vector2<Real> const &sp,
                        Vector2<Real> const &ep, int orientation) {
  Real a1 = angle(arcCenter, sp);
  Real a2 = angle(arcCenter, ep);
  Real absSweepAngle = std::abs(utils::deltaAngle(a1, a2));
  Real absBulge = std::tan(absSweepAngle / Real(4));
  if (orientation < 0) {
    return -absBulge;
  }

  return absBulge;
}

template <typename Real>
void lineToLineJoin(PlineOffsetSegment<Real> const &s1, PlineOffsetSegment<Real> const &s2,
                    int connectionArcOrientation, Polyline<Real> &result) {
  const auto &v1 = s1.v1;
  const auto &v2 = s1.v2;
  const auto &u1 = s2.v1;
  const auto &u2 = s2.v2;
  assert(v1.bulgeIsZero() && u1.bulgeIsZero() && "both pairs should be lines");

  auto connectUsingArc = [&] {
    auto const &arcCenter = s1.origV2Pos;
    auto const &sp = v2.pos();
    auto const &ep = u1.pos();
    Real bulge = bulgeForConnection(arcCenter, sp, ep, connectionArcOrientation);
    addOrReplaceIfSamePos(result, PlineVertex<Real>(sp, bulge));
    addOrReplaceIfSamePos(result, PlineVertex<Real>(ep, Real(0)));
  };

  if (s2.collapsedArc) {
    // connecting to collapsed arc, always connect using arc
    connectUsingArc();
  } else {
    auto intrResult = intrLineSeg2LineSeg2(v1.pos(), v2.pos(), u1.pos(), u2.pos());

    switch (intrResult.intrType) {
    case LineSeg2LineSeg2IntrType::None:
      // just join with straight line
      addOrReplaceIfSamePos(result, PlineVertex<Real>(v2.pos(), Real(0)));
      addOrReplaceIfSamePos(result, u1);
      break;
    case LineSeg2LineSeg2IntrType::True:
      addOrReplaceIfSamePos(result, PlineVertex<Real>(intrResult.point, Real(0)));
      break;
    case LineSeg2LineSeg2IntrType::Coincident:
      addOrReplaceIfSamePos(result, PlineVertex<Real>(v2.pos(), Real(0)));
      break;
    case LineSeg2LineSeg2IntrType::False:
      if (intrResult.t0 > Real(1) && falseIntersect(intrResult.t1)) {
        // extend and join the lines together using an arc
        connectUsingArc();
      } else {
        addOrReplaceIfSamePos(result, PlineVertex<Real>(v2.pos(), Real(0)));
        addOrReplaceIfSamePos(result, u1);
      }
      break;
    }
  }
}

template <typename Real>
void lineToArcJoin(PlineOffsetSegment<Real> const &s1, PlineOffsetSegment<Real> const &s2,
                   int connectionArcOrientation, Polyline<Real> &result) {

  const auto &v1 = s1.v1;
  const auto &v2 = s1.v2;
  const auto &u1 = s2.v1;
  const auto &u2 = s2.v2;
  assert(v1.bulgeIsZero() && !u1.bulgeIsZero() &&
         "first pair should be arc, second pair should be line");

  auto connectUsingArc = [&] {
    auto const &arcCenter = s1.origV2Pos;
    auto const &sp = v2.pos();
    auto const &ep = u1.pos();
    Real bulge = bulgeForConnection(arcCenter, sp, ep, connectionArcOrientation);
    addOrReplaceIfSamePos(result, PlineVertex<Real>(sp, bulge));
    addOrReplaceIfSamePos(result, u1);
  };

  const auto arc = arcRadiusAndCenter(u1, u2);

  auto processIntersect = [&](Real t, Vector2<Real> const &intersect) {
    const bool trueSegIntersect = !falseIntersect(t);
    const bool trueArcIntersect =
        pointWithinArcSweepAngle(arc.center, u1.pos(), u2.pos(), u1.bulge(), intersect);
    if (trueSegIntersect && trueArcIntersect) {
      // trim at intersect
      Real a = angle(arc.center, intersect);
      Real arcEndAngle = angle(arc.center, u2.pos());
      Real theta = utils::deltaAngle(a, arcEndAngle);
      Real bulge = std::tan(theta / Real(4));
      // ensure the sign matches (may get flipped if intersect is at the very end of the arc, in
      // which case we do not want to update the bulge)
      if ((bulge < Real(0)) == (u1.bulge() < Real(0))) {
        addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, bulge));
      } else {
        addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, u1.bulge()));
      }
    } else if (t > Real(1) && !trueArcIntersect) {
      connectUsingArc();
    } else if (s1.collapsedArc) {
      // collapsed arc connecting to arc, connect using arc
      connectUsingArc();
    } else {
      // connect using line
      addOrReplaceIfSamePos(result, PlineVertex<Real>(v2.pos(), Real(0)));
      addOrReplaceIfSamePos(result, u1);
    }
  };

  auto intrResult = intrLineSeg2Circle2(v1.pos(), v2.pos(), arc.radius, arc.center);
  if (intrResult.numIntersects == 0) {
    connectUsingArc();
  } else if (intrResult.numIntersects == 1) {
    Real x = segX(v1.pos(), v2.pos(), intrResult.t0);
    Real y = segY(v1.pos(), v2.pos(), intrResult.t0);
    processIntersect(intrResult.t0, Vector2<Real>(x, y));
  } else {
    assert(intrResult.numIntersects == 2);
    // always use intersect closest to original point
    Real x1 = segX(v1.pos(), v2.pos(), intrResult.t0);
    Real y1 = segY(v1.pos(), v2.pos(), intrResult.t0);
    Vector2<Real> i1 = Vector2<Real>(x1, y1);
    Real dist1 = distSquared(i1, s1.origV2Pos);
    Real x2 = segX(v1.pos(), v2.pos(), intrResult.t1);
    Real y2 = segY(v1.pos(), v2.pos(), intrResult.t1);
    Vector2<Real> i2 = Vector2<Real>(x2, y2);
    Real dist2 = distSquared(i2, s1.origV2Pos);

    if (dist1 < dist2) {
      processIntersect(intrResult.t0, i1);
    } else {
      processIntersect(intrResult.t1, i2);
    }
  }
}

template <typename Real>
void arcToLineJoin(PlineOffsetSegment<Real> const &s1, PlineOffsetSegment<Real> const &s2,
                   int connectionArcOrientation, Polyline<Real> &result) {

  const auto &v1 = s1.v1;
  const auto &v2 = s1.v2;
  const auto &u1 = s2.v1;
  const auto &u2 = s2.v2;
  assert(!v1.bulgeIsZero() && u1.bulgeIsZero() &&
         "first pair should be line, second pair should be arc");

  auto connectUsingArc = [&] {
    auto const &arcCenter = s1.origV2Pos;
    auto const &sp = v2.pos();
    auto const &ep = u1.pos();
    Real bulge = bulgeForConnection(arcCenter, sp, ep, connectionArcOrientation);
    addOrReplaceIfSamePos(result, PlineVertex<Real>(sp, bulge));
    addOrReplaceIfSamePos(result, u1);
  };

  const auto arc = arcRadiusAndCenter(v1, v2);

  auto processIntersect = [&](Real t, Vector2<Real> const &intersect) {
    if (utils::fuzzyEqual(t, Real(1))) {
      // NOTE: this fuzzy equals test with 1.0 is necessary to prevent the trim from flipping the
      // sign of the the arc at the end point/boundary case just add the vertex u2 (it's a true
      // intersect at the end of the line described by u1->u2 and next segment is a line so maintain
      // 0 bulge using u2)
      addOrReplaceIfSamePos(result, u2);
    } else {
      const bool trueSegIntersect = !falseIntersect(t);
      const bool trueArcIntersect =
          pointWithinArcSweepAngle(arc.center, v1.pos(), v2.pos(), v1.bulge(), intersect);
      if (trueSegIntersect && trueArcIntersect) {
        // modify previous bulge and trim at intersect
        PlineVertex<Real> &prevVertex = result.lastVertex();
        const bool updatePrevBulge =
            !prevVertex.bulgeIsZero() && !fuzzyEqual(prevVertex.pos(), v2.pos());

        if (updatePrevBulge) {
          Real a = angle(arc.center, intersect);
          auto prevArc = arcRadiusAndCenter(prevVertex, v2);
          Real prevArcStartAngle = angle(prevArc.center, prevVertex.pos());
          Real updatedPrevTheta = utils::deltaAngle(prevArcStartAngle, a);
          result.lastVertex().bulge() = std::tan(updatedPrevTheta / Real(4));
        }

        addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, Real(0)));

      } else {
        connectUsingArc();
      }
    }
  };

  auto intrResult = intrLineSeg2Circle2(u1.pos(), u2.pos(), arc.radius, arc.center);
  if (intrResult.numIntersects == 0) {
    connectUsingArc();
  } else if (intrResult.numIntersects == 1) {
    Real x = segX(u1.pos(), u2.pos(), intrResult.t0);
    Real y = segY(u1.pos(), u2.pos(), intrResult.t0);
    processIntersect(intrResult.t0, Vector2<Real>(x, y));
  } else {
    assert(intrResult.numIntersects == 2);
    const auto &origPoint = s2.collapsedArc ? u1.pos() : s1.origV2Pos;
    Real x1 = segX(u1.pos(), u2.pos(), intrResult.t0);
    Real y1 = segY(u1.pos(), u2.pos(), intrResult.t0);
    Vector2<Real> i1 = Vector2<Real>(x1, y1);
    Real dist1 = distSquared(i1, origPoint);
    Real x2 = segX(u1.pos(), u2.pos(), intrResult.t1);
    Real y2 = segY(u1.pos(), u2.pos(), intrResult.t1);
    Vector2<Real> i2 = Vector2<Real>(x2, y2);
    Real dist2 = distSquared(i2, origPoint);

    if (dist1 < dist2) {
      processIntersect(intrResult.t0, i1);
    } else {
      processIntersect(intrResult.t1, i2);
    }
  }
}

template <typename Real>
void arcToArcJoin(PlineOffsetSegment<Real> const &s1, PlineOffsetSegment<Real> const &s2,
                  int connectionArcOrientation, Polyline<Real> &result) {

  const auto &v1 = s1.v1;
  const auto &v2 = s1.v2;
  const auto &u1 = s2.v1;
  const auto &u2 = s2.v2;
  assert(!v1.bulgeIsZero() && !u1.bulgeIsZero() && "both pairs should be arcs");

  auto connectUsingArc = [&] {
    auto const &arcCenter = s1.origV2Pos;
    auto const &sp = v2.pos();
    auto const &ep = u1.pos();
    Real bulge = bulgeForConnection(arcCenter, sp, ep, connectionArcOrientation);
    addOrReplaceIfSamePos(result, PlineVertex<Real>(sp, bulge));
    addOrReplaceIfSamePos(result, u1);
  };

  const auto arc1 = arcRadiusAndCenter(v1, v2);
  const auto arc2 = arcRadiusAndCenter(u1, u2);

  auto processIntersect = [&](Vector2<Real> const &intersect) {
    if (fuzzyEqual(intersect, u1.pos())) {
      // NOTE: this fuzzy equals test with end point is necessary to prevent the trim from flipping
      // the sign of the the arc at the end point/boundary case
      addOrReplaceIfSamePos(result, u1);
    } else {
      const bool trueArcIntersect1 =
          pointWithinArcSweepAngle(arc1.center, v1.pos(), v2.pos(), v1.bulge(), intersect);
      const bool trueArcIntersect2 =
          pointWithinArcSweepAngle(arc2.center, u1.pos(), u2.pos(), u1.bulge(), intersect);

      if (trueArcIntersect1 && trueArcIntersect2) {
        // modify previous bulge and trim at intersect
        PlineVertex<Real> &prevVertex = result.lastVertex();
        const bool updatePrevBulge =
            !prevVertex.bulgeIsZero() && !fuzzyEqual(prevVertex.pos(), v2.pos());

        if (updatePrevBulge) {
          Real a = angle(arc1.center, intersect);
          auto prevArc = arcRadiusAndCenter(prevVertex, v2);
          Real prevArcStartAngle = angle(prevArc.center, prevVertex.pos());
          Real updatedPrevTheta = utils::deltaAngle(prevArcStartAngle, a);
          result.lastVertex().bulge() = std::tan(updatedPrevTheta / Real(4));
        }

        // add the vertex at our current trim/join point
        Real a = angle(arc2.center, intersect);
        Real endAngle = angle(arc2.center, u2.pos());
        Real theta = utils::deltaAngle(a, endAngle);
        addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, std::tan(theta / Real(4))));

      } else {
        connectUsingArc();
      }
    }
  };

  const auto intrResult = intrCircle2Circle2(arc1.radius, arc1.center, arc2.radius, arc2.center);
  switch (intrResult.intrType) {
  case Circle2Circle2IntrType::NoIntersect:
    connectUsingArc();
    break;
  case Circle2Circle2IntrType::OneIntersect:
    processIntersect(intrResult.point1);
    break;
  case Circle2Circle2IntrType::TwoIntersects: {
    Real dist1 = distSquared(intrResult.point1, s1.origV2Pos);
    Real dist2 = distSquared(intrResult.point2, s1.origV2Pos);
    if (dist1 < dist2) {
      processIntersect(intrResult.point1);
    } else {
      processIntersect(intrResult.point2);
    }
  } break;
  case Circle2Circle2IntrType::Coincident:
    // same constant arc radius and center, just add the vertex (nothing to trim/extend)
    addOrReplaceIfSamePos(result, u1);
    break;
  }
}

enum class PlineSegIntrType {
  NoIntersect,
  TangentIntersect,
  OneIntersect,
  TwoIntersects,
  SegmentOverlap,
  ArcOverlap
};

template <typename Real> struct IntrPlineSegsResult {
  PlineSegIntrType intrType;
  Vector2<Real> point1;
  Vector2<Real> point2;
};

template <typename Real>
IntrPlineSegsResult<Real> intrPlineSegs(PlineVertex<Real> const &v1, PlineVertex<Real> const &v2,
                                        PlineVertex<Real> const &u1, PlineVertex<Real> const &u2) {
  IntrPlineSegsResult<Real> result;
  const bool vIsLine = v1.bulgeIsZero();
  const bool uIsLine = u1.bulgeIsZero();

  // helper function to process line arc intersect
  auto processLineArcIntr = [&result](Vector2<Real> const &p0, Vector2<Real> const &p1,
                                      PlineVertex<Real> const &a1, PlineVertex<Real> const &a2) {

    auto arc = arcRadiusAndCenter(a1, a2);
    auto intrResult = intrLineSeg2Circle2(p0, p1, arc.radius, arc.center);

    // helper function to test and get point within arc sweep
    auto pointInSweep = [&](Real t) {
      if (t + utils::realThreshold<Real> < Real(0) || t > Real(1) + utils::realThreshold<Real>) {
        return std::make_pair(false, Vector2<Real>());
      }

      Vector2<Real> p(segX(p0, p1, t), segY(p0, p1, t));
      bool withinSweep = pointWithinArcSweepAngle(arc.center, a1.pos(), a2.pos(), a1.bulge(), p);
      return std::make_pair(withinSweep, p);
    };

    if (intrResult.numIntersects == 0) {
      result.intrType = PlineSegIntrType::NoIntersect;
    } else if (intrResult.numIntersects == 1) {
      auto p = pointInSweep(intrResult.t0);
      if (p.first) {
        result.intrType = PlineSegIntrType::OneIntersect;
        result.point1 = p.second;
      } else {
        result.intrType = PlineSegIntrType::NoIntersect;
      }
    } else {
      assert(intrResult.numIntersects == 2);
      auto p1 = pointInSweep(intrResult.t0);
      auto p2 = pointInSweep(intrResult.t1);

      if (p1.first && p2.first) {
        result.intrType = PlineSegIntrType::TwoIntersects;
        result.point1 = p1.second;
        result.point2 = p2.second;
      } else if (p1.first) {
        result.intrType = PlineSegIntrType::OneIntersect;
        result.point1 = p1.second;
      } else if (p2.first) {
        result.intrType = PlineSegIntrType::OneIntersect;
        result.point1 = p2.second;
      } else {
        result.intrType = PlineSegIntrType::NoIntersect;
      }
    }
  };

  if (vIsLine && uIsLine) {
    auto intrResult = intrLineSeg2LineSeg2(v1.pos(), v2.pos(), u1.pos(), u2.pos());
    switch (intrResult.intrType) {
    case LineSeg2LineSeg2IntrType::None:
      result.intrType = PlineSegIntrType::NoIntersect;
      break;
    case LineSeg2LineSeg2IntrType::True:
      result.intrType = PlineSegIntrType::OneIntersect;
      result.point1 = intrResult.point;
      break;
    case LineSeg2LineSeg2IntrType::Coincident:
      result.intrType = PlineSegIntrType::SegmentOverlap;
      // build points from parametric parameters for v1->v2
      result.point1 = Vector2<Real>(segX(v1.pos(), v2.pos(), intrResult.t0),
                                    segY(v1.pos(), v2.pos(), intrResult.t0));
      result.point2 = Vector2<Real>(segX(v1.pos(), v2.pos(), intrResult.t1),
                                    segY(v1.pos(), v2.pos(), intrResult.t1));
      break;
    case LineSeg2LineSeg2IntrType::False:
      result.intrType = PlineSegIntrType::NoIntersect;
      break;
    }

  } else if (vIsLine) {
    processLineArcIntr(v1.pos(), v2.pos(), u1, u2);
  } else if (uIsLine) {
    processLineArcIntr(u1.pos(), u2.pos(), v1, v2);
  } else {
    auto arc1 = arcRadiusAndCenter(v1, v2);
    auto arc2 = arcRadiusAndCenter(u1, u2);

    auto startAndSweepAngle = [](Vector2<Real> const &sp, Vector2<Real> const &center, Real bulge) {
      Real startAngle = utils::normalizeRadians(angle(center, sp));
      Real sweepAngle = Real(4) * std::atan(bulge);
      return std::make_pair(startAngle, sweepAngle);
    };

    auto bothArcsSweepPoint = [&](Vector2<Real> const &pt) {
      return pointWithinArcSweepAngle(arc1.center, v1.pos(), v2.pos(), v1.bulge(), pt) &&
             pointWithinArcSweepAngle(arc2.center, u1.pos(), u2.pos(), u1.bulge(), pt);
    };

    auto intrResult = intrCircle2Circle2(arc1.radius, arc1.center, arc2.radius, arc2.center);

    switch (intrResult.intrType) {
    case Circle2Circle2IntrType::NoIntersect:
      result.intrType = PlineSegIntrType::NoIntersect;
      break;
    case Circle2Circle2IntrType::OneIntersect:
      if (bothArcsSweepPoint(intrResult.point1)) {
        result.intrType = PlineSegIntrType::OneIntersect;
        result.point1 = intrResult.point1;
      } else {
        result.intrType = PlineSegIntrType::NoIntersect;
      }
      break;
    case Circle2Circle2IntrType::TwoIntersects: {
      const bool pt1InSweep = bothArcsSweepPoint(intrResult.point1);
      const bool pt2InSweep = bothArcsSweepPoint(intrResult.point2);
      if (pt1InSweep && pt2InSweep) {
        result.intrType = PlineSegIntrType::TwoIntersects;
        result.point1 = intrResult.point1;
        result.point2 = intrResult.point2;
      } else if (pt1InSweep) {
        result.intrType = PlineSegIntrType::OneIntersect;
        result.point1 = intrResult.point1;
      } else if (pt2InSweep) {
        result.intrType = PlineSegIntrType::OneIntersect;
        result.point1 = intrResult.point2;
      } else {
        result.intrType = PlineSegIntrType::NoIntersect;
      }
    } break;
    case Circle2Circle2IntrType::Coincident:
      // determine if arcs overlap along their sweep
      // start and sweep angles
      auto arc1StartAndSweep = startAndSweepAngle(v1.pos(), arc1.center, v1.bulge());
      auto arc2StartAndSweep = startAndSweepAngle(u1.pos(), arc2.center, u1.bulge());
      // end angles (start + sweep)
      auto arc1End = arc1StartAndSweep.first + arc1StartAndSweep.second;
      auto arc2End = arc2StartAndSweep.first + arc2StartAndSweep.second;

      if (utils::fuzzyEqual(arc1StartAndSweep.first, arc2End)) {
        // only end points touch at start of arc1
        result.intrType = PlineSegIntrType::OneIntersect;
        result.point1 = v1.pos();
      } else if (utils::fuzzyEqual(arc2StartAndSweep.first, arc1End)) {
        // only end points touch at start of arc2
        result.intrType = PlineSegIntrType::OneIntersect;
        result.point1 = u1.pos();
      } else {
        const bool arc2StartsInArc1Sweep = utils::angleIsWithinSweep(
            arc1StartAndSweep.first, arc1StartAndSweep.second, arc2StartAndSweep.first);
        const bool arc2EndsInArc1Sweep =
            utils::angleIsWithinSweep(arc1StartAndSweep.first, arc1StartAndSweep.second, arc2End);
        if (arc2StartsInArc1Sweep && arc2EndsInArc1Sweep) {
          // arc2 is fully overlapped by arc1
          result.intrType = PlineSegIntrType::ArcOverlap;
          result.point1 = u1.pos();
          result.point2 = u2.pos();
        } else if (arc2StartsInArc1Sweep) {
          // overlap from arc2 start to arc1 end
          result.intrType = PlineSegIntrType::ArcOverlap;
          result.point1 = u1.pos();
          result.point2 = v2.pos();
        } else if (arc2EndsInArc1Sweep) {
          // overlap from arc1 start to arc2 end
          result.intrType = PlineSegIntrType::ArcOverlap;
          result.point1 = v1.pos();
          result.point2 = u2.pos();
        } else {
          const bool arc1StartsInArc2Sweep = utils::angleIsWithinSweep(
              arc2StartAndSweep.first, arc2StartAndSweep.second, arc1StartAndSweep.first);
          if (arc1StartsInArc2Sweep) {
            result.intrType = PlineSegIntrType::ArcOverlap;
            result.point1 = v1.pos();
            result.point2 = v2.pos();
          } else {
            result.intrType = PlineSegIntrType::NoIntersect;
          }
        }
      }

      break;
    }
  }

  return result;
}
}

template <typename Real>
Polyline<Real> createRawOffsetPline(Polyline<Real> const &pline, Real offset,
                                    int plineCurveOrientation) {
  const bool inwardOffset = (offset > Real(0)) == (plineCurveOrientation > Real(0));
  const int connectionArcOrientation =
      inwardOffset ? -plineCurveOrientation : plineCurveOrientation;

  Polyline<Real> result;
  result.isClosed() = pline.isClosed();

  std::vector<PlineOffsetSegment<Real>> rawOffsets = createUntrimmedOffsetSegments(pline, offset);

  auto joinResultVisitor = [connectionArcOrientation](PlineOffsetSegment<Real> const &s1,
                                                      PlineOffsetSegment<Real> const &s2,
                                                      Polyline<Real> &result) {
    const bool s1IsLine = s1.v1.bulgeIsZero();
    const bool s2IsLine = s2.v1.bulgeIsZero();
    if (s1IsLine && s2IsLine) {
      detail::lineToLineJoin(s1, s2, connectionArcOrientation, result);
    } else if (s1IsLine) {
      detail::lineToArcJoin(s1, s2, connectionArcOrientation, result);
    } else if (s2IsLine) {
      detail::arcToLineJoin(s1, s2, connectionArcOrientation, result);
    } else {
      detail::arcToArcJoin(s1, s2, connectionArcOrientation, result);
    }
  };

  result.addVertex(rawOffsets[0].v1);

  for (std::size_t i = 1; i < rawOffsets.size(); ++i) {
    const auto &seg1 = rawOffsets[i - 1];
    const auto &seg2 = rawOffsets[i];
    joinResultVisitor(seg1, seg2, result);
  }

  if (pline.isClosed() && result.size() > 1) {
    // joining curves at vertex indexes (n, 0) and (0, 1)
    const auto &s1 = rawOffsets.back();
    const auto &s2 = rawOffsets[0];

    // temp polyline to capture results of joining (to avoid mutating result)
    Polyline<Real> closingPartResult;
    closingPartResult.addVertex(result.lastVertex());
    joinResultVisitor(s1, s2, closingPartResult);

    // update last vertexes
    result.lastVertex() = closingPartResult[0];
    for (std::size_t i = 1; i < closingPartResult.size(); ++i) {
      result.addVertex(closingPartResult[i]);
    }
    result.vertexes().pop_back();

    // update first vertex
    const Vector2<Real> &updatedFirstPos = closingPartResult.lastVertex().pos();
    if (result[0].bulgeIsZero()) {
      // just update position
      result[0].pos() = updatedFirstPos;
    } else {
      // update position and bulge
      const auto arc = arcRadiusAndCenter(result[0], result[1]);
      const Real a1 = angle(arc.center, updatedFirstPos);
      const Real a2 = angle(arc.center, result[1].pos());
      const Real updatedTheta = utils::deltaAngle(a1, a2);
      if ((updatedTheta < Real(0) && result[0].bulge() > Real(0)) ||
          (updatedTheta > Real(0) && result[0].bulge() < Real(0))) {
        // first vertex not valid, just update its position to be removed later
        result[0].pos() = updatedFirstPos;
      } else {
        // update position and bulge
        result[0].pos() = updatedFirstPos;
        result[0].bulge() = std::tan(updatedTheta / Real(4));
      }
    }

    // must do final singularity prune between first and second vertex after joining curves (n, 0)
    // and (0, 1)
    if (result.size() > 1) {
      if (fuzzyEqual(result[0].pos(), result[1].pos(), utils::realPrecision<Real>)) {
        result.vertexes().erase(result.vertexes().begin());
      }
    }
  } else {
    detail::addOrReplaceIfSamePos(result, rawOffsets.back().v2);
  }

  return result;
}

enum class PlineIntersectType { Simple, Tangent, Coincident };

template <typename Real> struct PlineSelfIntersect {
  // index of start vertex of first segment
  std::size_t sIndex1;
  // index of start vertex of second segment
  std::size_t sIndex2;
  // intersect position
  Vector2<Real> pos;
  // type of intersect
  PlineIntersectType intrType;
  PlineSelfIntersect() = default;
  PlineSelfIntersect(std::size_t si1, std::size_t si2, Vector2<Real> p, PlineIntersectType iType)
      : sIndex1(si1), sIndex2(si2), pos(p), intrType(iType) {}
};

// Finds all local self intersects of the polyline, local self intersects are defined as between two
// polyline segments that share a vertex.
// NOTES:
// - Singularities (repeating vertexes) are returned as coincident intersects
template <typename Real>
void localSelfIntersects(Polyline<Real> const &pline,
                         std::vector<PlineSelfIntersect<Real>> &output) {
  if (pline.size() < 2) {
    return;
  }

  if (pline.size() == 2) {
    if (pline.isClosed()) {
      // check if overlaps on itself from vertex 1 to vertex 2
      if (utils::fuzzyEqual(pline[0].bulge(), -pline[1].bulge())) {
        output.emplace_back(0, 1, pline[1].pos(), PlineIntersectType::Coincident);
        output.emplace_back(1, 0, pline[0].pos(), PlineIntersectType::Coincident);
      }
    }
    return;
  }

  auto testAndAddIntersect = [&](std::size_t i, std::size_t j, std::size_t k) {
    const PlineVertex<Real> &v1 = pline[i];
    const PlineVertex<Real> &v2 = pline[j];
    const PlineVertex<Real> &v3 = pline[k];
    // testing intersection between v1->v2 and v2->v3 segments

    if (fuzzyEqual(v1.pos(), v2.pos(), utils::realPrecision<Real>)) {
      // singularity
      output.emplace_back(i, j, v1.pos(), PlineIntersectType::Coincident);
    } else {
      using namespace detail;
      IntrPlineSegsResult<Real> intrResult = intrPlineSegs(v1, v2, v2, v3);
      switch (intrResult.intrType) {
      case PlineSegIntrType::NoIntersect:
        break;
      case PlineSegIntrType::TangentIntersect:
        if (!fuzzyEqual(intrResult.point1, v2.pos(), utils::realPrecision<Real>)) {
          output.emplace_back(i, j, intrResult.point1, PlineIntersectType::Tangent);
        }
        break;
      case PlineSegIntrType::OneIntersect:
        if (!fuzzyEqual(intrResult.point1, v2.pos(), utils::realPrecision<Real>)) {
          output.emplace_back(i, j, intrResult.point1, PlineIntersectType::Simple);
        }
        break;
      case PlineSegIntrType::TwoIntersects:
        if (!fuzzyEqual(intrResult.point1, v2.pos(), utils::realPrecision<Real>)) {
          output.emplace_back(i, j, intrResult.point1, PlineIntersectType::Simple);
        }
        if (!fuzzyEqual(intrResult.point2, v2.pos(), utils::realPrecision<Real>)) {
          output.emplace_back(i, j, intrResult.point2, PlineIntersectType::Simple);
        }
        break;
      case PlineSegIntrType::SegmentOverlap:
      case PlineSegIntrType::ArcOverlap:
        output.emplace_back(i, j, intrResult.point1, PlineIntersectType::Coincident);
        break;
      }
    }

  };

  for (std::size_t i = 2; i < pline.size(); ++i) {
    testAndAddIntersect(i - 2, i - 1, i);
  }

  if (pline.isClosed()) {
    // we tested for intersect between segments at indexes 0->1, 1->2 and everything up to and
    // including (count-3)->(count-2), (count-2)->(count-1), polyline is closed so now test
    // [(count-2)->(count-1), (count-1)->0] and [(count-1)->0, 0->1]
    testAndAddIntersect(pline.size() - 2, pline.size() - 1, 0);
    testAndAddIntersect(pline.size() - 1, 0, 1);
  }
}

// Finds all global self intersects of the polyline, global self intersects are defined as all
// intersects between polyline segments that DO NOT share a vertex (use the localSelfIntersects
// function to find those). A spatial index is used to minimize the intersect comparisons required,
// the spatial index should hold bounding boxes for all of the polyline's segments.
// NOTES:
// - We never include intersects at a segment's start point, the matching intersect from the
// previous segment's end point is included (no sense in including both)
template <typename Real, std::size_t N>
void globalSelfIntersects(Polyline<Real> const &pline,
                          std::vector<PlineSelfIntersect<Real>> &output,
                          StaticSpatialIndex<Real, N> const &spatialIndex) {
  if (pline.size() < 3) {
    return;
  }
  using namespace detail;

  std::unordered_set<std::pair<std::size_t, std::size_t>, IndexPairHash> visitedSegmentPairs;

  auto visitor = [&](std::size_t i, std::size_t j) {
    const PlineVertex<Real> &v1 = pline[i];
    const PlineVertex<Real> &v2 = pline[j];
    AABB<Real> envelope = createFastApproxBoundingBox(v1, v2);
    envelope.expand(utils::realThreshold<Real>);
    auto count = 0;
    auto indexVisitor = [&](std::size_t hitIndexStart) {
      std::size_t hitIndexEnd = utils::nextWrappingIndex(hitIndexStart, pline);
      // skip/filter already visited intersects
      // skip local segments
      if (i == hitIndexStart || i == hitIndexEnd || j == hitIndexStart || j == hitIndexEnd) {
        return true;
      }
      count++;
      // skip reversed segment order (would end up comparing the same segments)
      if (visitedSegmentPairs.find(std::make_pair(hitIndexStart, i)) != visitedSegmentPairs.end()) {
        return true;
      }

      // add the segment pair we're visiting now
      visitedSegmentPairs.emplace(i, hitIndexStart);

      const PlineVertex<Real> &u1 = pline[hitIndexStart];
      const PlineVertex<Real> &u2 = pline[hitIndexEnd];

      auto intrAtStartPt = [&](Vector2<Real> const &intr) {
        return fuzzyEqual(v1.pos(), intr) || fuzzyEqual(u1.pos(), intr);
      };

      IntrPlineSegsResult<Real> intrResult = intrPlineSegs(v1, v2, u1, u2);
      switch (intrResult.intrType) {
      case PlineSegIntrType::NoIntersect:
        break;
      case PlineSegIntrType::TangentIntersect:
        if (!intrAtStartPt(intrResult.point1)) {
          output.emplace_back(i, hitIndexStart, intrResult.point1, PlineIntersectType::Tangent);
        }
        break;
      case PlineSegIntrType::OneIntersect:
        if (!intrAtStartPt(intrResult.point1)) {
          output.emplace_back(i, hitIndexStart, intrResult.point1, PlineIntersectType::Simple);
        }
        break;
      case PlineSegIntrType::TwoIntersects:
        if (!intrAtStartPt(intrResult.point1)) {
          output.emplace_back(i, hitIndexStart, intrResult.point1, PlineIntersectType::Simple);
        }
        if (!intrAtStartPt(intrResult.point2)) {
          output.emplace_back(i, hitIndexStart, intrResult.point2, PlineIntersectType::Simple);
        }
        break;
      case PlineSegIntrType::SegmentOverlap:
      case PlineSegIntrType::ArcOverlap:
        if (!intrAtStartPt(intrResult.point1)) {
          output.emplace_back(i, hitIndexStart, intrResult.point1, PlineIntersectType::Coincident);
        }
        if (!intrAtStartPt(intrResult.point2)) {
          output.emplace_back(i, hitIndexStart, intrResult.point2, PlineIntersectType::Coincident);
        }
        break;
      }

      // visit the entire query
      return true;
    };

    spatialIndex.visitQuery(envelope.xMin, envelope.yMin, envelope.xMax, envelope.yMax,
                            indexVisitor);

    // visit all pline indexes
    return true;
  };

  iterateSegIndices(pline, visitor);
}

// Finds all self intersects of the polyline (equivalent to calling localSelfIntersects and
// globalSelfIntersects).
template <typename Real, std::size_t N>
void allSelfIntersects(Polyline<Real> const &pline, std::vector<PlineSelfIntersect<Real>> &output,
                       StaticSpatialIndex<Real, N> const &spatialIndex) {
  localSelfIntersects(pline, output);
  globalSelfIntersects(pline, output, spatialIndex);
}

template <typename Real, std::size_t N>
bool pointValidForOffset(Polyline<Real> const &pline, Real offset,
                         StaticSpatialIndex<Real, N> const &spatialIndex,
                         Vector2<Real> const &point, Real offsetTol = Real(1e-3)) {
  const Real absOffset = std::abs(offset) - offsetTol;
  const Real minDist = absOffset * absOffset;

  bool pointValid = true;

  auto visitor = [&](std::size_t i) {
    std::size_t j = utils::nextWrappingIndex(i, pline.vertexes());
    auto closestPoint = closestPointOnSeg(pline[i], pline[j], point);
    Real dist = distSquared(closestPoint, point);
    pointValid = dist > minDist;
    return pointValid;
  };

  spatialIndex.visitQuery(point.x() - absOffset, point.y() - absOffset, point.x() + absOffset,
                          point.y() + absOffset, visitor);
  return pointValid;
}

template <typename Real> struct OpenPolylineSlice {
  std::size_t intrStartIndex;
  Polyline<Real> pline;
  OpenPolylineSlice() = default;
  OpenPolylineSlice(std::size_t sIndex, Polyline<Real> slice)
      : intrStartIndex(sIndex), pline(std::move(slice)) {}
};

template <typename Real>
std::vector<OpenPolylineSlice<Real>> sliceAtIntersects(Polyline<Real> const &originalPline,
                                                       Polyline<Real> const &rawOffsetPline,
                                                       Real offset) {
  std::vector<OpenPolylineSlice<Real>> result;
  if (rawOffsetPline.size() < 2) {
    return result;
  }

  StaticSpatialIndex<Real> origPlineSpatialIndex = createApproxSpatialIndex(originalPline);
  StaticSpatialIndex<Real> rawOffsetPlineSpatialIndex = createApproxSpatialIndex(rawOffsetPline);

  std::vector<PlineSelfIntersect<Real>> selfIntersects;
  allSelfIntersects(rawOffsetPline, selfIntersects, rawOffsetPlineSpatialIndex);

  if (selfIntersects.size() == 0) {
    // no self intersects, test that all vertexes are valid distance from original polyline
    bool isValid = std::all_of(rawOffsetPline.vertexes().begin(), rawOffsetPline.vertexes().end(),
                               [&](PlineVertex<Real> const &v) {
                                 return pointValidForOffset(originalPline, offset,
                                                            origPlineSpatialIndex, v.pos());
                               });

    if (!isValid) {
      return result;
    }

    // copy and convert raw offset into open polyline
    result.emplace_back(std::numeric_limits<std::size_t>::max(), rawOffsetPline);
    result.back().pline.isClosed() = false;
    result.back().pline.addVertex(rawOffsetPline[0]);
    result.back().pline.lastVertex().bulge() = Real(0);
    return result;
  }

  std::unordered_map<std::size_t, std::vector<PlineSelfIntersect<Real>>> selfIntersectLookup;

  for (PlineSelfIntersect<Real> const &si : selfIntersects) {
    selfIntersectLookup[si.sIndex1].push_back(si);
    selfIntersectLookup[si.sIndex2].push_back(si);
  }

  // sort intersects by distance from start vertex
  for (auto &kvp : selfIntersectLookup) {
    Vector2<Real> startPos = rawOffsetPline[kvp.first].pos();
    auto cmp = [&](PlineSelfIntersect<Real> const &si1, PlineSelfIntersect<Real> const &si2) {
      return distSquared(si1.pos, startPos) < distSquared(si2.pos, startPos);
    };
    std::sort(kvp.second.begin(), kvp.second.end(), cmp);
  }

  auto intersectsOrigPline = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
    AABB<Real> approxBB = createFastApproxBoundingBox(v1, v2);
    bool intersects = false;
    auto visitor = [&](std::size_t i) {
      using namespace detail;
      std::size_t j = utils::nextWrappingIndex(i, originalPline);
      IntrPlineSegsResult<Real> intrResult =
          intrPlineSegs(v1, v2, originalPline[i], originalPline[j]);
      intersects = intrResult.intrType != PlineSegIntrType::NoIntersect;
      return !intersects;
    };

    origPlineSpatialIndex.visitQuery(approxBB.xMin, approxBB.yMin, approxBB.xMax, approxBB.yMax,
                                     visitor);

    return intersects;
  };

  for (auto const &kvp : selfIntersectLookup) {
    // start index for the slice we're about to build
    std::size_t sIndex = kvp.first;
    // self intersect list for this start index
    std::vector<PlineSelfIntersect<Real>> const &siList = kvp.second;

    const auto &startVertex = rawOffsetPline[sIndex];
    std::size_t nextIndex = utils::nextWrappingIndex(sIndex, rawOffsetPline);
    const auto &endVertex = rawOffsetPline[nextIndex];

    if (siList.size() != 1) {
      // build all the segments between the N intersects in siList (N > 1)
      SplitResult<Real> firstSplit = splitAtPoint(startVertex, endVertex, siList[0].pos);
      auto prevVertex = firstSplit.splitVertex;
      for (std::size_t i = 1; i < siList.size(); ++i) {
        SplitResult<Real> split = splitAtPoint(prevVertex, endVertex, siList[i].pos);
        // update prevVertex for next loop iteration
        prevVertex = split.splitVertex;
        // skip if they're ontop of each other
        if (fuzzyEqual(split.updatedStart.pos(), split.splitVertex.pos(),
                       utils::realPrecision<Real>)) {
          continue;
        }

        // test start point
        if (!pointValidForOffset(originalPline, offset, origPlineSpatialIndex,
                                 split.updatedStart.pos())) {
          continue;
        }

        // test end point
        if (!pointValidForOffset(originalPline, offset, origPlineSpatialIndex,
                                 split.splitVertex.pos())) {
          continue;
        }

        // test mid point
        auto midpoint = detail::segMidpoint(split.updatedStart, split.splitVertex);
        if (!pointValidForOffset(originalPline, offset, origPlineSpatialIndex, midpoint)) {
          continue;
        }

        // test intersection with original polyline
        if (intersectsOrigPline(split.updatedStart, split.splitVertex)) {
          continue;
        }

        result.emplace_back();
        result.back().intrStartIndex = sIndex;
        result.back().pline.addVertex(split.updatedStart);
        result.back().pline.addVertex(split.splitVertex);
      }
    }

    // build the segment between the last intersect in siList and the next intersect found

    // check that the first point is valid
    if (!pointValidForOffset(originalPline, offset, origPlineSpatialIndex, siList.back().pos)) {
      continue;
    }

    SplitResult<Real> split = splitAtPoint(startVertex, endVertex, siList.back().pos);
    Polyline<Real> currPline;
    currPline.addVertex(split.splitVertex);

    std::size_t index = nextIndex;
    bool isValidPline = true;
    while (true) {
      // check that vertex point is valid
      if (!pointValidForOffset(originalPline, offset, origPlineSpatialIndex,
                               rawOffsetPline[index].pos())) {
        isValidPline = false;
        break;
      }

      // check that the segment does not intersect original polyline
      if (intersectsOrigPline(currPline.lastVertex(), rawOffsetPline[index])) {
        isValidPline = false;
        break;
      }

      // add vertex
      detail::addOrReplaceIfSamePos(currPline, rawOffsetPline[index]);

      // check if segment that starts at vertex we just added has an intersect
      auto nextIntr = selfIntersectLookup.find(index);
      if (nextIntr != selfIntersectLookup.end()) {
        // there is an intersect, slice is done, check if final segment is valid

        // check intersect pos is valid (which will also be end vertex position)
        Vector2<Real> const &intersectPos = nextIntr->second[0].pos;
        if (!pointValidForOffset(originalPline, offset, origPlineSpatialIndex, intersectPos)) {
          isValidPline = false;
          break;
        }

        std::size_t nextIndex = utils::nextWrappingIndex(index, rawOffsetPline);
        SplitResult<Real> split =
            splitAtPoint(currPline.lastVertex(), rawOffsetPline[nextIndex], intersectPos);

        // check mid point is valid
        Vector2<Real> mp =
            detail::segMidpoint(split.updatedStart, PlineVertex<Real>(intersectPos, Real(0)));
        if (!pointValidForOffset(originalPline, offset, origPlineSpatialIndex, mp)) {
          isValidPline = false;
          break;
        }

        // trim last added vertex and add final intersect position
        currPline.lastVertex() = split.updatedStart;
        currPline.addVertex(PlineVertex<Real>(intersectPos, Real(0)));

        break;
      }
      // else there is not an intersect, increment index and continue
      index = utils::nextWrappingIndex(index, rawOffsetPline);
    }

    if (isValidPline && currPline.size() > 1) {
      result.emplace_back(sIndex, std::move(currPline));
    }
  }

  return result;
}

template <typename Real>
std::vector<Polyline<Real>> stitchSlicesTogether(std::vector<OpenPolylineSlice<Real>> const &slices,
                                                 std::size_t origMaxIndex,
                                                 Real joinThreshold = utils::realPrecision<Real>) {
  std::vector<Polyline<Real>> result;
  if (slices.size() == 0) {
    return result;
  }

  if (slices.size() == 1) {
    if (fuzzyEqual(slices[0].pline.lastVertex().pos(), slices[0].pline[0].pos(), joinThreshold)) {
      result.emplace_back(slices[0].pline);
      result.back().isClosed() = true;
      result.back().vertexes().pop_back();
    }

    return result;
  }

  // load spatial index with all start points
  StaticSpatialIndex<Real> spatialIndex(slices.size());

  for (const auto &slice : slices) {
    auto const &point = slice.pline[0].pos();
    spatialIndex.add(point.x() - joinThreshold, point.y() - joinThreshold,
                     point.x() + joinThreshold, point.y() + joinThreshold);
  }

  spatialIndex.finish();

  std::unordered_set<std::size_t> visitedIndexes;
  std::vector<std::size_t> queryResults;
  for (std::size_t i = 0; i < slices.size(); ++i) {
    auto insertResult = visitedIndexes.insert(i);
    if (!insertResult.second) {
      // already visited
      continue;
    }

    Polyline<Real> closedPline;
    closedPline.isClosed() = true;
    std::size_t currIndex = i;
    auto const &initialStartPoint = slices[i].pline[0].pos();

    while (true) {
      const std::size_t currLoopStartIndex = slices[currIndex].intrStartIndex;
      auto const &currSlice = slices[currIndex].pline;
      auto const &currEndPoint = slices[currIndex].pline.lastVertex().pos();
      closedPline.vertexes().insert(closedPline.vertexes().end(), currSlice.vertexes().begin(),
                                    currSlice.vertexes().end());
      queryResults.clear();
      spatialIndex.query(currEndPoint.x() - joinThreshold, currEndPoint.y() - joinThreshold,
                         currEndPoint.x() + joinThreshold, currEndPoint.y() + joinThreshold,
                         queryResults);

      queryResults.erase(std::remove_if(queryResults.begin(), queryResults.end(),
                                        [&](std::size_t index) {
                                          return visitedIndexes.find(index) != visitedIndexes.end();
                                        }),
                         queryResults.end());

      auto indexDistAndEqualInitial = [&](std::size_t index) {
        auto const &slice = slices[index];
        std::size_t indexDist;
        if (currLoopStartIndex <= slice.intrStartIndex) {
          indexDist = slice.intrStartIndex - currLoopStartIndex;
        } else {
          // forward wrapping distance (distance to end + distance to index)
          indexDist = origMaxIndex - currLoopStartIndex + slice.intrStartIndex;
        }

        bool equalToInitial = fuzzyEqual(slice.pline.lastVertex().pos(), initialStartPoint,
                                         utils::realPrecision<Real>);

        return std::make_pair(indexDist, equalToInitial);
      };

      std::sort(queryResults.begin(), queryResults.end(),
                [&](std::size_t index1, std::size_t index2) {
                  auto distAndEqualInitial1 = indexDistAndEqualInitial(index1);
                  auto distAndEqualInitial2 = indexDistAndEqualInitial(index2);
                  if (distAndEqualInitial1.first == distAndEqualInitial2.first) {
                    if (!distAndEqualInitial1.second) {
                      // return true (index1 < index2, to prioritize index1, we want the longest
                      // closed loop possible)
                      return true;
                    }

                    return false;
                  }

                  return distAndEqualInitial1.first < distAndEqualInitial2.first;
                });

      if (queryResults.size() == 0) {
        // we're done
        if (closedPline.size() > 1) {
          if (fuzzyEqual(closedPline[0].pos(), closedPline.lastVertex().pos(),
                         utils::realPrecision<Real>)) {
            closedPline.vertexes().pop_back();
          }
          result.emplace_back(std::move(closedPline));
        }
        break;
      }

      // else continue stitching
      visitedIndexes.insert(queryResults[0]);
      closedPline.vertexes().pop_back();
      currIndex = queryResults[0];
    }
  }

  return result;
}

template <typename Real>
std::vector<Polyline<Real>> paralleOffset(Polyline<Real> const &pline, Real offset,
                                          int orientation = 0) {
  if (orientation == 0) {
    orientation = area(pline) > Real(0) ? 1 : -1;
  }
  auto rawOffset = createRawOffsetPline(pline, offset, orientation);
  auto slices = sliceAtIntersects(pline, rawOffset, offset);
  return stitchSlicesTogether(slices, rawOffset.size() - 1);
}
}

#endif // POLYLINE_H
