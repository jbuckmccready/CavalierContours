#include "cavaliercontours.h"
#include "gmock/gmock.h"
#include "testhelpers.hpp"
#include "gtest/gtest.h"
#include <sstream>
#include <string>
#include <vector>

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

std::ostream &operator<<(std::ostream &os, cavc_vertex const &v) {
  os << '[' << v.x << "," << v.y << "," << v.bulge << ']';
  return os;
}

std::ostream &operator<<(std::ostream &os, cavc_point const &p) {
  os << '[' << p.x << "," << p.y << ']';
  return os;
}

struct cavc_plineTestCase {
  std::string name;
  cavc_pline *pline;
  cavc_real signedArea;
  cavc_real pathLength;
  cavc_real minX;
  cavc_real minY;
  cavc_real maxX;
  cavc_real maxY;

  std::vector<cavc_point> windingNumberTestPts;
  std::vector<int> windingNumberResults;

  void addWindingNumberTestPt(cavc_point testPt, int result) {
    windingNumberTestPts.push_back(testPt);
    windingNumberResults.push_back(result);
  }

  std::vector<cavc_point> closestPointTestPts;
  std::vector<uint32_t> closestPointIndexResults;
  std::vector<cavc_point> closestPointResults;
  std::vector<cavc_real> closestPointDistanceResults;

  void addClosestPointTestPt(cavc_point testPt, cavc_point closestPointResult,
                             cavc_real distanceResult,
                             uint32_t indexResult = std::numeric_limits<uint32_t>::max()) {
    closestPointTestPts.push_back(testPt);
    closestPointResults.push_back(closestPointResult);
    closestPointDistanceResults.push_back(distanceResult);
    closestPointIndexResults.push_back(indexResult);
  }

  std::vector<cavc_vertex> plineVertexes;

  bool isClosed() const { return cavc_pline_is_closed(pline); }
};

std::ostream &operator<<(std::ostream &os, cavc_plineTestCase const &c) {
  os << c.name;
  return os;
}

void addCircleCases(std::vector<cavc_plineTestCase> &cases, cavc_real circleRadius,
                    cavc_point circleCenter, bool reverse) {
  cavc_real expectedCircleArea = CAVC_PI * circleRadius * circleRadius;
  cavc_real expectedCircleLength = 2 * CAVC_PI * circleRadius;

  auto createName = [&](auto n) {
    std::stringstream ss;
    ss << n;
    if (reverse) {
      ss << "_rev";
    }
    ss << " (radius: " << circleRadius << ", center: " << circleCenter << ")";
    return ss.str();
  };

  auto addCase = [&](std::string name, std::vector<cavc_vertex> &vertexes, int direction) {
    cavc_plineTestCase test_case;
    test_case.name = createName(name);
    test_case.pline = cavc_pline_new(&vertexes[0], 2, 1);
    test_case.signedArea = direction * expectedCircleArea;
    test_case.pathLength = expectedCircleLength;
    test_case.minX = circleCenter.x - circleRadius;
    test_case.minY = circleCenter.y - circleRadius;
    test_case.maxX = circleCenter.x + circleRadius;
    test_case.maxY = circleCenter.y + circleRadius;

    // axis aligned and center cases
    test_case.addWindingNumberTestPt({test_case.minX - 0.01, circleCenter.y}, 0);
    test_case.addWindingNumberTestPt({test_case.maxX + 0.01, circleCenter.y}, 0);
    test_case.addWindingNumberTestPt({circleCenter.x, test_case.minY - 0.01}, 0);
    test_case.addWindingNumberTestPt({circleCenter.x, test_case.maxY + 0.01}, 0);
    test_case.addWindingNumberTestPt(circleCenter, direction);
    test_case.addWindingNumberTestPt({test_case.minX + 0.01, circleCenter.y}, direction);
    test_case.addWindingNumberTestPt({test_case.maxX - 0.01, circleCenter.y}, direction);
    test_case.addWindingNumberTestPt({circleCenter.x, test_case.minY + 0.01}, direction);
    test_case.addWindingNumberTestPt({circleCenter.x, test_case.maxY - 0.01}, direction);

    for (std::size_t i = 0; i < vertexes.size(); ++i) {
      auto const &v = vertexes[i];
      test_case.addClosestPointTestPt({v.x, v.y}, {v.x, v.y}, 0.0, static_cast<uint32_t>(i));
    }

    test_case.addClosestPointTestPt({circleCenter.x - 0.1, circleCenter.y},
                                    {circleCenter.x - circleRadius, circleCenter.y},
                                    circleRadius - 0.1);
    test_case.addClosestPointTestPt({circleCenter.x + 0.1, circleCenter.y},
                                    {circleCenter.x + circleRadius, circleCenter.y},
                                    circleRadius - 0.1);
    test_case.addClosestPointTestPt({circleCenter.x, circleCenter.y - 0.1},
                                    {circleCenter.x, circleCenter.y - circleRadius},
                                    circleRadius - 0.1);
    test_case.addClosestPointTestPt({circleCenter.x, circleCenter.y + 0.1},
                                    {circleCenter.x, circleCenter.y + circleRadius},
                                    circleRadius - 0.1);

    // points at 45 deg inside and outside
    std::vector<cavc_point> insideAt45Deg;
    insideAt45Deg.reserve(4);
    std::vector<cavc_point> outsideAt45Deg;
    outsideAt45Deg.reserve(4);
    cavc_real insideDist = 0.33 * circleRadius;
    cavc_real outsideDist = 1.5 * circleRadius;
    std::vector<cavc_point> onCirclAt45deg;
    onCirclAt45deg.reserve(4);
    for (std::size_t i = 0; i < 4; ++i) {
      cavc_real xCos = std::cos(CAVC_PI / 4 + i * CAVC_PI / 2);
      cavc_real ySin = std::sin(CAVC_PI / 4 + i * CAVC_PI / 2);

      cavc_real x = circleCenter.x + insideDist * xCos;
      cavc_real y = circleCenter.y + insideDist * ySin;
      insideAt45Deg.push_back({x, y});

      x = circleCenter.x + outsideDist * xCos;
      y = circleCenter.y + outsideDist * ySin;
      outsideAt45Deg.push_back({x, y});

      x = circleCenter.x + circleRadius * xCos;
      y = circleCenter.y + circleRadius * ySin;
      onCirclAt45deg.push_back({x, y});
    }

    for (std::size_t i = 0; i < insideAt45Deg.size(); ++i) {
      test_case.addWindingNumberTestPt(insideAt45Deg[i], direction);
      test_case.addWindingNumberTestPt(outsideAt45Deg[i], 0);
      test_case.addClosestPointTestPt(insideAt45Deg[i], onCirclAt45deg[i],
                                      circleRadius - insideDist);
      test_case.addClosestPointTestPt(outsideAt45Deg[i], onCirclAt45deg[i],
                                      outsideDist - circleRadius);
    }

    cases.push_back(std::move(test_case));
  };

  std::vector<cavc_vertex> vertexes = {{circleCenter.x - circleRadius, circleCenter.y, 1},
                                       {circleCenter.x + circleRadius, circleCenter.y, 1}};
  if (reverse) {
    std::reverse(vertexes.begin(), vertexes.end());
  }

  addCase("ccw_circle_x_aligned", vertexes, 1);

  for (auto &v : vertexes) {
    v = {v.x, v.y, -v.bulge};
  }
  addCase("cw_circle_x_aligned", vertexes, -1);

  vertexes = {{circleCenter.x, circleCenter.y - circleRadius, 1},
              {circleCenter.x, circleCenter.y + circleRadius, 1}};

  if (reverse) {
    std::reverse(vertexes.begin(), vertexes.end());
  }

  addCase("ccw_circle_y_aligned", vertexes, 1);
  for (auto &v : vertexes) {
    v = {v.x, v.y, -v.bulge};
  }
  addCase("cw_circle_y_aligned", vertexes, -1);
}

std::vector<cavc_plineTestCase> createCircleCases() {
  std::vector<cavc_plineTestCase> result;
  addCircleCases(result, 5.0, {1, 1}, false);
  addCircleCases(result, 5.0, {-1, 1}, false);
  addCircleCases(result, 5.0, {-1, -1}, false);
  addCircleCases(result, 5.0, {1, -1}, false);

  addCircleCases(result, 5.0, {1, 1}, true);
  addCircleCases(result, 5.0, {-1, 1}, true);
  addCircleCases(result, 5.0, {-1, -1}, true);
  addCircleCases(result, 5.0, {1, -1}, true);
  return result;
}

class cavc_plineTests : public ::testing::Test {
protected:
  void SetUp() override;

  void TearDown() override;

  std::vector<cavc_vertex> pline1Vertexes;
  cavc_pline *pline1;
  cavc_pline *pline2;
  std::size_t origPline1Size;
  uint32_t initialPline1Size() { return static_cast<uint32_t>(origPline1Size); }
};

void cavc_plineTests::SetUp() {
  pline1Vertexes = {{1, 2, 0.1}, {33, 3, 0.2}, {34, 35, 0.3}, {2, 36, 0.4}};
  origPline1Size = pline1Vertexes.size();
  pline1 = cavc_pline_new(&pline1Vertexes[0], initialPline1Size(), 0);
  pline2 = cavc_pline_new(nullptr, 0, 1);
}

void cavc_plineTests::TearDown() {
  cavc_pline_delete(pline1);
  cavc_pline_delete(pline2);
}

class cavc_plineFunctionTests : public testing::TestWithParam<cavc_plineTestCase> {
protected:
  void SetUp() override;
  void TearDown() override;

  static void SetUpTestSuite() {}
  static void TearDownTestSuite() {
    for (auto &c : circleCases) {
      cavc_pline_delete(c.pline);
    }
  }

public:
  static std::vector<cavc_plineTestCase> circleCases;
};
void cavc_plineFunctionTests::SetUp() {}
void cavc_plineFunctionTests::TearDown() {}

std::vector<cavc_plineTestCase> cavc_plineFunctionTests::circleCases = createCircleCases();

TEST_F(cavc_plineTests, cavc_pline_new) {

  // test capacity
  EXPECT_EQ(cavc_pline_capacity(pline1), initialPline1Size());

  // test is_closed
  EXPECT_EQ(cavc_pline_is_closed(pline1), 0);

  // test vertex_count
  auto pline1Count = cavc_pline_vertex_count(pline1);
  ASSERT_EQ(pline1Count, initialPline1Size());

  // test vertex_data
  std::vector<cavc_vertex> read_vertexes(pline1Count);
  cavc_pline_vertex_data(pline1, &read_vertexes[0]);
  ASSERT_THAT(pline1Vertexes, testing::Pointwise(VertexEqual(), read_vertexes));

  // test on empty pline
  // test capacity
  EXPECT_EQ(cavc_pline_capacity(pline2), 0);

  // test is_closed
  EXPECT_EQ(cavc_pline_is_closed(pline2), 1);

  // test vertex_count
  auto pline2Count = cavc_pline_vertex_count(pline2);
  ASSERT_EQ(pline2Count, 0);

  // test vertex_data
  cavc_pline_vertex_data(pline2, &read_vertexes[0]);
  // nothing should have been written to the buffer
  ASSERT_THAT(read_vertexes, testing::Pointwise(VertexEqual(), pline1Vertexes));
}

TEST_F(cavc_plineTests, cavc_pline_set_capacity) {
  // setting capacity less than current does nothing
  cavc_pline_set_capacity(pline1, 1);
  ASSERT_EQ(cavc_pline_capacity(pline1), 4);

  cavc_pline_set_capacity(pline1, 11);
  ASSERT_EQ(cavc_pline_capacity(pline1), 11);
}

TEST_F(cavc_plineTests, cavc_pline_set_vertex_data) {
  cavc_pline_set_vertex_data(pline2, &pline1Vertexes[0], initialPline1Size());
  ASSERT_EQ(cavc_pline_vertex_count(pline2), initialPline1Size());

  std::vector<cavc_vertex> readVertexes(initialPline1Size());
  cavc_pline_vertex_data(pline2, &readVertexes[0]);
  ASSERT_THAT(readVertexes, testing::Pointwise(VertexEqual(), pline1Vertexes));
}

TEST_F(cavc_plineTests, cavc_pline_add_vertex) {
  cavc_vertex v{555, 666, 0.777};
  cavc_pline_add_vertex(pline1, v);
  ASSERT_EQ(cavc_pline_vertex_count(pline1), initialPline1Size() + 1);

  std::vector<cavc_vertex> readVertexes(initialPline1Size() + 1);
  cavc_pline_vertex_data(pline1, &readVertexes[0]);
  pline1Vertexes.push_back(v);
  ASSERT_THAT(readVertexes, testing::Pointwise(VertexEqual(), pline1Vertexes));

  cavc_pline_add_vertex(pline2, v);
  ASSERT_EQ(cavc_pline_vertex_count(pline2), 1);

  readVertexes.resize(1);
  cavc_pline_vertex_data(pline2, &readVertexes[0]);
  ASSERT_THAT(readVertexes, testing::Pointwise(VertexEqual(), {v}));
}

TEST_F(cavc_plineTests, cavc_pline_remove_range) {
  // remove first vertex
  cavc_pline_remove_range(pline1, 0, 1);
  ASSERT_EQ(cavc_pline_vertex_count(pline1), initialPline1Size() - 1);

  std::vector<cavc_vertex> readVertexes(initialPline1Size() - 1);
  cavc_pline_vertex_data(pline1, &readVertexes[0]);
  pline1Vertexes.erase(pline1Vertexes.begin());
  ASSERT_THAT(readVertexes, testing::Pointwise(VertexEqual(), pline1Vertexes));

  // remove 2nd and 3rd vertex
  cavc_pline_remove_range(pline1, 1, 2);
  ASSERT_EQ(cavc_pline_vertex_count(pline1), initialPline1Size() - 3);
  readVertexes.resize(1);
  cavc_pline_vertex_data(pline1, &readVertexes[0]);
  pline1Vertexes.erase(pline1Vertexes.begin() + 1, pline1Vertexes.begin() + 3);
  ASSERT_THAT(readVertexes, testing::Pointwise(VertexEqual(), pline1Vertexes));

  // remove last vertex
  cavc_pline_remove_range(pline1, 0, 1);
  ASSERT_EQ(cavc_pline_vertex_count(pline1), initialPline1Size() - 4);

  readVertexes.resize(10);
  std::fill(readVertexes.begin(), readVertexes.end(), cavc_vertex{-1, -2, -3});
  auto copy = readVertexes;
  cavc_pline_vertex_data(pline1, &readVertexes[0]);
  // nothing should have been written to the buffer
  ASSERT_THAT(readVertexes, testing::Pointwise(VertexEqual(), copy));
}

TEST_F(cavc_plineTests, cavc_pline_clear) {
  cavc_pline_clear(pline1);
  ASSERT_EQ(cavc_pline_vertex_count(pline1), 0);

  cavc_pline_clear(pline2);
  ASSERT_EQ(cavc_pline_vertex_count(pline2), 0);
}

INSTANTIATE_TEST_SUITE_P(cavc_pline_circles, cavc_plineFunctionTests,
                         testing::ValuesIn(cavc_plineFunctionTests::circleCases));

TEST_P(cavc_plineFunctionTests, cavc_get_path_length) {
  cavc_plineTestCase const &test_case = GetParam();
  ASSERT_EQ(cavc_get_path_length(test_case.pline), test_case.pathLength);
}

TEST_P(cavc_plineFunctionTests, cavc_get_area) {
  cavc_plineTestCase const &test_case = GetParam();
  ASSERT_EQ(cavc_get_area(test_case.pline), test_case.signedArea);
}

TEST_P(cavc_plineFunctionTests, cavc_get_winding_number) {
  cavc_plineTestCase const &test_case = GetParam();
  std::vector<int> windingNumberResults;
  windingNumberResults.reserve(test_case.windingNumberTestPts.size());
  for (auto const &pt : test_case.windingNumberTestPts) {
    windingNumberResults.push_back(cavc_get_winding_number(test_case.pline, pt));
  }

  ASSERT_THAT(windingNumberResults,
              testing::Pointwise(testing::Eq(), test_case.windingNumberResults));
}

TEST_P(cavc_plineFunctionTests, cavc_get_extents) {
  cavc_plineTestCase const &test_case = GetParam();
  cavc_real minX;
  cavc_real minY;
  cavc_real maxX;
  cavc_real maxY;
  cavc_get_extents(test_case.pline, &minX, &minY, &maxX, &maxY);
  EXPECT_EQ(minX, test_case.minX);
  EXPECT_EQ(minY, test_case.minY);
  EXPECT_EQ(maxX, test_case.maxX);
  EXPECT_EQ(maxY, test_case.maxY);
}

TEST_P(cavc_plineFunctionTests, cavc_get_closest_point) {
  cavc_plineTestCase const &test_case = GetParam();
  std::size_t ptCount = test_case.closestPointTestPts.size();
  std::vector<uint32_t> indexResults;
  indexResults.reserve(ptCount);
  std::vector<cavc_point> pointResults;
  pointResults.reserve(ptCount);
  std::vector<cavc_real> distanceResults;
  distanceResults.reserve(ptCount);

  for (std::size_t i = 0; i < test_case.closestPointTestPts.size(); ++i) {
    auto const &pt = test_case.closestPointTestPts[i];
    uint32_t index;
    cavc_point p;
    cavc_real dist;
    cavc_get_closest_point(test_case.pline, pt, &index, &p, &dist);
    // set to max index so they compare equal (effectively skipping start index check)
    if (test_case.closestPointIndexResults[i] == std::numeric_limits<uint32_t>::max()) {
      indexResults.push_back(std::numeric_limits<uint32_t>::max());
    } else {
      indexResults.push_back(index);
    }
    pointResults.push_back(p);
    distanceResults.push_back(dist);
  }

  ASSERT_THAT(indexResults, testing::Pointwise(testing::Eq(), test_case.closestPointIndexResults));
  ASSERT_THAT(pointResults, testing::Pointwise(PointFuzzyEqual(), test_case.closestPointResults));
  ASSERT_THAT(distanceResults,
              testing::Pointwise(testing::DoubleEq(), test_case.closestPointDistanceResults));
}
