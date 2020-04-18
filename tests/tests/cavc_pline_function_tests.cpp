#include "c_api_test_helpers.hpp"
#include "cavaliercontours.h"
#include <iostream>
#include <string>
#include <vector>
struct cavc_plineFunctionsTestCase {

  cavc_plineFunctionsTestCase(std::string name, std::vector<cavc_vertex> vertexes, bool isClosed)
      : name(std::move(name)),
        pline(cavc_pline_new(&vertexes[0], static_cast<uint32_t>(vertexes.size()), isClosed)),
        plineVertexes(std::move(vertexes)) {}

  // simple name for the test case
  std::string name;

  // polyline for the test case
  cavc_pline *pline = nullptr;
  // vertexes of pline for easy access
  std::vector<cavc_vertex> plineVertexes;
  // helper method to create cavc_pline and set vertexes
  void setPline(std::vector<cavc_vertex> vertexes, bool isClosed) {
    plineVertexes = vertexes;
    if (pline) {
      cavc_pline_set_vertex_data(pline, &vertexes[0], static_cast<uint32_t>(vertexes.size()));
      cavc_pline_set_is_closed(pline, isClosed);
      return;
    }
    pline = cavc_pline_new(&vertexes[0], static_cast<uint32_t>(vertexes.size()), isClosed);
  }

  // expected signed area
  cavc_real signedArea = std::numeric_limits<cavc_real>::quiet_NaN();
  bool skipAreaTest() const { return std::isnan(signedArea); }

  // expected path length
  cavc_real pathLength = std::numeric_limits<cavc_real>::quiet_NaN();
  bool skipPathLengthTest() const { return std::isnan(signedArea); }

  // expected extents
  cavc_real minX = std::numeric_limits<double>::quiet_NaN();
  cavc_real minY = std::numeric_limits<double>::quiet_NaN();
  cavc_real maxX = std::numeric_limits<double>::quiet_NaN();
  cavc_real maxY = std::numeric_limits<double>::quiet_NaN();
  bool skipExtentsTest() const {
    return std::isnan(minX) && std::isnan(minY) && std::isnan(maxX) && std::isnan(maxY);
  }

  // winding number test input points
  std::vector<cavc_point> windingNumberTestPts;
  // expected winding number results for test pts
  std::vector<int> windingNumberResults;
  bool skipWindingNumberTest() const { return windingNumberTestPts.empty(); }

  // add a winding number test input point
  void addWindingNumberTestPt(cavc_point testPt, int result) {
    windingNumberTestPts.push_back(testPt);
    windingNumberResults.push_back(result);
  }

  // closest point test input points
  std::vector<cavc_point> closestPointTestPts;
  // expected results from closest point calculation
  std::vector<uint32_t> closestPointIndexResults;
  std::vector<cavc_point> closestPointResults;
  std::vector<cavc_real> closestPointDistanceResults;
  bool skipClosestPointTest() const { return closestPointTestPts.empty(); }

  // add a closestPoint test input point
  void addClosestPointTestPt(cavc_point testPt, cavc_point closestPointResult,
                             cavc_real distanceResult,
                             uint32_t indexResult = std::numeric_limits<uint32_t>::max()) {
    closestPointTestPts.push_back(testPt);
    closestPointResults.push_back(closestPointResult);
    closestPointDistanceResults.push_back(distanceResult);
    closestPointIndexResults.push_back(indexResult);
  }

  std::vector<cavc_real> offsetTestDeltas;
  std::vector<std::vector<cavc_vertex>> resultOffsetPlines;
  std::vector<bool> resultIsClosed;
  bool skipOffsetTest() const { return offsetTestDeltas.empty(); }

  // add an offset test input
  void addOffsetTest(cavc_real offsetDelta, std::vector<cavc_vertex> vertexes, bool isClosed) {
    offsetTestDeltas.push_back(offsetDelta);
    resultOffsetPlines.emplace_back(std::move(vertexes));
    resultIsClosed.push_back(isClosed);
  }

  std::vector<cavc_real> collapsedOffsetDeltas;
  void addCollapsedOffsetTest(cavc_real offsetDelta) {
    collapsedOffsetDeltas.push_back(offsetDelta);
  }
  bool skipCollapsedOffsetTest() const { return collapsedOffsetDeltas.empty(); }

  bool isClosed() const { return cavc_pline_is_closed(pline); }
};

std::ostream &operator<<(std::ostream &os, cavc_plineFunctionsTestCase const &c) {
  os << c.name;
  return os;
}

void addCircleCases(std::vector<cavc_plineFunctionsTestCase> &cases, cavc_real circleRadius,
                    cavc_point circleCenter, bool reverse) {
  cavc_real expectedCircleArea = PI() * circleRadius * circleRadius;
  cavc_real expectedCircleLength = 2 * PI() * circleRadius;

  auto createName = [&](auto n) {
    std::stringstream ss;
    ss << n;
    if (reverse) {
      ss << "_rev";
    }
    ss << " (radius: " << circleRadius << ", center: " << circleCenter << ")";
    return ss.str();
  };

  auto addCase = [&](std::string name, std::vector<cavc_vertex> const &vertexes, int direction) {
    cavc_plineFunctionsTestCase test_case(createName(name), vertexes, true);
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
      cavc_real xCos = std::cos(PI() / 4 + i * PI() / 2);
      cavc_real ySin = std::sin(PI() / 4 + i * PI() / 2);

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

    cavc_real offsetOutwardDelta = -direction * 0.25 * circleRadius;
    std::vector<cavc_vertex> offsetOutwardVertexes;

    cavc_real offsetInwardDelta = direction * 0.5 * circleRadius;
    std::vector<cavc_vertex> offsetInwardVertexes;

    for (auto const &v : test_case.plineVertexes) {
      cavc_point dirVec = {v.x - circleCenter.x, v.y - circleCenter.y};
      cavc_real dirVecLength = std::sqrt(dirVec.x * dirVec.x + dirVec.y * dirVec.y);
      cavc_point unitDirVec = {dirVec.x / dirVecLength, dirVec.y / dirVecLength};
      cavc_real outwardMagnitude = circleRadius + std::abs(offsetOutwardDelta);
      cavc_real inwardMagnitude = circleRadius - std::abs(offsetInwardDelta);
      offsetOutwardVertexes.push_back({outwardMagnitude * unitDirVec.x + circleCenter.x,
                                       outwardMagnitude * unitDirVec.y + circleCenter.y, v.bulge});
      offsetInwardVertexes.push_back({inwardMagnitude * unitDirVec.x + circleCenter.x,
                                      inwardMagnitude * unitDirVec.y + circleCenter.y, v.bulge});
    }
    test_case.addOffsetTest(offsetOutwardDelta, std::move(offsetOutwardVertexes), true);
    test_case.addOffsetTest(offsetInwardDelta, std::move(offsetInwardVertexes), true);

    // offsetting inward by full radius or more should result in no offsets
    test_case.addCollapsedOffsetTest(direction * circleRadius);
    test_case.addCollapsedOffsetTest(direction * 1.5 * circleRadius);
    test_case.addCollapsedOffsetTest(direction * 2.0 * circleRadius);

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

  vertexes = {{circleCenter.x + circleRadius * std::cos(PI() / 4),
               circleCenter.y + circleRadius * std::sin(PI() / 4), 1},
              {circleCenter.x + circleRadius * std::cos(5 * PI() / 4),
               circleCenter.y + circleRadius * std::sin(5 * PI() / 4), 1}};
  if (reverse) {
    std::reverse(vertexes.begin(), vertexes.end());
  }
  addCase("ccw_circle_not_axis_aligned", vertexes, 1);

  for (auto &v : vertexes) {
    v = {v.x, v.y, -v.bulge};
  }
  addCase("cw_circle_not_axis_aligned", vertexes, -1);
}

std::vector<cavc_plineFunctionsTestCase> createCircleCases() {
  std::vector<cavc_plineFunctionsTestCase> result;
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

class cavc_plineFunctionTests : public testing::TestWithParam<cavc_plineFunctionsTestCase> {
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
  static std::vector<cavc_plineFunctionsTestCase> circleCases;
};
void cavc_plineFunctionTests::SetUp() {}
void cavc_plineFunctionTests::TearDown() {}

std::vector<cavc_plineFunctionsTestCase> cavc_plineFunctionTests::circleCases = createCircleCases();

INSTANTIATE_TEST_SUITE_P(cavc_pline_circles, cavc_plineFunctionTests,
                         testing::ValuesIn(cavc_plineFunctionTests::circleCases));

TEST_P(cavc_plineFunctionTests, cavc_get_path_length) {
  cavc_plineFunctionsTestCase const &test_case = GetParam();
  if (test_case.skipPathLengthTest()) {
    GTEST_SKIP();
  }
  ASSERT_NEAR(cavc_get_path_length(test_case.pline), test_case.pathLength, TEST_EPSILON());
}

TEST_P(cavc_plineFunctionTests, cavc_get_area) {
  cavc_plineFunctionsTestCase const &test_case = GetParam();
  if (test_case.skipAreaTest()) {
    GTEST_SKIP();
  }
  ASSERT_NEAR(cavc_get_area(test_case.pline), test_case.signedArea, TEST_EPSILON());
}

TEST_P(cavc_plineFunctionTests, cavc_get_winding_number) {
  cavc_plineFunctionsTestCase const &test_case = GetParam();
  if (test_case.skipWindingNumberTest()) {
    GTEST_SKIP();
  }
  std::vector<int> windingNumberResults;
  windingNumberResults.reserve(test_case.windingNumberTestPts.size());
  for (auto const &pt : test_case.windingNumberTestPts) {
    windingNumberResults.push_back(cavc_get_winding_number(test_case.pline, pt));
  }

  ASSERT_THAT(windingNumberResults,
              testing::Pointwise(testing::Eq(), test_case.windingNumberResults));
}

TEST_P(cavc_plineFunctionTests, cavc_get_extents) {
  cavc_plineFunctionsTestCase const &test_case = GetParam();
  if (test_case.skipExtentsTest()) {
    GTEST_SKIP();
  }
  cavc_real minX;
  cavc_real minY;
  cavc_real maxX;
  cavc_real maxY;
  cavc_get_extents(test_case.pline, &minX, &minY, &maxX, &maxY);
  EXPECT_NEAR(minX, test_case.minX, TEST_EPSILON());
  EXPECT_NEAR(minY, test_case.minY, TEST_EPSILON());
  EXPECT_NEAR(maxX, test_case.maxX, TEST_EPSILON());
  EXPECT_NEAR(maxY, test_case.maxY, TEST_EPSILON());
}

TEST_P(cavc_plineFunctionTests, cavc_get_closest_point) {
  cavc_plineFunctionsTestCase const &test_case = GetParam();
  if (test_case.skipClosestPointTest()) {
    GTEST_SKIP();
  }
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
  ASSERT_THAT(distanceResults, testing::Pointwise(testing::DoubleNear(TEST_EPSILON()),
                                                  test_case.closestPointDistanceResults));
}

TEST_P(cavc_plineFunctionTests, cavc_parallel_offet) {
  cavc_plineFunctionsTestCase const &test_case = GetParam();
  if (test_case.skipOffsetTest() && test_case.skipCollapsedOffsetTest()) {
    GTEST_SKIP();
  }

  if (!test_case.skipOffsetTest()) {
    std::vector<cavc_pline_list *> results(test_case.offsetTestDeltas.size());
    for (std::size_t i = 0; i < test_case.offsetTestDeltas.size(); ++i) {
      cavc_parallel_offet(test_case.pline, test_case.offsetTestDeltas[i], &results[i], 0);
    }

    // test there is only one resulting offset pline
    std::vector<uint32_t> resultCounts;
    for (auto const &result : results) {
      resultCounts.push_back(cavc_pline_list_count(result));
    }
    ASSERT_THAT(resultCounts, testing::Each(1));

    // test is closed is correct
    std::vector<bool> isClosedResults;
    for (auto const &result : results) {
      isClosedResults.push_back(cavc_pline_is_closed(cavc_pline_list_get(result, 0)));
    }
    ASSERT_THAT(isClosedResults, testing::Pointwise(testing::Eq(), test_case.resultIsClosed));

    // test all vertex values
    std::vector<std::vector<cavc_vertex>> vertexResults;
    for (auto const &result : results) {
      cavc_pline *pline = cavc_pline_list_get(result, 0);
      vertexResults.push_back(std::vector<cavc_vertex>(cavc_pline_vertex_count(pline)));
      cavc_pline_vertex_data(pline, &vertexResults.back()[0]);
    }

    ASSERT_THAT(vertexResults,
                testing::Pointwise(VertexListsFuzzyEqual(), test_case.resultOffsetPlines));

    for (auto result : results) {
      cavc_pline_list_delete(result);
    }
  }

  if (!test_case.skipCollapsedOffsetTest()) {
    std::vector<cavc_pline_list *> results(test_case.collapsedOffsetDeltas.size());
    for (std::size_t i = 0; i < test_case.collapsedOffsetDeltas.size(); ++i) {
      cavc_parallel_offet(test_case.pline, test_case.collapsedOffsetDeltas[i], &results[i], 0);
    }
    ASSERT_THAT(results, testing::Each(
                             testing::Truly([](auto r) { return cavc_pline_list_count(r) == 0; })));
    for (auto result : results) {
      cavc_pline_list_delete(result);
    }
  }
}

TEST_P(cavc_plineFunctionTests, cavc_combine_plines) {
  cavc_plineFunctionsTestCase const &test_case = GetParam();
  cavc_pline_list *remaining = nullptr;
  cavc_pline_list *subtracted = nullptr;

  // test union with self is always self
  cavc_combine_plines(test_case.pline, test_case.pline, 0, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 1);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  cavc_pline *combineResult = cavc_pline_list_get(remaining, 0);
  std::vector<cavc_vertex> resultVertexes(cavc_pline_vertex_count(combineResult));
  cavc_pline_vertex_data(combineResult, &resultVertexes[0]);
  ASSERT_EQ(resultVertexes.size(), test_case.plineVertexes.size());
  ASSERT_THAT(resultVertexes, testing::Pointwise(VertexEqual(), test_case.plineVertexes));
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);

  // test exclude with self is always empty
  cavc_combine_plines(test_case.pline, test_case.pline, 1, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 0);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);

  // test intersect with self is always self
  cavc_combine_plines(test_case.pline, test_case.pline, 2, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 1);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  combineResult = cavc_pline_list_get(remaining, 0);
  resultVertexes.resize(cavc_pline_vertex_count(combineResult));
  cavc_pline_vertex_data(combineResult, &resultVertexes[0]);
  ASSERT_EQ(resultVertexes.size(), test_case.plineVertexes.size());
  ASSERT_THAT(resultVertexes, testing::Pointwise(VertexEqual(), test_case.plineVertexes));
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);

  // test XOR with self is always empty
  cavc_combine_plines(test_case.pline, test_case.pline, 3, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 0);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
