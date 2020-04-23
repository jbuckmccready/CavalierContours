#include "c_api_test_helpers.hpp"
#include "cavaliercontours.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <vector>

namespace t = testing;

struct CombinePlinesTestCase {
  std::string name;
  int combineMode;
  cavc_pline *plineA;
  cavc_pline *plineB;
  std::vector<PolylineProperties> expectedRemaining;
  std::vector<PolylineProperties> expectedSubtracted;

  CombinePlinesTestCase(std::string name, int combineMode,
                        std::vector<cavc_vertex> const &plineVertexesA,
                        std::vector<cavc_vertex> const &plineVertexesB,
                        std::vector<PolylineProperties> expectedRemaining,
                        std::vector<PolylineProperties> expectedSubtracted)
      : name(std::move(name)), combineMode(combineMode),
        plineA(
            cavc_pline_new(&plineVertexesA[0], static_cast<uint32_t>(plineVertexesA.size()), true)),
        plineB(
            cavc_pline_new(&plineVertexesA[0], static_cast<uint32_t>(plineVertexesB.size()), true)),
        expectedRemaining(std::move(expectedRemaining)),
        expectedSubtracted(std::move(expectedSubtracted)) {}

  CombinePlinesTestCase(std::string name, int combineMode, cavc_pline *plineA, cavc_pline *plineB,
                        std::vector<PolylineProperties> expectedRemaining,
                        std::vector<PolylineProperties> expectedSubtracted)
      : name(std::move(name)), combineMode(combineMode), plineA(plineA), plineB(plineB),

        expectedRemaining(std::move(expectedRemaining)),
        expectedSubtracted(std::move(expectedSubtracted)) {}
};

std::ostream &operator<<(std::ostream &os, CombinePlinesTestCase const &c) {
  os << "{ " << c.name << ", combineMode: " << c.combineMode << " }";
  return os;
}

static std::vector<CombinePlinesTestCase> createSimpleCases() {
  std::vector<CombinePlinesTestCase> cases;

  {
    // combining circle and rectangle
    std::vector<cavc_vertex> plineAVertexes = {{0, 1, 1}, {10, 1, 1}};
    std::vector<cavc_vertex> plineBVertexes = {{3, -10, 0}, {6, -10, 0}, {6, 10, 0}, {3, 10, 0}};
    cavc_pline *plineA = plineFromVertexes(plineAVertexes, true);
    cavc_pline *plineB = plineFromVertexes(plineBVertexes, true);
    // Union
    std::vector<PolylineProperties> expectedRemaining;
    expectedRemaining.emplace_back(10, 109.15381629282, 52.324068506275, 0, -10, 10, 10);
    std::vector<PolylineProperties> expectedSubtracted;
    cases.emplace_back("circle_rectangle_union", 0, plineA, plineB, expectedRemaining,
                       expectedSubtracted);

    // Exclude
    expectedRemaining.clear();
    expectedSubtracted.clear();
    expectedRemaining.emplace_back(3, 19.816835628274, 20.757946197186, 0, -4, 3, 5.5825756949558);
    expectedRemaining.emplace_back(3, 29.336980664548, 23.492343031178, 6, -3.8989794855664, 10, 6);
    cases.emplace_back("circle_rectangle_exclude", 1, plineA, plineB, expectedRemaining,
                       expectedSubtracted);

    // Intersect
    expectedRemaining.clear();
    expectedSubtracted.clear();
    expectedRemaining.emplace_back(4, 29.386000046924, 25.091858029623, 3, -4, 6, 6);
    cases.emplace_back("circle_rectangle_intersect", 2, plineA, plineB, expectedRemaining,
                       expectedSubtracted);

    // XOR
    expectedRemaining.clear();
    expectedRemaining.emplace_back(3, 19.816835628274, 20.757946197186, 0, -4, 3, 5.5825756949558);
    expectedRemaining.emplace_back(4, -18.306999976538, 18.582818653767, 3, -10, 6,
                                   -3.5825756949558);
    expectedRemaining.emplace_back(3, 29.336980664548, 23.492343031178, 6, -3.8989794855664, 10, 6);
    expectedRemaining.emplace_back(4, -12.306999976538, 14.582818653767, 3, 5.5825756949558, 6, 10);
    expectedSubtracted.clear();
    cases.emplace_back("circle_rectangle_xor", 3, plineA, plineB, expectedRemaining,
                       expectedSubtracted);
  }

  return cases;
}

static std::vector<CombinePlinesTestCase> simpleCases = createSimpleCases();

class cavc_combine_plinesTests : public t::TestWithParam<CombinePlinesTestCase> {};

INSTANTIATE_TEST_SUITE_P(simple_cases, cavc_combine_plinesTests, t::ValuesIn(simpleCases));

TEST_P(cavc_combine_plinesTests, combine_plines_test) {
  CombinePlinesTestCase const &testCase = GetParam();
  cavc_pline_list *remaining = nullptr;
  cavc_pline_list *subtracted = nullptr;
  cavc_combine_plines(testCase.plineA, testCase.plineB, testCase.combineMode, &remaining,
                      &subtracted);

  ASSERT_EQ(cavc_pline_list_count(remaining), testCase.expectedRemaining.size());
  ASSERT_EQ(cavc_pline_list_count(subtracted), testCase.expectedSubtracted.size());

  std::vector<PolylineProperties> remainingProperties;
  remainingProperties.reserve(testCase.expectedRemaining.size());
  for (uint32_t i = 0; i < testCase.expectedRemaining.size(); ++i) {
    cavc_pline *pline = cavc_pline_list_get(remaining, i);
    remainingProperties.emplace_back(pline);
  }
  ASSERT_THAT(remainingProperties, t::UnorderedPointwise(t::Eq(), testCase.expectedRemaining));

  std::vector<PolylineProperties> subtractedProperties;
  subtractedProperties.reserve(testCase.expectedSubtracted.size());
  for (uint32_t i = 0; i < testCase.expectedSubtracted.size(); ++i) {
    cavc_pline *pline = cavc_pline_list_get(subtracted, i);
    subtractedProperties.emplace_back(pline);
  }
  ASSERT_THAT(subtractedProperties, t::UnorderedPointwise(t::Eq(), testCase.expectedSubtracted));

  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);
}

TEST(cavc_combine_plinesTests, combine_with_self_invariants) {
  std::vector<cavc_vertex> plineVertexes = {{27.554688, 1, 0},    {27.554688, 0.75, 0.414214},
                                            {27.804688, 0.5, 0},  {32.195313, 0.5, 0.414214},
                                            {32.445313, 0.75, 0}, {32.445313, 1, 0.414214},
                                            {32.195313, 1.25, 0}, {31.5, 1.25, -0.414214},
                                            {31, 1.75, 0},        {29, 1.75, -0.414214},
                                            {28.5, 1.25, 0},      {27.804688, 1.25, 0.414214}};
  std::vector<cavc_vertex> revPlineVertexes = plineVertexes;
  reverseDirection(revPlineVertexes);
  cavc_pline *pline = plineFromVertexes(plineVertexes, true);
  cavc_pline *revPline = plineFromVertexes(revPlineVertexes, true);

  cavc_pline_list *remaining = nullptr;
  cavc_pline_list *subtracted = nullptr;

  // test union with self is always self
  cavc_combine_plines(pline, pline, 0, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 1);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  cavc_pline *combineResult = cavc_pline_list_get(remaining, 0);
  std::vector<cavc_vertex> resultVertexes(cavc_pline_vertex_count(combineResult));
  cavc_pline_vertex_data(combineResult, &resultVertexes[0]);
  ASSERT_THAT(resultVertexes, t::Pointwise(VertexEqual(), plineVertexes));
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);
  // reversed
  cavc_combine_plines(revPline, revPline, 0, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 1);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  combineResult = cavc_pline_list_get(remaining, 0);
  resultVertexes.clear();
  resultVertexes.resize(cavc_pline_vertex_count(combineResult));
  cavc_pline_vertex_data(combineResult, &resultVertexes[0]);
  ASSERT_THAT(resultVertexes, t::Pointwise(VertexEqual(), revPlineVertexes));
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);

  // test exclude with self is always empty
  cavc_combine_plines(pline, pline, 1, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 0);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);
  // reversed
  cavc_combine_plines(revPline, revPline, 1, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 0);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);
  // reversed and not reversed
  cavc_combine_plines(revPline, pline, 1, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 0);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);
  cavc_combine_plines(pline, revPline, 1, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 0);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);

  // test intersect with self is always self
  cavc_combine_plines(pline, pline, 2, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 1);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  combineResult = cavc_pline_list_get(remaining, 0);
  resultVertexes.clear();
  resultVertexes.resize(cavc_pline_vertex_count(combineResult));
  cavc_pline_vertex_data(combineResult, &resultVertexes[0]);
  ASSERT_THAT(resultVertexes, t::Pointwise(VertexEqual(), plineVertexes));
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);
  // reversed
  cavc_combine_plines(revPline, revPline, 2, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 1);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  combineResult = cavc_pline_list_get(remaining, 0);
  resultVertexes.clear();
  resultVertexes.resize(cavc_pline_vertex_count(combineResult));
  cavc_pline_vertex_data(combineResult, &resultVertexes[0]);
  ASSERT_THAT(resultVertexes, t::Pointwise(VertexEqual(), revPlineVertexes));
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);

  // test XOR with self is always empty
  cavc_combine_plines(pline, pline, 3, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 0);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);
  // reversed
  cavc_combine_plines(revPline, revPline, 3, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 0);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);
  // reversed and not reversed
  cavc_combine_plines(revPline, pline, 3, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 0);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);
  cavc_combine_plines(pline, revPline, 3, &remaining, &subtracted);
  ASSERT_EQ(cavc_pline_list_count(remaining), 0);
  ASSERT_EQ(cavc_pline_list_count(subtracted), 0);
  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);

  cavc_pline_delete(pline);
  cavc_pline_delete(revPline);
}

class cavc_combine_plinesNoModifySuite : public t::TestWithParam<CombinePlinesTestCase> {};
INSTANTIATE_TEST_SUITE_P(cavc_circle_rect_cases, cavc_combine_plinesNoModifySuite,
                         t::ValuesIn({simpleCases[0], simpleCases[1], simpleCases[2],
                                      simpleCases[3]}));

TEST_P(cavc_combine_plinesNoModifySuite, combine_plines_does_not_modify_input) {
  CombinePlinesTestCase const &testCase = GetParam();
  std::vector<cavc_vertex> vertexesBeforeA(cavc_pline_vertex_count(testCase.plineA));
  cavc_pline_vertex_data(testCase.plineA, &vertexesBeforeA[0]);
  std::vector<cavc_vertex> vertexesBeforeB(cavc_pline_vertex_count(testCase.plineB));
  cavc_pline_vertex_data(testCase.plineB, &vertexesBeforeB[0]);

  cavc_pline_list *remaining;
  cavc_pline_list *subtracted;
  cavc_combine_plines(testCase.plineA, testCase.plineB, testCase.combineMode, &remaining,
                      &subtracted);

  std::vector<cavc_vertex> vertexesAfterA(cavc_pline_vertex_count(testCase.plineA));
  cavc_pline_vertex_data(testCase.plineA, &vertexesAfterA[0]);
  ASSERT_THAT(vertexesAfterA, t::Pointwise(VertexEqual(), vertexesBeforeA));

  std::vector<cavc_vertex> vertexesAfterB(cavc_pline_vertex_count(testCase.plineB));
  cavc_pline_vertex_data(testCase.plineB, &vertexesAfterB[0]);
  ASSERT_THAT(vertexesAfterB, t::Pointwise(VertexEqual(), vertexesBeforeB));

  cavc_pline_list_delete(remaining);
  cavc_pline_list_delete(subtracted);
}

int main(int argc, char **argv) {
  t::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
