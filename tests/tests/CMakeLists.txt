cmake_minimum_required(VERSION 3.11)
include (googletest.cmake)

macro(cavc_add_test name)
  add_executable(${name} ${name}.cpp)

  target_include_directories(${name}
    PRIVATE include/)

  target_link_libraries(${name}
    PRIVATE PolylineFactory)

  target_link_libraries(${name}
    PRIVATE ${CAVC_C_API_LIB})

  target_link_libraries(${name}
    PRIVATE gtest_main gtest gmock)

  # note: adding executable as single test rather than discovering tests via gtest_discover_tests
  # because it adds parmeterized tests as individual tests and ctest runs them one by one which is
  # quite slow
  add_test(NAME ${name} COMMAND ${name})
endmacro()

cavc_add_test(cavc_pline_tests)
cavc_add_test(cavc_pline_function_tests)
cavc_add_test(cavc_parallel_offset_tests)
cavc_add_test(cavc_combine_plines_tests)
cavc_add_test(staticspatialindex_tests)
