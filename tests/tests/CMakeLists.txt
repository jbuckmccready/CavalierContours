cmake_minimum_required(VERSION 3.11)
include (googletest.cmake)

add_executable(c_api_tests
  main.cpp
  cavc_pline_tests.cpp
  cavc_pline_function_tests.cpp
)

target_include_directories(c_api_tests
  PRIVATE include/)

target_link_libraries(c_api_tests
  PRIVATE PolylineFactory)

target_link_libraries(c_api_tests
  PRIVATE ${CAVC_C_API_LIB})

target_link_libraries(c_api_tests
  PRIVATE gtest_main gtest gmock)

# note: adding executable as single test rather than discovering tests via gtest_discover_tests
# because it adds parmeterized tests as individual tests and ctest runs them one by one which is
# quite slow
add_test(NAME c_api_tests COMMAND c_api_tests)
