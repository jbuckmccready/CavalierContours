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
