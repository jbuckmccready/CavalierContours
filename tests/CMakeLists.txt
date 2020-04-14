cmake_minimum_required(VERSION 3.11)
include (googletest.cmake)

add_executable(c_api_tests main.cpp c_api_tests.cpp)
target_link_libraries(c_api_tests
  PRIVATE ${CAVC_C_API_LIB})

target_link_libraries(c_api_tests
  PRIVATE gtest_main gtest gmock)

add_test(NAME c_api_tests COMMAND c_api_tests)
