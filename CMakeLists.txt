cmake_minimum_required(VERSION 3.21)

# detect if CavalierContours is a subproject or not
# do not include examples, tests, and benchmarks in that case
if (NOT DEFINED PROJECT_NAME)
  set(NOT_SUBPROJECT ON)
endif()

project(CavalierContours VERSION 0.1)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Allow super-projects to override options
if (POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW)
endif ()

# Build options
option(CAVC_HEADER_ONLY "C++ header only library of CavalierContours, if ON then no library is built" OFF)
option(CAVC_BUILD_SHARED_LIB "Build the C API CavalierContours dynamic shared library (SET OFF for static library)" ON)

if(CAVC_HEADER_ONLY)
  set(CAVC_CPP_HEADER_ONLY_LIB ${PROJECT_NAME})
else()
  set(CAVC_CPP_HEADER_ONLY_LIB ${PROJECT_NAME}_header_only)
  set(CAVC_C_API_LIB ${PROJECT_NAME})
endif()

add_library(${CAVC_CPP_HEADER_ONLY_LIB} INTERFACE)
target_include_directories(${CAVC_CPP_HEADER_ONLY_LIB}
  INTERFACE include/
  ${CMAKE_SOURCE_DIR}/build/include
  )

if(NOT CAVC_HEADER_ONLY)
  if(CAVC_BUILD_SHARED_LIB)
    add_library(${CAVC_C_API_LIB}
      SHARED src/cavaliercontours.cpp)
  else()
    add_library(${CAVC_C_API_LIB}
      STATIC src/cavaliercontours.cpp)
    target_compile_definitions(${CAVC_C_API_LIB}
      PUBLIC CAVC_STATIC_LIB)
  endif()
  # compiler warnings, note: not using /W4 on MSVC since it emits warning due to default of /W3
  # from cmake
  if(MSVC)
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
      target_compile_options(${CAVC_C_API_LIB} PRIVATE
        /clang:-Werror /clang:-Wall /clang:-Wextra /clang:-pedantic-errors /clang:-Wsign-conversion /clang:-Wconversion)
    endif()
    target_compile_options(${CAVC_C_API_LIB} PRIVATE
      /WX)
  else()
    target_compile_options(${CAVC_C_API_LIB} PRIVATE
      -Werror -Wall -Wextra -pedantic-errors -Wsign-conversion -Wconversion)
  endif()
  target_compile_definitions(${CAVC_C_API_LIB}
    PRIVATE CAVC_EXPORTS)
  target_include_directories(${CAVC_C_API_LIB}
    PUBLIC c_api_include/
    )
  target_link_libraries(${CAVC_C_API_LIB}
    PRIVATE ${CAVC_CPP_HEADER_ONLY_LIB})
endif()



# Copy Headers to build/include/cavc
set(CAVC_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/build/include/cavc)
set(CAVC_C_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/build/include/cavc_c)
set(CAVC_INTERNAL_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/build/include/cavc/internal)

if (EXISTS ${CAVC_INCLUDE_DIR})
    file(REMOVE_RECURSE ${CAVC_INCLUDE_DIR})
endif()

set(CAVC_HEADERS
    include/intrlineseg2circle2.hpp
    include/intrlineseg2lineseg2.hpp
    include/mathutils.hpp
    include/intrcircle2circle2.hpp
    include/vector.hpp
    include/plinesegment.hpp
    include/polyline.hpp
    include/polylinecombine.hpp
    include/polylineintersects.hpp
    include/polylineoffset.hpp
    include/polylineoffsetislands.hpp
    include/staticspatialindex.hpp
    include/vector2.hpp
)

foreach(CAVC_HEADER ${CAVC_HEADERS})
    file(COPY ${CAVC_HEADER} DESTINATION ${CAVC_INCLUDE_DIR})
endforeach()

set(CAVC_INTERNAL_HEADERS
    include/internal/diagnostics.hpp
    include/internal/common.hpp
)

foreach(CAVC_INTERNAL_HEADER ${CAVC_INTERNAL_HEADERS})
    file(COPY ${CAVC_INTERNAL_HEADER} DESTINATION ${CAVC_INTERNAL_INCLUDE_DIR})
endforeach()
file(COPY src/cavaliercontours.h DESTINATION ${CAVC_C_INCLUDE_DIR})

if (NOT_SUBPROJECT AND NOT CAVC_HEADER_ONLY)
  # output all libraries and executables to the root build dir for simplicity in running tests with dynamic libs
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  add_subdirectory(examples)
  find_package(GTest REQUIRED)
  include(CTest)
  include(GoogleTest)
  enable_testing()
  add_subdirectory(tests)
endif()
