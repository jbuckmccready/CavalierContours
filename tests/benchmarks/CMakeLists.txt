cmake_minimum_required(VERSION 3.11)

include(clipper.cmake)
find_package(benchmark REQUIRED)

macro(add_benchmark name)
    add_executable(${name} ${name}.cpp)

    target_include_directories(${name}
    PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/include
        ${benchmark_INCLUDE_DIRS})

    target_link_libraries(${name}
    PRIVATE
        ${CAVC_CPP_HEADER_ONLY_LIB}
        PolylineFactory
        benchmark::benchmark)

    if (MSVC)
        target_link_options(${name} PRIVATE $<$<CONFIG:RELWITHDEBINFO>:/PROFILE>)
    endif()
endmacro()

add_benchmark(offsetbenchmarks)
add_benchmark(spatialindexbenchmarks)
add_benchmark(extentsbenchmarks)
add_benchmark(areabenchmarks)
add_benchmark(pathlengthbenchmarks)
add_benchmark(windingnumberbenchmarks)
add_benchmark(combinebenchmarks)

add_benchmark(clipperbenchmarks)
target_link_libraries(clipperbenchmarks
  PRIVATE clipper_static)
