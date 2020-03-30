cmake_minimum_required(VERSION 3.5)

add_library(CavalierContours INTERFACE)
target_include_directories(CavalierContours
  INTERFACE include/)

add_subdirectory(examples)