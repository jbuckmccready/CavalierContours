add_library(PolylineFactory src/polylinefactory.cpp)

target_include_directories(PolylineFactory
  PUBLIC include/)

target_link_libraries(PolylineFactory
  PUBLIC ${CAVC_CPP_HEADER_ONLY_LIB})

target_link_libraries(PolylineFactory
  PUBLIC ${CAVC_C_API_LIB})
