macro(add_example name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name}
    PRIVATE ${CAVC_CPP_HEADER_ONLY_LIB})
endmacro()

add_example(polylinecombine)
add_example(polylineoffset)
add_example(windingnumber)
add_example(basicpolylinefunctions)
add_example(polylineoffsetislands)
