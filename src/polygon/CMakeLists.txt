set(LIB_SRC
    bpc.cpp
    chull.cpp
    diagonal2.cpp
    dp_approx.cpp
    poly_approx.cpp
    polygon.cpp
    polyline.cpp
    net_analyzer.cpp
    bbox2d.cpp
)

include_directories(
    "${CMAKE_SOURCE_DIR}/src"
)

add_library(polygon ${LIB_SRC})
