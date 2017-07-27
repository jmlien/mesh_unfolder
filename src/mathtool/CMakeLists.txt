set(LIB_SRC
    Basic.cpp
    Quaternion.cpp
    Geometry.cpp
    Box.cpp)

include_directories(
    "${CMAKE_SOURCE_DIR}/src"
)

add_library(mathtool ${LIB_SRC})