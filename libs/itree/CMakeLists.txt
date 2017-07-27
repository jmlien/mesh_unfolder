set(LIB_SRC
    EndPoint.cpp
    Interval.cpp
    RectKD.cpp
)

include_directories(
    "${CMAKE_SOURCE_DIR}/itree"
)

add_library(itree ${LIB_SRC})
