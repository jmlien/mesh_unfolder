set(LIB_SRC
ArcLengthCurveSignature.cpp
CurveSignature.cpp
CurveCSS.cpp
img2ply.cpp
)

include_directories(
    "${CMAKE_SOURCE_DIR}/src"
)

add_library(CurveMatching ${LIB_SRC})
