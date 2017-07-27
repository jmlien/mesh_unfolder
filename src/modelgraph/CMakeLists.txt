set(LIB_SRC
    ModelEdge.cpp
    ModelGraph.cpp
    ModelNode.cpp)

include_directories(
    "${CMAKE_SOURCE_DIR}/src"
)

add_library(model_graph ${LIB_SRC})
