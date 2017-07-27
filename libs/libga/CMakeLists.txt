set(LIB_SRC
    Individual.cpp
    Problem.cpp
    Species.cpp
    Breeder.cpp
    util/StreamHelper.cpp
    util/Random.cpp
)

include_directories(
    "${CMAKE_SOURCE_DIR}/libga"
)

add_library(ga ${LIB_SRC})