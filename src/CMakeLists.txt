#####################################################
#                COMMON Libraries                   #
#####################################################

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(OpenCV 3.0 REQUIRED)

include_directories(
    ${OPENGL_INCLUDE_DIRS}
    ${GLUT_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS})

link_directories(
    "${CMAKE_SOURCE_DIR}/libs"
    "${CMAKE_SOURCE_DIR}/libs/GL"
    "${CMAKE_SOURCE_DIR}/libs/libga"
    "${CMAKE_SOURCE_DIR}/libs/itree"
    "${CMAKE_SOURCE_DIR}/libs/glpk")

set(libunfolder_deps
  ga
  itree
  mathtool
  polygon
  CurveMatching
  model_graph
  ${OpenCV_LIBS}
  ${GLPK_LIBRARY})

if(APPLE)
    set(libunfolder_deps
      ${libunfolder_deps}
      libglpk.osx.a
      libRAPID.osx.a)
elseif(UNIX)
    set(libunfolder_deps
      ${libunfolder_deps}
      z
      gmp
      colamd
      dl
      rt
      glpk
      ltdl
      libRAPID.linux64.a)
endif()

if(WIN32)
    set(libunfolder_deps
      ${libunfolder_deps}
      RAPID.lib
      glut32.lib
      glpk.lib)
endif()

########## add subdirectories ############
#add_subdirectory(app)
add_subdirectory(GL)
add_subdirectory(mathtool)
add_subdirectory(modelgraph)
add_subdirectory(util)
add_subdirectory(polygon)
add_subdirectory(CurveMatching)

#############################################
# Libunfolder
#############################################

set(LIBUNFOLDER_SRC
    ${LIBUNFOLDER_SRC}
    model.cpp
    main.cpp
    net.cpp
    netset.cpp
    netsurgent.cpp
    MeshCutter.cpp
    Splitter.cpp
    CD.cpp
    ClusterUnfolding.cpp
    LPUnfolding.cpp
    spatial_hash.cpp
    objReader.cpp
    unfolder.cpp
    UnfoldingEvaluator.cpp
    UnfoldingProblem.cpp
    UnfoldingState.cpp
    OverlappingChecker.cpp
    CompactProblem.cpp)

add_library(libunfolder STATIC ${LIBUNFOLDER_SRC})
target_link_libraries(libunfolder ${libunfolder_deps})

#################################################################
# Unfolder GUI
#################################################################

set(APP_NAME unfolder)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

include_directories(
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}/src/GL"
    "${CMAKE_SOURCE_DIR}/src/mathtool"
    "${CMAKE_SOURCE_DIR}/src/modelgraph"
    "${CMAKE_SOURCE_DIR}/libs/glpk"
    "${CMAKE_SOURCE_DIR}/libs")

link_directories(
  "${CMAKE_SOURCE_DIR}/libs")

set(UNFOLDER_SRC
  ${UNFOLDER_SRC}
  main.cpp
)

add_executable(${APP_NAME} ${UNFOLDER_SRC})
add_dependencies(${APP_NAME} libunfolder)

target_link_libraries(${APP_NAME}
                    ${libunfolder_deps}
                    libunfolder
                    ${OPENGL_LIBRARIES}
                    ${GLUT_LIBRARY})
