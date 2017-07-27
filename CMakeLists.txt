cmake_minimum_required(VERSION 2.6)

project(unfolder)

set(PACKAGE_VERSION "0.1.0")


########## add packages ########## 
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

if(NOT WIN32)
endif()

########## set compiler flags ##########
if(NOT WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
endif()

if(WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /D _CRT_SECURE_NO_WARNINGS")
endif()

########## add subdirectories ############
add_subdirectory(src)
add_subdirectory(libs/libga)
add_subdirectory(libs/itree)
