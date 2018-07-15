project(libShapeOp)

set(CMAKE_CONFIGURATION_TYPES Release Debug CACHE TYPE INTERNAL FORCE)
  
set(OPENMP OFF CACHE BOOL "OpenMP")
set(SWIG_PYTHON OFF CACHE BOOL "SWIG Python")
set(SWIG_JAVA OFF CACHE BOOL "SWIG Java")
set(SWIG_CSHARP OFF CACHE BOOL "SWIG CSharp")

# Tell CMake he'll be able to include ".cmake" configurations 
# files in the folder where the current file is located
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

# Find Eigen3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
if(NOT EIGEN3_FOUND)
    message(ERROR " EIGEN not found.")
else()
    message("EIGEN found.")
endif() 

# Make a release library build by default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

message(${CMAKE_BUILD_TYPE})

# Where are the headers found?
include_directories("src/")
include_directories("api/")

# h files
file(GLOB SHAPEOP_HEADERS "src/*.h")
# cpp files
file(GLOB SHAPEOP_SOURCES "src/*.cpp")
# h files
file(GLOB API_HEADERS "api/*.h")
# cpp files
file(GLOB API_SOURCES "api/*.cpp")
# Shared library
set(HEADERS ${SHAPEOP_HEADERS} ${API_HEADERS})
set(SOURCES ${SHAPEOP_SOURCES} ${API_SOURCES})
add_library(ShapeOp SHARED ${HEADERS} ${SOURCES})

#Adding subdirectories
add_subdirectory(bindings)
