cmake_minimum_required(VERSION 2.8)
project(ShapeOp)

# Configure a header file to pass some of the CMake settings to the Doxygen documentation
configure_file("${PROJECT_SOURCE_DIR}/doc/DoxyShapeOp.in" "${PROJECT_BINARY_DIR}/doc/DoxyShapeOp")
configure_file("${PROJECT_SOURCE_DIR}/doc/DoxyShapeOp.bib" "${PROJECT_BINARY_DIR}/doc/DoxyShapeOp.bib" COPYONLY)

# Compiler flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    message("Clang.")
    set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "-O4 -DNDEBUG")
    add_definitions(-DSHAPEOP_CLANG)
    add_definitions(-DSHAPEOP_EXPORT)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("GCC.")
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
    add_definitions(-DSHAPEOP_GCC)
    add_definitions(-DSHAPEOP_EXPORT)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    message("ICC.")
    if(WIN32)
        set(CMAKE_CXX_FLAGS "/Qstd=c++11")
        set(CMAKE_CXX_FLAGS_RELEASE "/O3 /DNDEBUG")
        add_definitions(/DSHAPEOP_ICC)
        add_definitions(/DSHAPEOP_EXPORT)
    else()
        set(CMAKE_CXX_FLAGS "-std=c++11")
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
        add_definitions(-DSHAPEOP_ICC)
        add_definitions(-DSHAPEOP_EXPORT)
    endif()
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message("MSVC.")
    set(CMAKE_CXX_FLAGS_DEBUG "/MTd")
    set(CMAKE_CXX_FLAGS_RELEASE "/MT /Ox /Ot /fp:fast /GS- /GL /DNDEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG /INCREMENTAL:NO /OPT:REF")
    add_definitions(/DSHAPEOP_MSVC)
    add_definitions(/DSHAPEOP_EXPORT)
endif()

#OpenMP
if(OPENMP)
  FIND_PACKAGE(OpenMP QUIET)
  if(OPENMP_FOUND)
      message("OpenMP found. OpenMP activated in release.")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${OpenMP_EXE_LINKER_FLAGS}")
      add_definitions(-DSHAPEOP_OPENMP)
  else()
      message("OpenMP not found.")
  endif()
endif()

add_subdirectory(libShapeOp)
add_subdirectory(applications)
