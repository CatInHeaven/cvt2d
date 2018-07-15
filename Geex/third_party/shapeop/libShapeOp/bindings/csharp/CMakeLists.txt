find_package(SWIG QUIET)
if(SWIG_FOUND)
    message("Swig found. Building CSharp bindings.")
    project(libShapeOpCSharp)
    include(${SWIG_USE_FILE})

    # Tell CMake he'll be able to include ".cmake" configurations
    # files in the folder where the current file is located
    set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

    find_package(CSharp)
    include_directories(${CSHARP_INCLUDE_PATH})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})

    set(CMAKE_SWIG_FLAGS "")

    set_source_files_properties(shapeop.i PROPERTIES CPLUSPLUS ON)
    set_source_files_properties(shapeop.i PROPERTIES SWIG_FLAGS "-includeall")

    swig_add_module(shapeopCSharp csharp shapeop.i ${SOURCES})
    swig_link_libraries(shapeopCSharp ${CSHARP_LIBRARIES})

    configure_file(runme.cs runme.cs COPYONLY)
else(SWIG_FOUND)
        message("Swig not found. Not building CSharp bindings.")
endif(SWIG_FOUND)
