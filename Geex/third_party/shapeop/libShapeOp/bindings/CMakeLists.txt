cmake_minimum_required(VERSION 2.8)

add_subdirectory(c)
add_subdirectory(c++)
if(SWIG_JAVA)
    add_subdirectory(java)
endif()
if(SWIG_PYTHON)
    add_subdirectory(python)
endif()
if(SWIG_CSHARP)
    add_subdirectory(csharp)
endif()

