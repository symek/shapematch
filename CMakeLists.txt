# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.4 )
project( shapematch )

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Debug
            CACHE STRING "Choose the type of build, options are: ${CMAKE_CONFIGURATION_TYPES}" FORCE)
    set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${CMAKE_CONFIGURATION_TYPES})
endif ()

set (CMAKE_BUILD_TYPE DEBUG)
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    message("Buiding in debug mode")
endif (CMAKE_BUILD_TYPE MATCHES DEBUG)

list( APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
list( APPEND CMAKE_PREFIX_PATH "$ENV{HT}/cmake")

find_package( Houdini REQUIRED )
message (STATUS "Houdini HDK: ${Houdini_VERSION}")

find_package(Cpd REQUIRED)
message (STATUS "Cpd lib: ${Cpd_VERSION}")

find_package(Eigen3 REQUIRED)
message(STATUS "Eigen3: ${Eigen3_VERSION}")

if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
endif()

add_subdirectory (src/SOP_ShapeMatch)
add_subdirectory (src/SOP_Pca)





