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

#find_package(Fgt REQUIRED)
#message (STATUS "Fgt lib: ${Houdini_VERSION}")

find_package(Cpd REQUIRED)
message (STATUS "Cpd lib: ${Houdini_VERSION}")

find_package(Eigen3 REQUIRED)
message(STATUS "Eigen3 version: ${EIGEN3_VERSION}")

# Main sources
list(APPEND _sources
        src/SOP_ShapeMatch.cpp
        src/SOP_ShapeMatch.h
)

add_library( SOP_ShapeMatch SHARED
        ${_sources}
        )

target_include_directories(SOP_ShapeMatch
    PUBLIC ${EIGEN3_INCLUDE_DIR})

target_link_libraries(SOP_ShapeMatch
        PUBLIC
        Cpd::Library-C++
        )

if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

endif()
target_link_libraries(SOP_ShapeMatch PUBLIC OpenMP::OpenMP_CXX)


target_link_libraries( SOP_ShapeMatch  PUBLIC Houdini )
houdini_configure_target( SOP_ShapeMatch )



