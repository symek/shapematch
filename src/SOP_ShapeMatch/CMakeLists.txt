add_library( SOP_ShapeMatch SHARED
        SOP_ShapeMatch.cpp
        SOP_ShapeMatch.h
        )

target_include_directories(SOP_ShapeMatch
        PUBLIC ${EIGEN3_INCLUDE_DIR})

target_link_libraries(SOP_ShapeMatch
        PUBLIC
        Cpd::Library-C++
        )

target_link_libraries(SOP_ShapeMatch PUBLIC OpenMP::OpenMP_CXX)

target_link_libraries( SOP_ShapeMatch PUBLIC Houdini )
houdini_configure_target( SOP_ShapeMatch )