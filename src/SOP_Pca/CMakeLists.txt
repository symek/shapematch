add_library(SOP_Pca SHARED
        SOP_Pca.cpp
        SOP_Pca.hpp
        pca.cpp
        pca.hpp)

target_include_directories(SOP_Pca
        PUBLIC ${EIGEN3_INCLUDE_DIR})


target_link_libraries( SOP_Pca PUBLIC Houdini )
houdini_configure_target( SOP_Pca )
