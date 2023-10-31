set(AMGCL_ROOT ${OPEN321_DOWNLOADS_ROOT}/amgcl)
ExternalProject_Add(
        amgcl
        SOURCE_DIR ${AMGCL_ROOT}
        BINARY_DIR "${CMAKE_BINARY_DIR}/thirdparty/amgcl"
        GIT_REPOSITORY  "https://github.com/ddemidov/amgcl"
        GIT_TAG master
        UPDATE_COMMAND ""
        CMAKE_ARGS
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
        -DCMAKE_INSTALL_PREFIX:PATH=${OPEN321_DEPENDENCIES_INSTALL_DIR}/amgcl
        INSTALL_DIR ${OPEN321_DEPENDENCIES_INSTALL_DIR}/amgcl
)
set(AMGCL_INCLUDE_DIRS "${OPEN321_DEPENDENCIES_INSTALL_DIR}/amgcl/include/")
set(AMGCL_LIBRARY_DIRS "")
include_directories(${AMGCL_INCLUDE_DIRS})
add_definitions(-DAMGCL_NO_BOOST)
