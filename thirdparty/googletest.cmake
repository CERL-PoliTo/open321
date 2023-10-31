find_package(GTest)

if(NOT ${GTEST_FOUND})
    set(GTEST_ROOT ${OPEN321_DOWNLOADS_ROOT}/googletest)
    ExternalProject_Add(
        googletest
        SOURCE_DIR ${GTEST_ROOT}
        BINARY_DIR "${CMAKE_BINARY_DIR}/thirdparty/googletest"
        GIT_REPOSITORY "https://github.com/google/googletest.git"
        GIT_TAG "release-1.10.0"
        GIT_SHALLOW TRUE
        UPDATE_COMMAND ""
        CMAKE_ARGS
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
	    -DCMAKE_INSTALL_PREFIX:PATH=${OPEN321_DEPENDENCIES_INSTALL_DIR}/googletest
            -DBUILD_GMOCK:BOOL=OFF
            -DINSTALL_GTEST:BOOL=ON
            -Dgtest_force_shared_crt:BOOL=ON # for windows use C runtime DLL even in static libraries
	INSTALL_DIR "${OPEN321_DEPENDENCIES_INSTALL_DIR}/googletest"
        )
    set(GTEST_INCLUDE_DIRS "${OPEN321_DEPENDENCIES_INSTALL_DIR}/googletest/include")
    set(GTEST_LIBRARY_DIRS "${OPEN321_DEPENDENCIES_INSTALL_DIR}/googletest/lib"
                           "${OPEN321_DEPENDENCIES_INSTALL_DIR}/googletest/lib64"
                           "${OPEN321_DEPENDENCIES_INSTALL_DIR}/googletest/lib/$<CONFIG>")

    set(GTEST_LIBRARIES "gtest$<$<CONFIG:Debug>:d>" "gtest_main$<$<CONFIG:Debug>:d>")
endif()

