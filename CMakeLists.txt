cmake_minimum_required(VERSION 3.18)

enable_language(CXX)
enable_testing()

# Options
option(BUILD_SHARED_LIBS "Build a shared library(ON) or a static library (OFF)" ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_STANDARD 17)

project(open321 CXX)

# Include directory
set(Open321LibsIncludeDir "${CMAKE_SOURCE_DIR}/open321/include")
include_directories("${CMAKE_BINARY_DIR}/open321/include")

# Output
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Data
set(OPEN321_TEST_DATA_ROOT "${CMAKE_SOURCE_DIR}/Data/test")
set(OPEN321_EXAMPLE_DATA_ROOT "${CMAKE_SOURCE_DIR}/Data/example")

# Thirdparty
set(OPEN321_DOWNLOADS_ROOT "${CMAKE_SOURCE_DIR}/downloads")
file(MAKE_DIRECTORY ${OPEN321_DOWNLOADS_ROOT})
set(OPEN321_DEPENDENCIES_INSTALL_DIR "${CMAKE_BINARY_DIR}/thirdparty_install")

include(ExternalProject)
include(thirdparty/eigen.cmake)
include(thirdparty/amgcl.cmake)
include(thirdparty/googletest.cmake)
include_directories(OPEN321_DOWNLOADS_ROOT)

# Compiler/linker warnings
if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # Can't use /Wall for MSVC because their headers generate multiple warnings...
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /W4 /std:c++17")
    # ignore the std::tr1 deprecation in gtest
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4996")
    # ignore the "conditional is a constant" that triggers false positive on if that are actually if constexpr
    # The day we move to C++17 we should remove this ignore warning and change the corresponding if to if constexpr
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4127")
    # ignore <complex> header that triggers a "possible division by 0" warning
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4723")
    # ignore the "multiple copy constructors specified" warning
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4521")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4715")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4459")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd /Zi")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD /O2 /Ot /GL /GS- /GR-")
    # The set of options bellow are here for reference: They should be faster but might
    # give incorrect results.
    #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD /O2 /Ot /GL /GS- /Qpar /fp:fast /GR-")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS}")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND WIN32)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_DEBUG} /MD /MP /Qstd=c++17")
     set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD /Ox /MP /Qstd=c++17")
     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
     set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")
     set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS}")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++17 -Wall -Wextra -Werror -pedantic -fdiagnostics-color=auto")
     set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3")
     set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")
     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
     if(PROFILING)
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3")
         set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lprofiler -ltcmalloc")
     endif()
endif()

# Libraries
add_subdirectory(${CMAKE_SOURCE_DIR}/open321 "open321")
add_subdirectory(${CMAKE_SOURCE_DIR}/examples "examples")
add_subdirectory(${CMAKE_SOURCE_DIR}/tests "tests")
