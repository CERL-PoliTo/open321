file(GLOB_RECURSE qh_preconditioning_test_files "./*")

# Data
set(Open321_TEST_INC_IN ${CMAKE_CURRENT_SOURCE_DIR}/qh_test_config.hpp.in)
set(Open321_TEST_INC_OUT ${CMAKE_CURRENT_BINARY_DIR}/qh_test_config.hpp)
configure_file(${Open321_TEST_INC_IN} ${Open321_TEST_INC_OUT} @ONLY)

add_executable(open321_qh_tests ${qh_preconditioning_test_files})
target_include_directories(open321_qh_tests PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(open321_qh_tests PRIVATE open321Preconditioner)

if(NOT ${GTEST_FOUND})
    add_dependencies(open321_qh_tests googletest)
    target_include_directories(open321_qh_tests PRIVATE "${GTEST_INCLUDE_DIRS}")
    target_link_directories(open321_qh_tests PRIVATE "${GTEST_LIBRARY_DIRS}")
    target_link_libraries(open321_qh_tests PRIVATE "${GTEST_LIBRARIES}")
else()
    target_link_libraries(open321_qh_tests PRIVATE GTest::GTest GTest::Main)
endif()
