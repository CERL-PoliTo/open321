file(GLOB_RECURSE qh_preconditioning_example_files "./*")
add_library(qh_preconditioning_example SHARED "${qh_preconditioning_example_files}")
target_include_directories(qh_preconditioning_example PRIVATE "${Open321LibsIncludeDir}")

# Data
set(Open321_EXAMPLE_INC_IN ${CMAKE_CURRENT_SOURCE_DIR}/qh_example_config.hpp.in)
set(Open321_EXAMPLE_INC_OUT ${CMAKE_CURRENT_BINARY_DIR}/qh_example_config.hpp)
configure_file(${Open321_EXAMPLE_INC_IN} ${Open321_EXAMPLE_INC_OUT} @ONLY)

target_include_directories(qh_preconditioning_example PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

add_executable(main_simple_exe main.cpp ${qh_preconditioning_example_hpp_files})
target_include_directories(main_simple_exe PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(main_simple_exe qh_preconditioning_example open321Preconditioner)
