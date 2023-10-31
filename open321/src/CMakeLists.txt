file(GLOB_RECURSE open321Preconditioner_files
    "${Open321LibsIncludeDir}/*"
    "./*")

add_library(open321Preconditioner ${open321Preconditioner_files})
target_include_directories(open321Preconditioner PUBLIC "${Open321LibsIncludeDir}")

if(TARGET Eigen3::Eigen)
   target_link_libraries(open321Preconditioner PUBLIC Eigen3::Eigen)
else()
    add_dependencies(open321Preconditioner Eigen3)
    target_include_directories(open321Preconditioner PUBLIC ${EIGEN_INCLUDE_DIRS})
    target_link_libraries(open321Preconditioner PUBLIC ${EIGEN_LIBRARIES})
endif()

add_dependencies(open321Preconditioner amgcl)
target_include_directories(open321Preconditioner PUBLIC ${AMGCL_INCLUDE_DIRS})