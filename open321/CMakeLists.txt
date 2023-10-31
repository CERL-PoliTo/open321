function(Open321_declareLibrary libName)
        add_subdirectory(src/)
        set_property(TARGET ${libName} PROPERTY FOLDER "open321")
endfunction(Open321_declareLibrary)

# Project sub-folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#Configure the test header to provide usefull values
set(open321_CONFIG_INC_IN "${CMAKE_CURRENT_SOURCE_DIR}/include/open321_config.hpp.in")
set(open321_CONFIG_INC_OUT "${CMAKE_CURRENT_BINARY_DIR}/include/open321_config.hpp")
configure_file(${open321_CONFIG_INC_IN} ${open321_CONFIG_INC_OUT} @ONLY)

Open321_declareLibrary(open321Preconditioner)
