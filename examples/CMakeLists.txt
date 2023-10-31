function(Open321_declareExample exampleName)
        add_subdirectory(${exampleName})
        set_property(TARGET ${exampleName} PROPERTY FOLDER "Examples")
endfunction(Open321_declareExample)

# Project sub-folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

Open321_declareExample(qh_preconditioning_example)

