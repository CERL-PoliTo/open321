function(Open321_declareTest testName exeTargetName)
        add_subdirectory(${exeTargetName})
        set_property(TARGET ${exeTargetName} PROPERTY FOLDER "Tests")
        add_test(NAME ${testName} COMMAND ${exeTargetName})
endfunction(Open321_declareTest)

# Project sub-folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

Open321_declareTest(open321_qh_validation open321_qh_tests)
