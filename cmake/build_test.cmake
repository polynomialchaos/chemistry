####################################################################################################################################
# Chemistry - Finite rate chemistry library and solver
# (c) 2020 | Florian Eigentler
####################################################################################################################################
set( test_target "test_${target}" )

####################################################################################################################################
# Add subdirecotries
#-----------------------------------------------------------------------------------------------------------------------------------
file( GLOB SOURCES
    test/*/*.c
    test/*.c
)

####################################################################################################################################
# External libraries
#-----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Build target
#-----------------------------------------------------------------------------------------------------------------------------------
add_executable( ${test_target} ${SOURCES} )
add_dependencies( ${test_target} ${PROJECT_NAME} )
target_link_libraries( ${test_target} ${PROJECT_NAME} )
set_target_properties( ${test_target} PROPERTIES LINKER_LANGUAGE C )
set_target_properties( ${test_target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/test" )

####################################################################################################################################
# Setup test environment
#-----------------------------------------------------------------------------------------------------------------------------------
add_test(
    NAME ${test_target}
    COMMAND ${test_target}
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/test"
)