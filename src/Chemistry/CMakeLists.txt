####################################################################################################################################
# Chemistry - Finite rate chemistry library and solver
# (c) 2020 | Florian Eigentler
####################################################################################################################################
set( target "${PROJECT_NAME}" )

####################################################################################################################################
# Add subdirecotries
#-----------------------------------------------------------------------------------------------------------------------------------
file( GLOB SOURCES
    src/*/*.c
    src/*.c
)

####################################################################################################################################
# External libraries
#-----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Build target
#-----------------------------------------------------------------------------------------------------------------------------------
if( SHARED )
    add_library( ${target} SHARED ${SOURCES} )
else()
    add_library( ${target} STATIC ${SOURCES} )
endif()
target_include_directories( ${target} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src" )
target_link_libraries( ${target} C ${external_libraries} )
set_target_properties( ${target} PROPERTIES LINKER_LANGUAGE C )
set_target_properties( ${target} PROPERTIES COMPILE_FLAGS "${external_library_compile_flags}" )
set_target_properties( ${target} PROPERTIES LINK_FLAGS "${external_library_link_flags}" )

####################################################################################################################################
# Install target
#-----------------------------------------------------------------------------------------------------------------------------------
install( DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h" )
install( TARGETS ${target} EXPORT Chemistry
    RUNTIME DESTINATION bin INCLUDES DESTINATION include LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )

####################################################################################################################################
# Build test target
#-----------------------------------------------------------------------------------------------------------------------------------
if( BUILD_TESTING )
    include( build_test )
endif()