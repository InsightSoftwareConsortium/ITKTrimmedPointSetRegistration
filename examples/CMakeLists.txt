add_executable( TrimmedPointSetRegistrationExample TrimmedPointSetRegistrationExample.cxx )
 
target_link_libraries( TrimmedPointSetRegistrationExample ${ITK_LIBRARIES})

target_include_directories( TrimmedPointSetRegistrationExample PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../include" )

set_target_properties( TrimmedPointSetRegistrationExample  
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/examples"
)
