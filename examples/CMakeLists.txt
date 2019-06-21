add_executable( TrimmedPointSetRegistrationExample TrimmedPointSetRegistrationExample.cxx )
 
target_link_libraries( TrimmedPointSetRegistrationExample ${ITK_LIBRARIES})

set_target_properties( TrimmedPointSetRegistrationExample  
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/examples"
)
