function(elff_basilisk_add_executable SOURCE_FILE)
  get_filename_component(source_name ${SOURCE_FILE} NAME_WE)
  set(output_c "${CMAKE_CURRENT_BINARY_DIR}/_${source_name}.c")
  
  file(GLOB_RECURSE basilisk_headers
    CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/basilisk/*.h"
  )

  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/_${source_name}.c"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${SOURCE_FILE}" "${CMAKE_CURRENT_BINARY_DIR}/${source_name}.c"
    COMMAND $<TARGET_FILE:basilisk::qcc>
      "${source_name}.c"
      -I"${CMAKE_SOURCE_DIR}"
      -I"${CMAKE_SOURCE_DIR}/basilisk" 
      -DTRACE=2
      -source
    DEPENDS ${SOURCE_FILE} ${basilisk_headers}
    BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/_${source_name}.c"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" 
  )
  
  add_executable(${source_name} "_${source_name}.c")

  target_link_libraries(${source_name}
    PRIVATE
      HDF5::HDF5
      m
  )

  set_target_properties(${source_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR}"
  )
endfunction()

