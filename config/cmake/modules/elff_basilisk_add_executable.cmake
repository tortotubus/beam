function(elff_basilisk_add_executable SOURCE_FILE)
  get_filename_component(source_name ${SOURCE_FILE} NAME_WE)
  set(output_c "${CMAKE_CURRENT_BINARY_DIR}/_${source_name}.c")
  set(qcc_options ${ARGN})
  
  file(GLOB_RECURSE basilisk_headers
    CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/basilisk/*.h"
    "${CMAKE_SOURCE_DIR}/basilisk/templates/*.c"
  )

  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/_${source_name}.c"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${SOURCE_FILE}" "${CMAKE_CURRENT_BINARY_DIR}/${source_name}.c"
    COMMAND $<TARGET_FILE:basilisk::qcc>
      ${qcc_options}
      -DTRACE=3
      "${source_name}.c"
      -I"${CMAKE_SOURCE_DIR}/basilisk" 
      -I"${CMAKE_BINARY_DIR}/include"
      -source
    DEPENDS ${SOURCE_FILE} ${basilisk_headers}
    BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/_${source_name}.c"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" 
  )
  
  add_executable(${source_name} "_${source_name}.c")

  if(ELFF_USE_MPI) 
    target_link_libraries(${source_name}
      PUBLIC
        MPI::MPI_C
    )    
  endif()

  if(ELFF_USE_HDF5) 
    target_link_libraries(${source_name}
      PUBLIC
        hdf5::hdf5
        hdf5::hdf5_hl  
    )    
  endif()

  target_link_libraries(${source_name}
    PUBLIC
      ELFF
      m
  )

  set_target_properties(${source_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR}"
    BUILD_RPATH "${CMAKE_BINARY_DIR}"
  )

  install(TARGETS ${source_name}
    RUNTIME DESTINATION ${ELFF_INSTALL_BINDIR}
    COMPONENT ElFF_Runtime
  )
endfunction()
