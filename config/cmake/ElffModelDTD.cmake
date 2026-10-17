# Expand the root DTD's native includes at configure time. Runtime validation
# consumes only the flattened, embedded schema and needs no entity resolver.
set_property(GLOBAL PROPERTY ELFF_DTD_BUNDLER
  "${CMAKE_CURRENT_LIST_DIR}/../scripts/bundle_dtd.py")

function(elff_bundle_dtd root output)
  find_package(Python3 3.7 REQUIRED COMPONENTS Interpreter)
  get_property(bundler GLOBAL PROPERTY ELFF_DTD_BUNDLER)
  execute_process(
    COMMAND "${Python3_EXECUTABLE}" "${bundler}" "${root}" "${output}"
    RESULT_VARIABLE status OUTPUT_VARIABLE dependencies ERROR_VARIABLE diagnostic
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "DTD bundling failed:\n${diagnostic}")
  endif()
  string(REPLACE ";" "\\;" dependencies "${dependencies}")
  string(REPLACE "\n" ";" dependencies "${dependencies}")
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${bundler}" ${dependencies})
endfunction()
