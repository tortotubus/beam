find_package(Doxygen REQUIRED)

set(ELFF_DOXYGEN_DOT_PATH "")
if(DOXYGEN_DOT_FOUND)
  get_filename_component(ELFF_DOXYGEN_DOT_PATH "${DOXYGEN_DOT_EXECUTABLE}" DIRECTORY)
endif()

# Keep the checked-in styling while writing generated files into the build tree.
set(ELFF_DOXYGEN_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
file(GENERATE OUTPUT "${ELFF_DOXYGEN_CONFIG}" CONTENT
"@INCLUDE = \"${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile\"
OUTPUT_DIRECTORY = \"${CMAKE_CURRENT_BINARY_DIR}\"
HAVE_DOT = $<IF:$<BOOL:${DOXYGEN_DOT_FOUND}>,YES,NO>
DOT_PATH = \"${ELFF_DOXYGEN_DOT_PATH}\"
")

add_custom_target(doxygen ALL
  COMMAND Doxygen::doxygen "${ELFF_DOXYGEN_CONFIG}"
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  COMMENT "Generating ELFF API documentation"
  VERBATIM
)
