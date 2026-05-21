find_program(LATEXMK_EXECUTABLE NAMES latexmk REQUIRED)
find_program(BIBER_EXECUTABLE NAMES biber REQUIRED)

add_subdirectory(figure)

file(GLOB_RECURSE THESIS_SOURCES
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.tex"
  "${CMAKE_CURRENT_SOURCE_DIR}/admin/*.tex"
  "${CMAKE_CURRENT_SOURCE_DIR}/header/*.tex"
  "${CMAKE_CURRENT_SOURCE_DIR}/body/*.tex"
  "${CMAKE_CURRENT_SOURCE_DIR}/appendix/*.tex"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.bib"
  "${CMAKE_CURRENT_SOURCE_DIR}/figure/*.pdf"
  "${CMAKE_CURRENT_SOURCE_DIR}/figure/*.pdf_tex"
  "${CMAKE_CURRENT_SOURCE_DIR}/figure/*.svg"
)

list(APPEND THESIS_SOURCES
  ${THESIS_FIGURE_INPUTS}
  ${THESIS_GENERATED_FIGURES}
)

set(THESIS_PDF "${CMAKE_CURRENT_BINARY_DIR}/main.pdf")

add_custom_command(
  OUTPUT "${THESIS_PDF}"
  COMMAND "${LATEXMK_EXECUTABLE}"
    -xelatex
    -interaction=nonstopmode
    -file-line-error
    -synctex=1
    "-output-directory=${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/main.tex"
  DEPENDS ${THESIS_SOURCES}
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  COMMENT "Building thesis PDF"
  VERBATIM
)

add_custom_target(thesis ALL
  DEPENDS "${THESIS_PDF}"
)

add_dependencies(thesis thesis_figures)


# install(
#   FILES "${THESIS_PDF}"
#   DESTINATION "${CMAKE_INSTALL_DOCDIR}"
#   COMPONENT ELFF_Documentation
#   OPTIONAL
# )
