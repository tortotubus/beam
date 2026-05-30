cmake_minimum_required(VERSION 3.12.0...4.0.0)

project(ELFFBeamerDocumentation NONE)

find_program(LATEXMK_EXECUTABLE NAMES latexmk REQUIRED)
find_program(BIBER_EXECUTABLE NAMES biber REQUIRED)

if(EXISTS "${CMAKE_SOURCE_DIR}/third_party/ubcbeamer/ubcbeamer_submodule")
  set(BEAMER_REPOSITORY_SOURCE_DIR "${CMAKE_SOURCE_DIR}")
else()
  get_filename_component(BEAMER_REPOSITORY_SOURCE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../.."
    ABSOLUTE
  )
endif()

add_subdirectory(figure)

file(GLOB_RECURSE BEAMER_SOURCES
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.tex"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.bib"
  "${CMAKE_CURRENT_SOURCE_DIR}/figure/vector/*.pdf"
  "${CMAKE_CURRENT_SOURCE_DIR}/figure/vector/*.pdf_tex"
  "${CMAKE_CURRENT_SOURCE_DIR}/figure/vector/*.tex"
  "${CMAKE_CURRENT_SOURCE_DIR}/figure/raster/*.png"
  "${BEAMER_REPOSITORY_SOURCE_DIR}/third_party/ubcbeamer/ubcbeamer_submodule/tex/latex/ubcbeamer/*.sty"
  "${BEAMER_REPOSITORY_SOURCE_DIR}/third_party/ubcbeamer/ubcbeamer_submodule/tex/latex/ubcbeamer/UBC-images/*.pdf"
)

list(APPEND BEAMER_SOURCES
  ${BEAMER_FIGURE_INPUTS}
  ${BEAMER_GENERATED_FIGURES}
)

set(BEAMER_PDF "${CMAKE_CURRENT_BINARY_DIR}/main.pdf")

add_custom_command(
  OUTPUT "${BEAMER_PDF}"
  COMMAND "${LATEXMK_EXECUTABLE}"
    -xelatex
    -interaction=nonstopmode
    -file-line-error
    -synctex=1
    "-output-directory=${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/main.tex"
  DEPENDS ${BEAMER_SOURCES}
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  COMMENT "Building beamer PDF"
  VERBATIM
)

add_custom_target(beamer ALL
  DEPENDS "${BEAMER_PDF}"
)

add_dependencies(beamer beamer_figures)
