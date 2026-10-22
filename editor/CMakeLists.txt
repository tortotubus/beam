cmake_minimum_required(VERSION 3.12)
project(structured_xml_editor LANGUAGES CXX)
include(GNUInstallDirs)
if(NOT ELFF_INSTALL_BINDIR)
  set(ELFF_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}")
endif()
find_package(LibXml2 REQUIRED)
set(CURSES_NEED_NCURSES TRUE)
find_package(Curses REQUIRED)
# Keep the editor's single schema in sync with the repository at build time.
set(EDITOR_DTD "${CMAKE_CURRENT_SOURCE_DIR}/../io/DTD/ELFF/ELFF-v2.dtd")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${EDITOR_DTD}")
file(READ "${EDITOR_DTD}" EDITOR_DTD_TEXT)
configure_file(schema.hpp.in "${CMAKE_CURRENT_BINARY_DIR}/generated/schema.hpp" @ONLY)
add_executable(elff-xml main.cpp document.cpp constraints.cpp terminal.cpp editor.cpp projection.cpp cursor.cpp)
target_compile_features(elff-xml PRIVATE cxx_std_17)
target_include_directories(elff-xml PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated" ${CURSES_INCLUDE_DIRS})
set_target_properties(elff-xml PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(elff-xml PRIVATE LibXml2::LibXml2 ${CURSES_LIBRARIES})
install(TARGETS elff-xml
  RUNTIME DESTINATION "${ELFF_INSTALL_BINDIR}"
  COMPONENT ELFF_Runtime
)
target_compile_options(elff-xml PRIVATE -Wall -Wextra -Wpedantic)
# Direct DTD model inspection is intentional; libxml2 2.15 marks these
# structures deprecated without providing an equivalent model inspection API.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(elff-xml PRIVATE -Wno-deprecated-declarations)
endif()
