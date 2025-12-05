option(BUILD_BASILISK_INTERFACE "Build basilisk interface" ON)
if (BUILD_BASILISK_INTERFACE)
  add_subdirectory(interface)
endif()

include(ExternalProject)

ExternalProject_Add(basilisk
	GIT_REPOSITORY        https://github.com/tortotubus/beam-basilisk
  UPDATE_DISCONNECTED   1
	SOURCE_DIR            ${CMAKE_BINARY_DIR}/_ext/basilisk-src/
	BINARY_DIR            ${CMAKE_BINARY_DIR}/_ext/basilisk-bin/
	INSTALL_DIR           ${CMAKE_BINARY_DIR}/_ext/basilisk-install/

	CMAKE_ARGS
	-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
	-DCMAKE_BUILD_TYPE=Debug

	BUILD_COMMAND 
	${CMAKE_COMMAND}
	--build <BINARY_DIR>
	--config Debug
	--target qcc

	INSTALL_COMMAND
	${CMAKE_COMMAND}
	--install <BINARY_DIR>
	--prefix <INSTALL_DIR>
	--config ${CMAKE_BUILD_TYPE}
	--component QCC
)

ExternalProject_Get_Property(basilisk install_dir)
ExternalProject_Get_Property(basilisk binary_dir)

set(QCC_EXECUTABLE "${install_dir}/bin/qcc" CACHE FILEPATH "Path to built qcc")


#
#
#

option(BUILD_BASILISK_EXAMPLES "Build basilisk example programs" ON)
if(BUILD_BASILISK_EXAMPLES)
  add_subdirectory(examples/serial)
  add_subdirectory(examples/mpi)
endif()

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
