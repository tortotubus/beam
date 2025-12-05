list(APPEND SRCS 
  ibm/immersedboundary.cpp
  ibm/spring.cpp
  ibm/beam.cpp
  # df/directforcing.cpp
)

list(APPEND HDRS 
  ibm/immersedboundary.h
  ibm/spring.h
  ibm/beam.h
  # df/directforcing.h
)

convert_filenames_to_full_paths(SRCS)
convert_filenames_to_full_paths(HDRS)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)