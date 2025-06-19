list(APPEND SRCS
  basilisk_interface.cpp
)

list(APPEND HDRS
  basilisk_interface.h
)

convert_filenames_to_full_paths(SRCS)
convert_filenames_to_full_paths(HDRS)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)