list(APPEND SRCS
  error.cpp
  text.cpp
  globals.cpp 
)

list(APPEND HDRS
  error.hpp
  text.hpp
  globals.hpp 
)

convert_filenames_to_full_paths(SRCS)
convert_filenames_to_full_paths(HDRS)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)