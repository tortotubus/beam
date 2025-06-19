list(APPEND SRCS
  mesh.cpp
  point.cpp
  segment.cpp
  triangle.cpp
  vertex.cpp
)

list(APPEND HDRS
  mesh.hpp
  point.hpp
  segment.hpp
  triangle.hpp
  vertex.hpp
)

convert_filenames_to_full_paths(SRCS)
convert_filenames_to_full_paths(HDRS)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)