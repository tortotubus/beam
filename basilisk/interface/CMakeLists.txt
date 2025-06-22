list(APPEND SRCS
  dummy.cpp
  immersedboundary.cpp
  spring.cpp
  #vertex.cpp
)

list(APPEND HDRS
  dummy.h
  immersedboundary.h
  spring.h
  #vertex.h
)

convert_filenames_to_full_paths(SRCS)
convert_filenames_to_full_paths(HDRS)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)