set(SRCS
  bilinearform.cpp
  coefficient.cpp
  fe_base.cpp
  fe_coll.cpp
  fespace.cpp
  geom.cpp
)

set(HDRS
  bilinearform.hpp
  coefficient.hpp
  fe_base.hpp
  fe_coll.hpp
  fespace.hpp
  geom.hpp
)

convert_filenames_to_full_paths(SRCS)
convert_filenames_to_full_paths(HDRS)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)