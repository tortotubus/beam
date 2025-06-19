set(HDRS
  immersedboundary.hpp
  spring/Spring.hpp
  #beam/EulerBeam.hpp
  #beam/InextensibleEulerBeam.hpp
  #beam/StaticEulerBeam.hpp
)

set(SRCS
  immersedboundary.cpp
  spring/Spring.cpp
  #beam/EulerBeam.cpp
  #beam/InextensibleEulerBeam.cpp
  #beam/StaticEulerBeam.cpp
)

convert_filenames_to_full_paths(HDRS)
convert_filenames_to_full_paths(SRCS)

set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)