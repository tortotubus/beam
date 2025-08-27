set(HDRS
  immersedboundary.hpp
  spring/Spring.hpp
  beam/EulerBeamStatic.hpp
  beam/EulerBeamStaticInextensible.hpp
  beam/EulerBeamStaticInextensibleADDM.hpp
  beam/EulerBeamStaticInextensibleADDMSparse.hpp
  beam/Shapes.hpp
)

set(SRCS
  immersedboundary.cpp
  spring/Spring.cpp
  beam/EulerBeamStatic.cpp
  beam/EulerBeamStaticInextensible.cpp
  beam/EulerBeamStaticInextensibleADDM.cpp
  beam/EulerBeamStaticInextensibleADDMSparse.cpp
)

convert_filenames_to_full_paths(HDRS)
convert_filenames_to_full_paths(SRCS)

set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)