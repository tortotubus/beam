set(HDRS
  immersedboundary.hpp
  spring/Spring.hpp
  beam/EulerBeamDynamicInextensibleADDM.hpp
  beam/EulerBeamDynamicInextensibleMoM.hpp
  beam/EulerBeamStatic.hpp
  beam/EulerBeamStaticInextensibleADDM.hpp
  beam/EulerBeamStaticInextensibleMoM.hpp
  beam/EulerBeamStaticInextensibleMoMSparse.hpp
  beam/Shapes.hpp
)

set(SRCS
  immersedboundary.cpp
  spring/Spring.cpp
  beam/EulerBeamStatic.cpp
  beam/EulerBeamStaticInextensibleADDM.cpp
  beam/EulerBeamStaticInextensibleMoM.cpp
  beam/EulerBeamStaticInextensibleMoMSparse.cpp
)

convert_filenames_to_full_paths(HDRS)
convert_filenames_to_full_paths(SRCS)

set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)