set(HDRS
  ibm/immersedboundary.hpp
  spring/Spring.hpp
  beam/IBEulerBeam.hpp
  beam/EulerBeamDynamicInextensibleADDM.hpp
  beam/EulerBeamDynamicInextensibleAugKKT.hpp
  beam/EulerBeamDynamicInextensibleMoM.hpp
  beam/EulerBeamDynamicInextensibleMoMSparse.hpp
  beam/EulerBeamStatic.hpp
  beam/EulerBeamStaticInextensibleADDM.hpp
  beam/EulerBeamStaticInextensibleAugKKT.hpp
  beam/EulerBeamStaticInextensibleMoM.hpp
  beam/EulerBeamStaticInextensibleMoMSparse.hpp
  beam/Shapes.hpp
)

set(SRCS
  ibm/immersedboundary.cpp
  spring/Spring.cpp
  beam/IBEulerBeam.hpp
  beam/EulerBeamDynamicInextensibleADDM.cpp
  beam/EulerBeamDynamicInextensibleAugKKT.cpp
  beam/EulerBeamDynamicInextensibleMoM.cpp
  beam/EulerBeamDynamicInextensibleMoMSparse.cpp
  beam/EulerBeamStatic.cpp
  beam/EulerBeamStaticInextensibleADDM.cpp
  beam/EulerBeamStaticInextensibleAugKKT.cpp
  beam/EulerBeamStaticInextensibleMoM.cpp
  beam/EulerBeamStaticInextensibleMoMSparse.cpp
)

convert_filenames_to_full_paths(HDRS)
convert_filenames_to_full_paths(SRCS)

set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)