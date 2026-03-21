set(HDRS
  beam/EulerBeam.hpp
  beam/EulerBeamDynamicInextensibleADDM.hpp
  beam/EulerBeamDynamicInextensibleAugKKT.hpp
  beam/EulerBeamDynamicInextensibleMoM.hpp
  beam/EulerBeamDynamicInextensibleMoMSparse.hpp
  beam/EulerBeamMesh.hpp
  beam/EulerBeamStatic.hpp
  beam/EulerBeamStaticInextensibleADDM.hpp
  beam/EulerBeamStaticInextensibleAugKKT.hpp
  beam/EulerBeamStaticInextensibleMoM.hpp
  beam/EulerBeamStaticInextensibleMoMSparse.hpp
  beam/IBEulerBeam.hpp
  beam/Shapes.hpp
  ibm/IBForceCoupled.hpp
  ibm/IBMesh.hpp
  ibm/IBRuntime.hpp
  ibm/IBVelocityCoupled.hpp
  spring/IBSpringCircle.hpp
)

set(SRCS
  beam/EulerBeam.cpp
  beam/EulerBeamDynamicInextensibleADDM.cpp
  beam/EulerBeamDynamicInextensibleAugKKT.cpp
  beam/EulerBeamDynamicInextensibleMoM.cpp
  beam/EulerBeamDynamicInextensibleMoMSparse.cpp
  beam/EulerBeamMesh.cpp
  beam/EulerBeamStatic.cpp
  beam/EulerBeamStaticInextensibleADDM.cpp
  beam/EulerBeamStaticInextensibleAugKKT.cpp
  beam/EulerBeamStaticInextensibleMoM.cpp
  beam/EulerBeamStaticInextensibleMoMSparse.cpp
  beam/IBEulerBeam.cpp
  # beam/Shapes.cpp
  ibm/IBForceCoupled.cpp
  ibm/IBMesh.cpp
  ibm/IBVelocityCoupled.cpp
  spring/IBSpringCircle.cpp
)

convert_filenames_to_full_paths(HDRS)
convert_filenames_to_full_paths(SRCS)

set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)