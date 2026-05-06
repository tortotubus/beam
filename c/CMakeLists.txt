list(APPEND SRCS
  models/beam/IBEulerBeam.cpp
  models/ibm/IBForceCoupled.cpp
  models/ibm/IBMesh.cpp
  models/ibm/IBVelocityCoupled.cpp
  models/ibm/IBRuntime.cpp
  models/rigidbody/IBRigidBody.cpp
)

list(APPEND HDRS 
  models/beam/IBEulerBeam.h
  models/ibm/IBForceCoupled.h
  models/ibm/IBMesh.h
  models/ibm/IBVelocityCoupled.h
  models/ibm/IBRuntime.h
  models/rigidbody/IBRigidBody.h
)

convert_filenames_to_full_paths(SRCS)
convert_filenames_to_full_paths(HDRS)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
