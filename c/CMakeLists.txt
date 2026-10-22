list(APPEND SRCS
  models/ibm/IBForceCoupled.cpp
  models/ibm/IBMesh.cpp
  models/ibm/IBVelocityCoupled.cpp
  models/ModelRuntime.cpp
  models/ModelRuntimeCase.cpp
  models/ibm/IBCoupling.cpp
  models/rigidbody/IBRigidBody.cpp
)

list(APPEND HDRS 
  models/ibm/IBForceCoupled.h
  models/ibm/IBMesh.h
  models/ibm/IBVelocityCoupled.h
  models/ModelRuntime.h
  models/ibm/IBCoupling.h
  models/rigidbody/IBRigidBody.h
)

convert_filenames_to_full_paths(SRCS)
convert_filenames_to_full_paths(HDRS)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
