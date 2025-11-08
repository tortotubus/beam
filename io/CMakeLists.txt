list(APPEND SRCS
  VTKPolyData.cpp
  HDF.cpp
  VTKHDF.cpp
  VTKHDFPolyData.cpp
)

list(APPEND HDRS
  VTKPolyData.hpp
  HDF.hpp
  VTKHDF.hpp
  VTKHDFPolyData.hpp
)

convert_filenames_to_full_paths(SRCS)
convert_filenames_to_full_paths(HDRS)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)