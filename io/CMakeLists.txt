list(APPEND SRCS
  C/vtkHDF.cpp
  C/vtkHDFPolyData.cpp
  C/vtkPolyData.cpp
  CXX/vtkHDFPolyData.cpp
  CXX/vtkPolyData.cpp
)

list(APPEND HDRS
  C/vtkHDF.hpp
  C/vtkHDFPolyData.hpp
  C/vtkPolyData.hpp
  CXX/vtkHDFPolyData.hpp
  CXX/vtkPolyData.hpp
)

convert_filenames_to_full_paths(SRCS)
convert_filenames_to_full_paths(HDRS)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)