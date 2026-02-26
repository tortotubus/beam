set(HDRS
    CommHandle.hpp
    Policy.hpp
    Runtime.hpp
)

set(SRCS 
    CommHandle.cpp
    Runtime.cpp
)

convert_filenames_to_full_paths(HDRS)
convert_filenames_to_full_paths(SRCS)

set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)