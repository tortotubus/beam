set(HDRS
 config.hpp
)

convert_filenames_to_full_paths(HDRS)

set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)