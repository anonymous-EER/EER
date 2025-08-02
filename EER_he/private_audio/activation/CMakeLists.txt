project(HEaaNTest LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(HEaaN REQUIRED)
# 또는
find_library(HEAAN_LIB NAMES HEaaN heaan PATHS /usr/local/lib)

add_executable(activation activation.cpp activation.hpp)
target_link_libraries(activation HEaaN HEaaN-math /usr/local/include)