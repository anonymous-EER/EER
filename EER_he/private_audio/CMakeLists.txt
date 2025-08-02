project(HEaaNTest LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

find_package(HEaaN REQUIRED)
# 또는
find_library(HEAAN_LIB NAMES HEaaN heaan PATHS /usr/local/lib)

add_executable(inference inference.cpp json.hpp model.hpp /home/user/EER/private_audio/activation /home/user/EER/private_audio/ccmm)
target_link_libraries(inference HEaaN HEaaN-math /usr/local/include )

add_executable(training training.cpp json.hpp model.hpp  /home/user/EER/private_audio/activation /home/user/EER/private_audio/ccmm)
target_link_libraries(training HEaaN HEaaN-math /usr/local/include )

