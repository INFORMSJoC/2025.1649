# Configuration
cmake_minimum_required(VERSION 3.26)
project(CircleCovering)

# Compiler 

# Files
file(GLOB_RECURSE SRCS src/*.cpp)
file(GLOB_RECURSE HEADERS include/*)

# Search directories
include_directories(include/)
include_directories(include/clipper)
include_directories(include/cg_descent)


# Executable
set(EXECUTABLE_NAME "CircleCovering")
add_executable(${EXECUTABLE_NAME} ${SRCS} ${HEADERS})