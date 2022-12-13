cmake_minimum_required (VERSION 3.0)
# Creating the project
project(SCClient LANGUAGES Fortran)

# Set the source code directory

set(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(BUILD_TYPE SHARED)
set(LIB static)

message("Configuring SCClient source files: " ${SRC_DIR})
add_subdirectory(${SRC_DIR})

