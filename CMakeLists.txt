cmake_minimum_required(VERSION 3.17)
project(roko VERSION 0.1.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")

find_package(PkgConfig REQUIRED)
find_package(Python REQUIRED 
  COMPONENTS 
    Interpreter 
    Development 
    NumPy
)

# TODO: try finding htslib localy and then make a git pull
list(APPEND CMAKE_PREFIX_PATH "${PROJECT_BINARY_DIR}/htslib")
pkg_check_modules(HTSLIB REQUIRED htslib)

add_custom_target(
  htslib_dep
  COMMAND autoreconf -i
  COMMAND ./configure
  COMMAND make
  COMMAND make install prefix=${PROJECT_BINARY_DIR}/htslib
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib
)


set(${PROJECT_NAME}_SOURCES
  ${PROJECT_SOURCE_DIR}/src/gen.cpp
  ${PROJECT_SOURCE_DIR}/src/types.cpp
  ${PROJECT_SOURCE_DIR}/src/models.cpp
  ${PROJECT_SOURCE_DIR}/src/generate.cpp
)

Python_add_library(${PROJECT_NAME} MODULE ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} htslib_dep)


target_include_directories(${PROJECT_NAME} 
  PRIVATE #TODO: Generator expressions?
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${HTSLIB_INCLUDE_DIRS}
    # ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib/htslib
)

target_link_libraries(${PROJECT_NAME}
  PRIVATE
    Python::NumPy
    ${HTSLIB_LIBRARIES}
)
