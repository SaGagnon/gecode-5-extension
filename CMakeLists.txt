cmake_minimum_required(VERSION 3.7)
project(gecode-5.0.0-extension)


set(CMAKE_PREFIX_PATH /home/sam/gecode-5.0.0-cbs)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(Gecode REQUIRED)
include_directories(${Gecode_INCLUDE_DIRS})

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cbs)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp cbs/cbs.hpp)
add_executable(gecode-5.0.0-extension ${SOURCE_FILES})
target_link_libraries(gecode-5.0.0-extension ${Gecode_LIBRARIES})
