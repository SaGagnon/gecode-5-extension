cmake_minimum_required(VERSION 3.7)
project(gecode-5.0.0-extension)

# Path to custom Gecode with cbs support
set(CMAKE_PREFIX_PATH /home/sam/gecode-5.0.0-cbs)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(Gecode REQUIRED)
include_directories(${Gecode_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 11)

# This is only for having syntax coloring in Clion. Compiling cbs won't actually
# work.
add_executable(cbs cbs.hpp)
set_target_properties(cbs PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(cbs ${Gecode_LIBRARIES})

add_executable(latinsquare latinsquare.cpp)
target_link_libraries(latinsquare ${Gecode_LIBRARIES})

add_executable(nonogram nonogram.cpp)
target_link_libraries(nonogram ${Gecode_LIBRARIES})

add_executable(langford-number langford-number.cpp)
target_link_libraries(langford-number ${Gecode_LIBRARIES})
