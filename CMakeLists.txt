execute_process(
        COMMAND "make"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

cmake_minimum_required(VERSION 3.6)
project(gecode-5.0.0-extension)

# Path to custom Gecode with cbs support
set(CMAKE_PREFIX_PATH /home/sam/gecode-5.0.0-cbs)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(Gecode REQUIRED)
include_directories(${Gecode_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

# This is only for having syntax coloring in Clion. Compiling cbs won't actually
# work.
add_executable(cbs cbs.hpp)
set_target_properties(cbs PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(cbs ${Gecode_LIBRARIES})


#set(EXAMPLES
#        "latinsquare"
#        "nonogram"
#        "langford-number"
#        "ortho-latin"
#        "pentominoes"
#        "qcp"
#        )

set(SQL_SRC "")
option(BUILD_SQL "Build SQLite profiling." OFF)

if(${BUILD_SQL})
    find_package(Threads REQUIRED)
    include_directories(sql)
    set(SQL_SRC sql/sqlite3.c sql/sql-interface.cpp)
    add_definitions(-DSQL)
endif()

include_directories(./)

foreach(ex "latinsquare" "nonogram" "qcp" "magic-square" "minimum-deficiency"
        "bin-packing" "crowded-chess")
    if(${ex} STREQUAL "bin-packing")
        find_package(Threads REQUIRED)
    endif()
    add_executable(${ex} models/${ex}.cpp ${SQL_SRC})
    target_link_libraries(${ex} ${Gecode_LIBRARIES})
    if(${BUILD_SQL} OR ${ex} MATCHES "bin-packing")
        target_link_libraries(${ex} Threads::Threads ${CMAKE_DL_LIBS})
    endif()
endforeach(ex)
