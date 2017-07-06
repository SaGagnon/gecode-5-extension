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

# This is only for having syntax coloring in Clion. Compiling cbs won't actually
# work.
add_executable(cbs cbs.hpp)
set_target_properties(cbs PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(cbs ${Gecode_LIBRARIES})


set(EXAMPLES
        "latinsquare"
        "nonogram"
        "langford-number"
        "ortho-latin"
        "pentominoes"
        "qcp"
        )

#add_custom_target(
#        sql_structs ALL
#        COMMAND "{CMAKE_CURRENT_SOURCE_DIR}/sql/sql-interface.hh.sh"
#        DEPENDS {CMAKE_CURRENT_SOURCE_DIR}/schema.txt
#)

set(SQL_SRC "")
option(BUILD_SQL "Build SQLite profiling." OFF)
if(${BUILD_SQL})
    find_package(Threads REQUIRED)
    include_directories(sql)
    set(SQL_SRC sql/sqlite3.c sql/sql-interface.cpp)
    add_definitions(-DSQL)
endif()

include_directories(./)

foreach(ex "latinsquare" "nonogram")
    add_executable(${ex} models/${ex}.cpp ${SQL_SRC})
    target_link_libraries(${ex} ${Gecode_LIBRARIES})
    if(${BUILD_SQL})
        target_link_libraries(${ex} Threads::Threads ${CMAKE_DL_LIBS})
    endif()
endforeach(ex)
