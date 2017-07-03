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

#set(OTHER_SRC sql-interface.cpp)

set(EXAMPLES
        "latinsquare"
        "nonogram"
        "langford-number"
        "ortho-latin"
        "pentominoes"
        "qcp"
        )


set(SQL_SRC "")
option(BUILD_SQL "Build SQLite profiling." ON)
if(${BUILD_SQL})
    find_package(Threads REQUIRED)
    set(SQL_SRC sqlite3.c sql-interface.cpp)
    add_definitions(-DSQL)
endif()

foreach(ex "latinsquare" "nonogram" "qcp" "sudoku")
    add_executable(${ex} ${ex}.cpp ${SQL_SRC})
    target_link_libraries(${ex} ${Gecode_LIBRARIES})
    if(${BUILD_SQL})
        target_link_libraries(${ex} Threads::Threads ${CMAKE_DL_LIBS})
    endif()
endforeach(ex)


#add_executable(latinsquare latinsquare.cpp ${OTHER_SRC})
#target_link_libraries(latinsquare ${Gecode_LIBRARIES})
#
#add_executable(nonogram nonogram.cpp)
#target_link_libraries(nonogram ${Gecode_LIBRARIES})
#
#add_executable(langford-number langford-number.cpp ${OTHER_SRC})
#target_link_libraries(langford-number ${Gecode_LIBRARIES})
#
#add_executable(ortho-latin ortho-latin.cpp ${OTHER_SRC})
#target_link_libraries(ortho-latin ${Gecode_LIBRARIES})
#
#add_executable(pentominoes pentominoes.cpp ${OTHER_SRC})
#target_link_libraries(pentominoes ${Gecode_LIBRARIES})
#
#add_executable(qcp qcp.cpp ${OTHER_SRC})
#target_link_libraries(qcp ${Gecode_LIBRARIES})
