cmake_minimum_required(VERSION 3.22)
project(MONI_K)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "g++")

include_directories("~/lib")
include_directories("~/include")
include_directories("~/include/sdsl")

link_directories("~/lib/")

link_libraries(sdsl divsufsort divsufsort64)

add_executable(construct_tables construct_tables.cpp construct_tables.h)
add_executable(MONI_K MONI_K.cpp Tables.cpp Tables.h)
add_executable(bin_table_test Tables.cpp Tables.h)
