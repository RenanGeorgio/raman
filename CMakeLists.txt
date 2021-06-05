cmake_minimum_required (VERSION 3.1)
project (raman)

set (CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE "Release")

#find_package(opencv 4 REQUIRED)

include_directories (dependencies)
add_subdirectory (dependencies)

include_directories (comunication)
add_subdirectory (comunication)




