cmake_minimum_required (VERSION 3.1)

set (CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE "Release")

set ( SOURCES
	main.cpp
	)
	

add_executable (comunication main.cpp)
target_link_libraries(comunication) 

