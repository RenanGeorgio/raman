cmake_minimum_required (VERSION 3.1)
project (raman)

set (CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE "Release")

#find_package(opencv 4 REQUIRED)

if(POLICY CMP0020)
	cmake_policy(SET CMP0020 NEW) 
endif()

if(MSVC)
	ADD_DEFINITIONS(/D _CRT_SECURE_NO_WARNINGS)
	ADD_DEFINITIONS(/D _SCL_SECURE_NO_WARNINGS)
	ADD_DEFINITIONS(/D _CRT_NONSTDC_NO_DEPRECATE)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0 ${OpenMP_CXX_FLAGS}")
endif(MSVC)

include_directories (include)
add_subdirectory (include)

include_directories (comunication)
add_subdirectory (comunication)




