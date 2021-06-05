cmake_minimum_required (VERSION 3.1)

set (CMAKE_CXX_STANDARD 11)

set ( SOURCES
	comunication_controller.cpp
	TLCCS.c
	)
	
set ( MOC_HEADERS
	comunication_controller.hpp
	visa.h
	TLCCS.h
	visatype.h
	vpptype.h
	)

add_library (dependencies ${SOURCES} ${MOC_HEADERS})
target_link_libraries(dependencies) 


