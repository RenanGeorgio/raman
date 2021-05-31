cmake_minimum_required (VERSION 3.1)

set (CMAKE_CXX_STANDARD 11)

set ( SOURCES
	visa.cpp
	)
	
set ( MOC_HEADERS
	TLCCS.h
	visa.h
	visatype.h
	vpptype.h
	)

add_library (dependencies visa.cpp TLCCS.h visa.h visatype.h vpptype.h)
target_link_libraries(dependencies) 


