cmake_minimum_required (VERSION 3.1)

set (CMAKE_CXX_STANDARD 11)

set ( SOURCES
	main.c
	)
	
set ( MOC_HEADERS
	TLCCS.h
	visa.h
	TLCCS.h
	visatype.h
	vpptype.h
	)

add_executable (raman ${SOURCES} ${MOC_HEADERS})
target_link_libraries(raman) 


