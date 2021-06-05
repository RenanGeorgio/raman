set (CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
set (LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib")
set (EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")

#find_package (opencv 4 REQUIRED)

set ( SOURCES
	main.cpp
	comunication_controller.cpp
	)
	
set ( MOC_HEADERS
	comunication_controller.hpp
	TLCCS.h
	)

add_executable (comunication main.cpp comunication_controller.cpp comunication_controller.hpp)
target_link_libraries(comunication) 

