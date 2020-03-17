cmake_minimum_required(VERSION 2.8)

project(Stellar_Library)
set(stellar_suite_VERSION_MAJOR 1)
set(stellar_suite_VERSION_MINOR 0)

SET(CMAKE_MODULE_PATH ${PROJECT_BINARY_DIR}/cmake/Modules/)
#~ MESSAGE(STATUS "local cmake root: " ${CMAKE_MODULE_PATH})

set(LANGUAGE C++)
#set(CMAKE_VERBOSE_MAKEFILE on)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )

set(CMAKE_BUILD_TYPE Release)

find_package(GMP REQUIRED)
include_directories(${GMP_INCLUDE_DIR})
find_package(BITMAGIC REQUIRED)
include_directories(${BM_INCLUDE_DIR})
#~ MESSAGE(STATUS "GMP LIB: " ${GMP_INCLUDE_DIR})

#if(CMAKE_COMPILER_IS_GNUCXX) #do not work on mac
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "-O3 -std=c++11 -fopenmp") ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "")
#    set(CMAKE_EXE_LINKER_FLAGS "-s -lrt") ## Strip binary 
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
#    set(CMAKE_CXX_FLAGS "-O3 -march=native -std=c++11 -fopenmp -Wall") ## Optimize
    set(CMAKE_CXX_FLAGS "-O3 -march=native -std=c++11 -fopenmp") ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s ")
endif()

#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
#CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
#if(COMPILER_SUPPORTS_CXX11)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#elseif(COMPILER_SUPPORTS_CXX0X)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#else()
#        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
#endif()

################################################################## prova non funzionante
#add_custom_target(clean-cmake-files
#   COMMAND ${CMAKE_COMMAND} -P clean-all.cmake
#)
## clean-all.cmake
#set(cmake_generated ${CMAKE_BINARY_DIR}/CMakeCache.txt
#                    ${CMAKE_BINARY_DIR}/cmake_install.cmake
#                    ${CMAKE_BINARY_DIR}/Makefile
#                    ${CMAKE_BINARY_DIR}/CMakeFiles
#)
#foreach(file ${cmake_generated})
#  if (EXISTS ${file})
#     file(REMOVE_RECURSE ${file})
#  endif()
#endforeach(file)
#add_custom_target(clean-all
#   COMMAND ${CMAKE_BUILD_TOOL} clean
#   COMMAND ${CMAKE_COMMAND} -P clean-all.cmake
#)
############################################################### ===============================

include_directories(sources)
add_subdirectory(sources)

### MAIN OF THE LIBRARY
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/dist )
#add_executable(stellar_suite
#		examples/stellar_suite/main_utility_functions.h
#		examples/stellar_suite/main_utility_functions.cpp
#		examples/stellar_suite/main_simplicial_functions.h
#		examples/stellar_suite/main_simplicial_functions.cpp
#		examples/stellar_suite/main_cp_functions.h
#		examples/stellar_suite/main_cp_functions.cpp
#		examples/stellar_suite/main.cpp
#	)
#target_link_libraries(stellar_suite stellar ${GMP_LIBRARY})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/dist )
add_executable(stellar_lite
                examples/stellar_suite/main_utility_functions.h
                examples/stellar_suite/main_utility_functions.cpp
                examples/stellar_suite/main_simplicial_functions.h
                examples/stellar_suite/main_simplicial_functions.cpp
                examples/stellar_suite/main_cp_functions.h
                examples/stellar_suite/main_cp_functions.cpp
                examples/stellar_suite/main_lite.cpp
        )
target_link_libraries(stellar_lite stellar ${GMP_LIBRARY})

### TEST THE LIBRARY
#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/test )
#~ add_executable(test_index_generation examples/test_spatial_index_generation.cpp)
#~ target_link_libraries(test_index_generation triangle_trees_lib)

#~ add_executable(test_curvature examples/test_curvature.cpp)
#~ target_link_libraries(test_curvature triangle_trees_lib)

#~ add_executable(test_morse examples/test_morse.cpp)
#~ target_link_libraries(test_morse triangle_trees_lib)

#~ add_executable(test_morse_simplification examples/test_morse_simplification.cpp)
#~ target_link_libraries(test_morse_simplification triangle_trees_lib)

#~ add_executable(test_soup examples/test_soup.cpp)
#~ target_link_libraries(test_soup triangle_trees_lib)

#~ add_executable(test_batched_queries examples/test_batched_queries.cpp)
#~ target_link_libraries(test_batched_queries triangle_trees_lib)

#~ add_executable(test_terrain_feature_extraction examples/test_terrain_feature_extraction.cpp)
#~ target_link_libraries(test_terrain_feature_extraction triangle_trees_lib)

#~ add_executable(test_spatial_queries examples/test_spatial_queries.cpp)
#~ target_link_libraries(test_spatial_queries triangle_trees_lib)
