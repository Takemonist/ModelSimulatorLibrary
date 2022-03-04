cmake_minimum_required(VERSION 3.0.0)
project(MSL VERSION 0.1.0)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

include(CTest)
enable_testing()

add_executable(MSL main.cpp )

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

find_package(Python3 COMPONENTS Development NumPy)
target_include_directories(MSL PRIVATE ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS} ./MSL ./MSL/NMPC ./MSL/NMPC/matrix ./MSL/NMPC/solver )
target_link_libraries(MSL Python3::Python Python3::NumPy)