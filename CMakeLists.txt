cmake_minimum_required(VERSION 3.6)
project(matrix)

include_directories(include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -Wall -Wextra -Weffc++")

set(SOURCE_FILES src/main.cpp src/rational.cpp src/parser.cpp)
add_executable(matrix ${SOURCE_FILES})
