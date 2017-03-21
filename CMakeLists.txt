cmake_minimum_required(VERSION 3.6)
project(matrix)

include_directories(include)
include_directories(/usr/include/readline)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -Wall -Wextra -Weffc++")

set(SOURCE_FILES src/main.cpp src/rational.cpp src/parser.cpp src/biginteger.cpp src/matrix.cpp src/polynom.cpp)
add_executable(matrix ${SOURCE_FILES})
target_link_libraries(matrix readline)
