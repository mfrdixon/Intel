cmake_minimum_required(VERSION 3.1)
project(pnf)

set(CMAKE_CXX_STANDARD 17)

include_directories(${PROJECT_SOURCE_DIR})

set(SOURCE_FILES main.cpp)
add_executable(pnf ${SOURCE_FILES})
