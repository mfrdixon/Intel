cmake_minimum_required(VERSION 2.8)
project(thrust_svd)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -mkl=parallel -tbb -g -DPSTL_USE_PARALLEL_POLICIES")
find_package(Boost 1.53 COMPONENTS program_options REQUIRED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/config")
find_package(MKL REQUIRED)

include_directories(${BOOST_INCLUDE_DIRS}
                    ${MKL_INCLUDE_DIRS}
                    "${PROJECT_SOURCE_DIR}/include")

link_directories(${MKL_LIBRARIES} ${Boost_LIBRARIES})

add_subdirectory(src)



