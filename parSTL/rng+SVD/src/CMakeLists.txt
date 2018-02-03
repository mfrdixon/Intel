project(svd)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${MKL_LIBRARIES} pthread)
