cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
project(QUtil)

find_package(GSL REQUIRED)

add_subdirectory(lib/googletest)
add_subdirectory(${PROJECT_SOURCE_DIR}/test)