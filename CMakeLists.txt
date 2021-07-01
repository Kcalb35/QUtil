cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
project(QUtil)

set(SOURCES
        src/hello.cpp)

find_package(GSL REQUIRED)

add_subdirectory(lib/googletest)
add_subdirectory(${PROJECT_SOURCE_DIR}/test)

add_library(${PROJECT_NAME} STATIC ${SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)
