cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)
project(QUtil)

set(SOURCES
        src/hello.cpp)

add_library(${PROJECT_NAME} STATIC ${SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)
