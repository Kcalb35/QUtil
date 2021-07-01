enable_testing()

add_executable(TestGslextra TestGslextra.cpp)
target_include_directories(TestGslextra PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/lib/googletest/googletest/include/gtest)
target_link_libraries(TestGslextra gtest_main GSL::gsl)
add_test(
        NAME TestGslextra
        COMMAND TestGslextra
)
