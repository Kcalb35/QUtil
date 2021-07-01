#include "gslextra.h"
#include "gtest/gtest.h"

TEST(shared_ptr, matrix) {
    auto p =QUtil::gslextra::make_shared_matrix_ptr(2, 2);
}

TEST(shared_ptr, vector) {
    auto p = QUtil::gslextra::make_shared_vector_ptr(5);
}

TEST(shared_ptr, matrix_complex) {
    auto p = QUtil::gslextra::make_shared_matrix_complex_ptr(2, 2);
}

TEST(shared_ptr,vector_complex){
    auto p = QUtil::gslextra::make_shared_vector_complex_ptr(1);
}