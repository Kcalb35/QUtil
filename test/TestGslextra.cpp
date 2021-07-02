#include "gslextra.hpp"
#include "gtest/gtest.h"

TEST(shared_ptr, matrix) {
    auto p =QUtil::gslextra::make_shared_matrix_ptr(2, 2);
    gsl_matrix_set(p.get(),1,1,1.0);
}

TEST(shared_ptr, vector) {
    auto p = QUtil::gslextra::make_shared_vector_ptr(5);
    gsl_vector_set(p.get(),1,1.0);
}

TEST(shared_ptr, matrix_complex) {
    auto p = QUtil::gslextra::make_shared_matrix_complex_ptr(2, 2);
    gsl_matrix_complex_set(p.get(),1,1,gsl_complex{1,1});
}

TEST(shared_ptr,vector_complex){
    auto p = QUtil::gslextra::make_shared_vector_complex_ptr(2);
    gsl_vector_complex_set(p.get(),1,gsl_complex{1,1});
}