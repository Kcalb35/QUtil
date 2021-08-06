#include "QUtil.hpp"
#include "gtest/gtest.h"

using namespace QUtil::gslextra;
using namespace QUtil::QMath;

TEST(shared_ptr, matrix) {
    auto p = make_shared_matrix_ptr(2, 2);
    gsl_matrix_set(p.get(), 1, 1, 1.0);
}

TEST(shared_ptr, vector) {
    auto p = make_shared_vector_ptr(5);
    gsl_vector_set(p.get(), 1, 1.0);
}

TEST(shared_ptr, matrix_complex) {
    auto p = make_shared_matrix_complex_ptr(2, 2);
    gsl_matrix_complex_set(p.get(), 1, 1, gsl_complex{1, 1});
}

TEST(shared_ptr, vector_complex) {
    auto p = make_shared_vector_complex_ptr(2);
    gsl_vector_complex_set(p.get(), 1, gsl_complex{1, 1});
}

TEST(math, sign) {
    EXPECT_EQ(sign(1), 1);
    EXPECT_EQ(sign(1.0), 1);
    EXPECT_EQ(sign(-1), -1);
    EXPECT_EQ(sign(-1.0), -1);
    EXPECT_EQ(sign(0), 0);
    EXPECT_EQ(sign(0.0), 0);
}

TEST(log, matrix) {
    auto p = make_shared_matrix_ptr(2, 2);
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            gsl_matrix_set(p.get(), i, j, i + 0.5 * j);
        }
    }
    fmt::print(format_matrix(p.get(), 3));
}

TEST(log, matrix_complex) {
    auto p = make_shared_matrix_complex_ptr(2, 2);
    gsl_matrix_complex_set(p.get(), 0, 0, gsl_complex{1, 2.5});
    gsl_matrix_complex_set(p.get(), 0, 1, gsl_complex{-1, -2.5});
    gsl_matrix_complex_set(p.get(), 1, 0, gsl_complex{0, 1.5});
    gsl_matrix_complex_set(p.get(), 1, 1, gsl_complex{1, 0});
    fmt::print(format_matrix(p.get(), 3));
}

TEST(log, vecotr) {
    auto p = make_shared_vector_ptr(5);
    for (int i = 0; i < p->size; ++i) {
        gsl_vector_set(p.get(), i, 2.0 * i);
    }
    fmt::print(format_vector(p.get(), 3));
}

TEST(log, vector_complex) {
    auto p = make_shared_vector_complex_ptr(4);
    gsl_vector_complex_set(p.get(), 0, gsl_complex{1, 2.5});
    gsl_vector_complex_set(p.get(), 1, gsl_complex{-1, -2.5});
    gsl_vector_complex_set(p.get(), 2, gsl_complex{0, 1.5});
    gsl_vector_complex_set(p.get(), 3, gsl_complex{1, 0});
    fmt::print(format_vector(p.get(), 3));
}

TEST(math, integral) {
    auto l = make_shared_vector_ptr(3);
    auto wb = make_shared_vector_ptr(3);
    auto m = make_shared_matrix_ptr(3, 3);
    gsl_matrix_set_all(m.get(), 1);
    gsl_vector_set_all(l.get(), 1);

    EXPECT_EQ(9, integral(l.get(), m.get(), l.get(), wb.get()));
}