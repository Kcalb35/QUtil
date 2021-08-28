#include "QUtil.hpp"
#include "Model.hpp"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_vector.h"
#include "gtest/gtest.h"
#include "fstream"
#include "fmt/ostream.h"

using namespace QUtil::QMath;
using namespace QUtil::gslextra;

TEST(Model, Models) {
    NumericalModel *models[]{new SAC(), new DAC(), new ECR(), new DBG(), new DAG(), new DRN()};
    std::string names[]{"SAC", "DAC", "ECR", "DBG", "DAG", "DRN"};

    auto h = make_shared_matrix_ptr(2, 2);
    auto e_value = make_shared_vector_ptr(2);
    auto e_vector = make_shared_matrix_ptr(2, 2);
    auto nac = make_shared_matrix_ptr(2, 2);
    auto v = make_vectors(2, 2);
    auto tmp_v = make_vectors(2, 2);
    auto wb = gsl_eigen_symmv_alloc(2);
    double e[2];

    for (int i = 0; i < 6; ++i) {
        auto m = models[i];
        std::ofstream file(names[i] + ".txt");

        for (int k = 0; k < 100; ++k) {
            double x = (m->right - m->left) / 100 * k + m->left;
            m->hamitonian_cal(h.get(), x);
            diagonalize(h.get(), tmp_v, e, e_value.get(), e_vector.get(), wb);
            if (k > 0) {
                // correct
                for (int j = 0; j < 2; ++j) {
                    correct_wave_function(v[j], tmp_v[j]);
                }
            }
            for (int j = 0; j < 2; ++j) {
                gsl_vector_memcpy(v[j], tmp_v[j]);
            }
            m->d_hamitonian_cal(h.get(), x);
            set_NAC_m(nac.get(), h.get(), v, e, e_value.get());
            fmt::print(file, "{:.5e} {:.5e} {:.5e} {:.5e}\n", x, e[0], e[1], gsl_matrix_get(nac.get(), 0, 1));
        }
        file.close();
    }

    delete_vectors(v, 2);
    delete_vectors(tmp_v, 2);
    gsl_eigen_symmv_free(wb);
}

TEST(math, integral) {
    auto l = make_shared_vector_ptr(3);
    auto wb = make_shared_vector_ptr(3);
    auto m = make_shared_matrix_ptr(3, 3);
    gsl_matrix_set_all(m.get(), 1);
    gsl_vector_set_all(l.get(), 1);

    EXPECT_EQ(9, integral(l.get(), m.get(), l.get(), wb.get()));
}

TEST(math, integral_z) {
    auto l = make_shared_vector_complex_ptr(2);
    auto m = make_shared_matrix_complex_ptr(2, 2);
    auto wb = make_shared_vector_complex_ptr(2);
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            gsl_matrix_complex_set(m.get(), i, j, gsl_complex{i + 1.0, j + 1.0});
        }
    }
    for (int i = 0; i < 2; ++i) {
        gsl_vector_complex_set(l.get(), i, gsl_complex{1.0 + i, 2.0 + i});
    }
    EXPECT_EQ(56, integral(l.get(), m.get(), l.get(), wb.get()));
}

TEST(math, inner_product) {
    auto l = make_shared_vector_complex_ptr(2);
    for (int i = 0; i < 2; ++i) {
        gsl_vector_complex_set(l.get(), i, gsl_complex{1.0 + i, 2.0 + i});
    }
    EXPECT_EQ(18, inner_product(l.get(), l.get()));
}