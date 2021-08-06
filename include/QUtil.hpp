#ifndef QUTIL_HPP
#define QUTIL_HPP

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_blas.h"
#include "memory"
#include "vector"
#include "fmt/core.h"
#include "fmt/ostream.h"
#include "random"

namespace QUtil {

    namespace gslextra {
        inline std::shared_ptr<gsl_matrix> make_shared_matrix_ptr(const size_t n1, const size_t n2) {
            return std::shared_ptr<gsl_matrix>(gsl_matrix_alloc(n1, n2), [](auto *p) { gsl_matrix_free(p); });
        }

        inline std::shared_ptr<gsl_vector> make_shared_vector_ptr(const size_t n) {
            return std::shared_ptr<gsl_vector>(gsl_vector_alloc(n), [](auto *p) { gsl_vector_free(p); });
        }

        inline std::shared_ptr<gsl_matrix_complex> make_shared_matrix_complex_ptr(const size_t n1, const size_t n2) {
            return std::shared_ptr<gsl_matrix_complex>(gsl_matrix_complex_alloc(n1, n2),
                                                       [](auto *p) { gsl_matrix_complex_free(p); });
        }

        inline std::shared_ptr<gsl_vector_complex> make_shared_vector_complex_ptr(const size_t n) {
            return std::shared_ptr<gsl_vector_complex>(gsl_vector_complex_alloc(n),
                                                       [](auto p) { gsl_vector_complex_free(p); });
        }

        inline gsl_vector **make_vectors(const size_t n, const size_t l) {
            auto v = new gsl_vector *[l];
            for (int i = 0; i < l; ++i) {
                v[i] = gsl_vector_alloc(n);
            }
            return v;
        }

        inline void delete_vectors(gsl_vector **v, const size_t l) {
            for (int i = 0; i < l; ++i) {
                gsl_vector_free(v[i]);
            }
            delete[] v;
        }

        template<typename T>
        inline int sign(T val) {
            return (0 < val) - (0 > val);
        }

        inline std::string format_complex(const gsl_complex c, const int precision) {
            return fmt::format("{1:.{0}f}{2}{3:.{0}f}i", precision, c.dat[0], c.dat[1] < 0 ? "" : "+", c.dat[1]);
        }

        inline std::string format_matrix(gsl_matrix *m, const int precision = 5) {
            auto s = fmt::memory_buffer();
            for (int i = 0; i < m->size1; ++i) {
                for (int j = 0; j < m->size2; ++j) {
                    fmt::format_to(std::back_inserter(s), "{:<{}.{}f}", gsl_matrix_get(m, i, j), precision + 5,
                                   precision);
                }
                fmt::format_to(std::back_inserter(s), "\n");
            }
            return fmt::to_string(s);
        }

        inline std::string format_matrix(gsl_matrix_complex *m, const int precision = 5) {
            auto s = fmt::memory_buffer();
            for (int i = 0; i < m->size1; ++i) {
                for (int j = 0; j < m->size2; ++j) {
                    fmt::format_to(std::back_inserter(s), "{:<{}}",
                                   format_complex(gsl_matrix_complex_get(m, i, j), precision), 2 * precision + 10);
                }
                fmt::format_to(std::back_inserter(s), "\n");
            }
            return fmt::to_string(s);
        }

        inline std::string format_vector(gsl_vector *v, const int precision = 5) {
            auto s = fmt::memory_buffer();
            for (int i = 0; i < v->size; ++i) {
                fmt::format_to(std::back_inserter(s), "{:<{}.{}f}", gsl_vector_get(v, i), precision + 5, precision);
            }
            fmt::format_to(std::back_inserter(s), "\n");
            return fmt::to_string(s);
        }

        inline std::string format_vector(gsl_vector_complex *v, const int precision = 5) {
            auto s = fmt::memory_buffer();
            for (int i = 0; i < v->size; ++i) {
                fmt::format_to(std::back_inserter(s), "{:<{}}", format_complex(gsl_vector_complex_get(v, i), precision),
                               2 * precision + 10);
            }
            fmt::format_to(std::back_inserter(s), "\n");
            return fmt::to_string(s);
        }

        inline int step_function(double x) {
            return x > 0 ? 1 : 0;
        }

        template<typename matrix, typename f, typename ... Args>
        inline void free_gsl(f func, matrix m, Args...args) {
            func(m);
            if constexpr (sizeof...(args) > 0) {
                free_gsl(func, args...);
            }
        }

    }

    namespace QMath {
        inline double integral(gsl_vector *left, gsl_matrix *op, gsl_vector *right, gsl_vector *wb) {
            gsl_vector_set_zero(wb);
            gsl_blas_dgemv(CblasNoTrans, 1, op, right, 0, wb);
            double result;
            gsl_blas_ddot(left, wb, &result);
            return result;
        }

        inline double cal_NAC(gsl_matrix *dh, gsl_vector *s1, gsl_vector *s2, double e1, double e2, gsl_vector *wb) {
            return QUtil::QMath::integral(s1, dh, s2, wb) / (e2 - e1);
        }

        /// set nac matrix
        /// \param nac
        /// \param dh
        /// \param v
        /// \param e
        /// \param wb
        inline void set_NAC_m(gsl_matrix *nac, gsl_matrix *dh, gsl_vector **v, const double e[], gsl_vector *wb) {
            for (int i = 0; i < nac->size1; ++i) {
                gsl_matrix_set(nac, i, i, 0);
            }
            double nac_x = 0;
            for (int i = 0; i < nac->size1; ++i) {
                for (int j = 0; j < i; ++j) {
                    nac_x = QUtil::QMath::cal_NAC(dh, v[i], v[j], e[i], e[j], wb);
                    gsl_matrix_set(nac, i, j, nac_x);
                    gsl_matrix_set(nac, j, i, -nac_x);
                }
            }
        }

        inline void
        diagonalize(gsl_matrix *hamitonian, gsl_vector **v, double e[], gsl_vector *e_value, gsl_matrix *e_vector,
                    gsl_eigen_symmv_workspace *wb) {
            gsl_eigen_symmv(hamitonian, e_value, e_vector, wb);
            gsl_eigen_symmv_sort(e_value, e_vector, GSL_EIGEN_SORT_VAL_ASC);
            for (int i = 0; i < hamitonian->size1; ++i) {
                gsl_matrix_get_col(v[i], e_vector, i);
                e[i] = gsl_vector_get(e_value, i);
            }
        }

        inline void correct_wave_function(gsl_vector *ref, gsl_vector *now) {
            bool flag = false;
            for (int i = 0; i < ref->size; ++i) {
                if (gsl_vector_get(ref, i) * gsl_vector_get(now, i) < 0)
                    flag = true;
            }
            if (flag) gsl_vector_scale(now, -1);
        }
    }

    namespace rng {

        inline double norm_dist(const double avg, const double sigma) {
            static thread_local std::mt19937 generator(std::random_device{}());
            std::normal_distribution<double> distribution(avg, sigma);
            return distribution(generator);
        }

        inline double uni_dist(const double min, const double max) {
            static thread_local std::mt19937 generator(std::random_device{}());
            std::uniform_real_distribution<double> distribution(min, max);
            return distribution(generator);
        }
    }
}
#endif