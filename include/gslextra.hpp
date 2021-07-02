#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "memory"
#include "vector"
#include "fmt/core.h"
#include "fmt/ostream.h"

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

        template<typename T>
        inline int sign(T val) {
            return (0 < val) - (0 > val);
        }

        inline std::string format_complex(const gsl_complex c, const int precision) {
            return fmt::format("{1:.{0}f}{2}{3:.{0}f}i", precision, c.dat[0], c.dat[1] < 0 ? "" : "+", c.dat[1]);
        }

        std::string format_matrix(gsl_matrix *m, const int precision = 5) {
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

        std::string format_matrix(gsl_matrix_complex *m, const int precision = 5) {
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

        std::string format_vector(gsl_vector *v, const int precision = 5) {
            auto s = fmt::memory_buffer();
            for (int i = 0; i < v->size; ++i) {
                fmt::format_to(std::back_inserter(s), "{:<{}.{}f}", gsl_vector_get(v, i), precision + 5, precision);
            }
            fmt::format_to(std::back_inserter(s), "\n");
            return fmt::to_string(s);
        }

        std::string format_vector(gsl_vector_complex *v, const int precision = 5) {
            auto s = fmt::memory_buffer();
            for (int i = 0; i < v->size; ++i) {
                fmt::format_to(std::back_inserter(s), "{:<{}}", format_complex(gsl_vector_complex_get(v, i), precision),
                               2 * precision + 10);
            }
            fmt::format_to(std::back_inserter(s), "\n");
            return fmt::to_string(s);
        }

    }

}