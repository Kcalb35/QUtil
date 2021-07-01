#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "memory"

namespace QUtil {

    namespace gslextra {
        std::shared_ptr<gsl_matrix> make_shared_matrix_ptr(const size_t n1, const size_t n2) {
            return std::shared_ptr<gsl_matrix>(gsl_matrix_alloc(n1, n2), [](auto *p) { gsl_matrix_free(p); });
        }

        std::shared_ptr<gsl_vector> make_shared_vector_ptr(const size_t n) {
            return std::shared_ptr<gsl_vector>(gsl_vector_alloc(n), [](auto *p) { gsl_vector_free(p); });
        }

        std::shared_ptr<gsl_matrix_complex> make_shared_matrix_complex_ptr(const size_t n1, const size_t n2) {
            return std::shared_ptr<gsl_matrix_complex>(gsl_matrix_complex_alloc(n1, n2),
                                                       [](auto *p) { gsl_matrix_complex_free(p); });
        }

        std::shared_ptr<gsl_vector_complex> make_shared_vector_complex_ptr(const size_t n) {
            return std::shared_ptr<gsl_vector_complex>(gsl_vector_complex_alloc(n),
                                                       [](auto p) { gsl_vector_complex_free(p); });
        }
    } // namespace gslextra

} // namespace QUtil