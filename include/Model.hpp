#ifndef MODEL_HPP
#define MODEL_HPP

#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

class NumericalModel {
public:
    virtual void hamitonian_cal(gsl_matrix *m, double x) = 0;

    virtual void d_hamitonian_cal(gsl_matrix *m, double x) = 0;

    virtual double sigma_x(double k) = 0;

    virtual double sigma_p(double k) = 0;

    double x0{}, left{}, right{};
    int DoF{};
};

class ECR : public NumericalModel {
public:
    void hamitonian_cal(gsl_matrix *m, double x) override {
        double h12 = x < 0 ? 0.1 * exp(0.9 * x) : 0.1 * (2 - exp(-0.9 * x));
        gsl_matrix_set(m, 0, 0, 6e-4);
        gsl_matrix_set(m, 1, 1, -6e-4);
        gsl_matrix_set(m, 0, 1, h12);
        gsl_matrix_set(m, 1, 0, h12);
    };

    void d_hamitonian_cal(gsl_matrix *m, double x) override {
        double d12 = 0.1 * 0.9 * exp((x > 0 ? -1 : 1) * 0.9 * x);
        gsl_matrix_set(m, 0, 0, 0);
        gsl_matrix_set(m, 1, 1, 0);
        gsl_matrix_set(m, 1, 0, d12);
        gsl_matrix_set(m, 0, 1, d12);
    };

    double sigma_x(const double k) override {
        return 10 / k;
    }

    double sigma_p(const double k) override {
        return k / 20;
    }

    ECR() {
        x0 = -17.5;
        left = -15;
        right = 15;
        DoF = 2;
    }

};

class SAC : public NumericalModel {
public:
    void hamitonian_cal(gsl_matrix *m, double x) override {
        int flag = (x > 0 ? 1 : -1);
        double h11 = flag * 0.01 * (1 - exp(-flag * 1.6 * x));
        double h12 = 0.005 * exp(-x * x);
        gsl_matrix_set(m, 0, 0, h11);
        gsl_matrix_set(m, 1, 1, -h11);
        gsl_matrix_set(m, 0, 1, h12);
        gsl_matrix_set(m, 1, 0, h12);
    }

    void d_hamitonian_cal(gsl_matrix *m, double x) override {
        double h11 = 0.01 * 1.6 * exp((x < 0 ? 1 : -1) * 1.6 * x);
        double h12 = -2 * 0.005 * x * exp(-x * x);
        gsl_matrix_set(m, 0, 0, h11);
        gsl_matrix_set(m, 1, 1, -h11);
        gsl_matrix_set(m, 1, 0, h12);
        gsl_matrix_set(m, 0, 1, h12);
    };

    double sigma_x(double k) override {
        return 10 / k;
    };

    double sigma_p(double k) override {
        return k / 20;
    }

    SAC() {
        x0 = -17.5;
        left = -10;
        right = 10;
        DoF = 2;
    }
};

class DAC : public NumericalModel {
public:
    void hamitonian_cal(gsl_matrix *m, double x) override {
        double h12 = 0.015 * exp(-0.06 * x * x);
        gsl_matrix_set(m, 0, 0, 0);
        gsl_matrix_set(m, 1, 1, -0.1 * exp(-0.28 * x * x) + 0.05);
        gsl_matrix_set(m, 0, 1, h12);
        gsl_matrix_set(m, 1, 0, h12);
    };

    void d_hamitonian_cal(gsl_matrix *m, double x) override {
        double d12 = -2 * 0.015 * 0.06 * x * exp(-0.06 * x * x);
        gsl_matrix_set(m, 0, 0, 0);
        gsl_matrix_set(m, 1, 1, 2 * 0.1 * 0.28 * x * exp(-0.28 * x * x));
        gsl_matrix_set(m, 1, 0, d12);
        gsl_matrix_set(m, 0, 1, d12);
    };

    double sigma_x(double k) override {
        return 10 / k;
    };

    double sigma_p(double k) override {
        return k / 20;
    };

    DAC() {
        x0 = -17.5;
        left = -15;
        right = 15;
        DoF = 2;
    }
};

class DBG : public NumericalModel {
public:
    void d_hamitonian_cal(gsl_matrix *m, double x) override {
        gsl_matrix_set(m, 0, 0, 0);
        gsl_matrix_set(m, 1, 1, 0);
        double d, z = 10, c = 0.9, b = 0.1;
        if (x < -z) {
            d = b * c * exp(c * (x - z)) - b * c * exp(c * (x + z));
        } else if (x < z) {
            d = b * c * exp(c * (x - z)) - b * c * exp(-c * (x + z));
        } else {
            d = -b * c * exp(-c * (x + z)) + b * c * exp(-c * (x - z));
        }
        gsl_matrix_set(m, 0, 1, d);
        gsl_matrix_set(m, 1, 0, d);
    };

    void hamitonian_cal(gsl_matrix *m, double x) override {
        gsl_matrix_set(m, 0, 0, 6e-4);
        gsl_matrix_set(m, 1, 1, -6e-4);
        double h, z = 10, c = 0.9, b = 0.1;
        if (x < -z) {
            h = b * exp(c * (x - z)) + b * (2 - exp(c * (x + z)));
        } else if (x < z) {
            h = b * exp(c * (x - z)) + b * exp(-c * (x + z));
        } else {
            h = b * exp(-c * (x + z)) + b * (2 - exp(-c * (x - z)));
        }
        gsl_matrix_set(m, 1, 0, h);
        gsl_matrix_set(m, 0, 1, h);
    };

    double sigma_x(double k) override {
        return 3 * sqrt(2) / 2;
    };

    double sigma_p(double k) override {
        return 1 / 3.0 / sqrt(2);
    };

    DBG() {
        x0 = -22.5;
        left = -20;
        right = 20;
        DoF = 2;
    };
};

class DAG : public NumericalModel {
public:

    void hamitonian_cal(gsl_matrix *m, double x) override {
        gsl_matrix_set(m, 0, 0, 6e-4);
        gsl_matrix_set(m, 1, 1, -6e-4);
        double h, b = 0.1, c = 0.9, z = 4;
        if (x < -z) {
            h = -b * exp(c * (x - z)) + b * exp(c * (x + z));
        } else if (x < z) {
            h = -b * exp(c * (x - z)) - b * exp(-c * (x + z)) + 2 * b;
        } else {
            h = b * exp(-c * (x - z)) - b * exp(-c * (x + z));
        }
        gsl_matrix_set(m, 0, 1, h);
        gsl_matrix_set(m, 1, 0, h);
    };

    void d_hamitonian_cal(gsl_matrix *m, double x) override {
        gsl_matrix_set(m, 0, 0, 0);
        gsl_matrix_set(m, 1, 1, 0);
        double d, b = 0.1, c = 0.9, z = 4;
        if (x < -z) {
            d = -b * c * exp(c * (x - z)) + b * c * exp(c * (x + z));
        } else if (x < z) {
            d = -b * c * exp(c * (x - z)) + b * c * exp(-c * (x + z));
        } else {
            d = -b * c * exp(-c * (x - z)) + b * c * exp(-c * (x + z));
        }
        gsl_matrix_set(m, 0, 1, d);
        gsl_matrix_set(m, 1, 0, d);
    };


    double sigma_x(double k) override {
        return 2;
    };

    double sigma_p(double k) override {
        return 0.25;
    };

    DAG() {
        x0 = -27.5;
        left = -20;
        right = 20;
        DoF = 2;
    }
};

class DRN : public NumericalModel {
public:
    void d_hamitonian_cal(gsl_matrix *m, double x) override {
        gsl_matrix_set(m, 0, 0, 0);
        gsl_matrix_set(m, 1, 1, 0);
        double d =
                0.03 * (-2 * 3.2 * ((x - 2) * exp(-3.2 * (x - 2) * (x - 2)) + (x + 2) * exp(-3.2 * (x + 2) * (x + 2))));
        gsl_matrix_set(m, 0, 1, d);
        gsl_matrix_set(m, 1, 0, d);
    };

    void hamitonian_cal(gsl_matrix *m, double x) override {
        gsl_matrix_set(m, 0, 0, 0);
        gsl_matrix_set(m, 1, 1, 0.01);
        double h = 0.03 * (exp(-3.2 * (x - 2) * (x - 2)) + exp(-3.2 * (x + 2) * (x + 2)));
        gsl_matrix_set(m, 0, 1, h);
        gsl_matrix_set(m, 1, 0, h);
    };

    double sigma_x(double k) override {
        return 0.5;
    };

    double sigma_p(double k) override {
        return 1.0;
    };

    DRN() {
        x0 = -12.5;
        left = -10;
        right = 10;
        DoF = 2;
    }
};

#endif
