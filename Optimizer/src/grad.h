#ifndef EIGEN_MATRIX_H
#include <Eigen/Dense>
#endif

#ifndef OPTIMIZE_GRAD_H
#define OPTIMIZE_GRAD_H
#endif

Eigen::Vector2d derivative (double (*func)(double), double x) {
    double delta = 1e-5;
    double fx1 = func(x + delta);
    double fx2 = func(x - delta);
    double fx = func(x);

    Eigen::Vector2d res;
    double a = (fx1 - fx2) / (2 * delta);
    double b = (fx1 + fx2 - 2 * fx) / pow(delta, 2);
    res << a, b;

    return res;
}
