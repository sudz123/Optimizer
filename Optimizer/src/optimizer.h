#ifndef EIGEN_MATRIX_H
#include <Eigen/Dense>
#endif


#ifndef OPTIMIZE_GRAD_H
#include "grad.h"
#endif


Eigen::Vector2d boundingPhase (double (*obj_func)(double), double x) {
    double delta = 0.5;
    int k = 0, M = 500;
    double fx1 = obj_func(x + delta), fx2 = obj_func(x), fx3 = obj_func(x - delta);
    if (fx3 <= fx2 && fx2 <= fx1)
        delta = -delta;

    double x1 = x + delta, x2 = x, x3 = x - delta;
    fx1 = obj_func(x1);

    while (k < M && fx1 < fx2) {
        k = k + 1;

        x3 = x2; x2 = x1;
        fx2 = fx1;
        x1 = x1 + pow(2, k) * delta;
        fx1 = obj_func(x1);
    }

    Eigen::Vector2d res;
    if (x1 < x3)
        res << x1, x3;
    else
        res << x3, x1;

    return res;
}

double newtonRapshon (double (*obj_func)(double), Eigen::Vector2d range) {
    double x = (range(0) + range(1)) / 2;
    double epsilon = 1e-5;
    int it = 1, M = 500;

    Eigen::Vector2d dfx = derivative(obj_func, x);

    while (std::abs(dfx(0)) > epsilon && dfx(1) != 0 && it < M) {
        x -= dfx(0) / dfx(1);
        dfx = derivative(obj_func, x);
        ++it;
    }
    return x;
}
