// Macro for eigen
#ifndef EIGEN_MATRIX_H
#include <Eigen/Dense>
#endif

// Macro for Optimizer
#ifndef OPTIMIZE_GRAD_H
#include "grad.h"
#endif

// The Bounding Phase method
// Input is a function pointer(of the objective function) and the variable of type double
// This variable is the point at which we start the algorithm
// Output is a Eigen::Vector2d
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

// The Newton Raphson method
// Input is a function pointer(of the objective function) and an Eigen::Vector2d
// This vector has the range over which the algorithm is evaluated
// Output is a the Optimum point with +- epsilon accuracy
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
