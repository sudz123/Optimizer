// Macro for Eigen
#ifndef EIGEN_MATRIX_H
#include <Eigen/Dense>
#endif

// Macro for Optimizer
#ifndef OPTIMIZE_GRAD_H
#define OPTIMIZE_GRAD_H
#endif

// Function to calculate derivative using central differencing
// Inputs are a function pointer and a double variable
// This variable holds the value at which we have to apply central differencing
// Returns an Eigen::Vector2d
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
