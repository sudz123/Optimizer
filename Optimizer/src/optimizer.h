// Macro for eigen
#ifndef EIGEN_MATRIX_H
#include <Eigen/Dense>
#endif

// Macro for Optimizer
#ifndef OPTIMIZE_GRAD_H
#include "grad.h"
#endif

Eigen::Vector2d BoundingPhase (std::function<double (double)> obj_func, double x) {
// The Bounding Phase method
// Input is a function pointer(of the objective function) and the variable of type double
// This variable is the point at which we start the algorithm
// Output is a Eigen::Vector2d
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

double NewtonRapshon (std::function<double (double)> obj_func, Eigen::Vector2d range) {
// The Newton Raphson method
// Input is a function pointer(of the objective function) and an Eigen::Vector2d
// This vector has the range over which the algorithm is evaluated
// Output is a the Optimum point with +- epsilon accuracy
    double x = (range(0) + range(1)) / 2;
    double epsilon = 1e-5;
    int it = 1, M = 500;

    Eigen::Vector2d dfx = Derivative(obj_func, x);

    while (std::abs(dfx(0)) > epsilon && dfx(1) != 0 && it < M) {
        x -= dfx(0) / dfx(1);
        dfx = Derivative(obj_func, x);
        ++it;
    }
    return x;
}

double SVOptimize (std::function<double (double)> obj_func, double x) {
    Eigen::Vector2d range = BoundingPhase(obj_func, x);
    return NewtonRapshon(obj_func, range);
}

Eigen::VectorXd DFP (std::function<double (Eigen::VectorXd)> obj_func, Eigen::VectorXd x) {
    int n = x.size();
    int it = 0, M = 1000;
    double epsilon = 1e-5;

    Eigen::VectorXd x_prev = x;
    Eigen::VectorXd G = Gradient(obj_func, x), G_prev = G;
    if (G.norm() < epsilon)
        return x;

    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n, n);
    Eigen::VectorXd S = -G;

    while (it < M) {
        std::function<double (double)> func = [obj_func, x, S] (double a)->double { return obj_func(x + a * S); };
        double alpha = SVOptimize(func, 0.0);
        x_prev = x;
        x += alpha * S;
        ++it;

        if ((x - x_prev).norm() / x_prev.norm() < epsilon)
            break;

        G_prev = G;
        G = Gradient(obj_func, x);
        if (G.norm() < epsilon)
            break;

        Eigen::VectorXd delta_x = x - x_prev;
        Eigen::VectorXd delta_G = G - G_prev;
        A = A + (delta_x * delta_x.transpose())/(delta_x.transpose() * delta_G) - (A * delta_G * (A * delta_G).transpose())/(delta_G.transpose() * A * delta_G);

        S = -A * G;
    }

    return x;
}
