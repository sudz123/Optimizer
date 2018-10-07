#ifndef EIGEN_MATRIX_H
#include <Eigen/Dense>
#endif

#ifndef OPTIMIZE_GRAD_H
#define OPTIMIZE_GRAD_H
#endif

Eigen::Vector2d Derivative (std::function<double (double)> func, double x) {
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

Eigen::VectorXd Gradient (std::function<double (Eigen::VectorXd)> func, Eigen::VectorXd x) {
    int n = x.size();
    Eigen::VectorXd G(n), delta_xi = Eigen::VectorXd::Zero(n);
    double delta = 1e-5;
    for (int i = 0; i < n; ++i) {
        delta_xi(i) = delta;
        G(i) = (func(x + delta_xi) - func(x - delta_xi)) / (2 * delta);
        delta_xi(i) = 0;
    }

    return G;
}
