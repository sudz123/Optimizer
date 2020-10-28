// Macro for Optimizer
#ifndef OPTIMIZE_GRAD_H
#define OPTIMIZE_GRAD_H

// Macro for Eigen
#include <Eigen/Dense>
namespace Optimizer {

    Eigen::Vector2d Derivative (std::function<double (double)> func, double x) {
        // Function to calculate derivative using central differencing
        // Inputs are a function pointer and a double variable
        // This variable holds the value at which we have to apply central differencing
        // Returns an Eigen::Vector2d, which first element is single derivative of a func
        // and the second element is second derivative
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

    Eigen::MatrixXd Hessian (std::function<double (Eigen::VectorXd)> func, Eigen::VectorXd x) {

        int n = x.size();
        Eigen::MatrixXd H(n, n), delta_xi = Eigen::VectorXd::Zero(n), delta_xj = Eigen::VectorXd::Zero(n);
        double delta = 1e-5;
        double fx = func(x);

        for(int i = 0; i < n; ++i){
            for(int j = 0; j < n; ++j){

                delta_xi(i) = delta;
                delta_xj(j) = delta;
                
                if(i == j){
                    H(i, j) = (func(x + delta_xi) + func(x - delta_xi) - 2*fx ) / pow(delta, 2); 
                }
                else{
                    H(i, j) =  (func(x + delta_xi + delta_xj) + func(x - delta_xi - delta_xj) - func(x + delta_xi - delta_xj) - func(x - delta_xi + delta_xj) ) / (4 * pow(delta, 2) ) ; 
                }

                delta_xj(j) = 0;

            }

            delta_xi(i) = 0;

        }

        return H;
    }



}

#endif
