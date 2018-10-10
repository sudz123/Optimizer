// Macro for eigen
#ifndef EIGEN_MATRIX_H
#include <Eigen/Dense>
#endif

// Macro for Optimizer
#ifndef OPTIMIZE_GRAD_H
#include "grad.h"
#endif

namespace Optimizer {
    
    Eigen::Vector2d BoundingPhase (std::function<double (double)> obj_func, double x) {
        // The Bounding Phase method
        // Input is a std::function(of the objective function) and the variable of type double
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
        // Input is a std::function(of the objective function) and an Eigen::Vector2d
        // This vector has the range over which the algorithm is evaluated
        // Output is a the Optimum point, type double

        double x = (range(0) + range(1)) / 2;
        double epsilon = 1e-5;
        int it = 0, M = 500;

        Eigen::Vector2d dfx = Derivative(obj_func, x);

        while (std::abs(dfx(0)) > epsilon && dfx(1) != 0 && it < M) {
            x -= dfx(0) / dfx(1);
            dfx = Derivative(obj_func, x);
            ++it;
        }
        return x;
    }

    double GoldenSection (std::function<double (double)> obj_func, Eigen::Vector2d range) {
        // The Golden Section search method
        // Input is a std::function(of the objective function) and an Eigen::Vector2d
        // This vector has the range over which the algorithm is evaluated
        // Output is a the Optimum point, type double

        double epsilon = 1e-5;
        int it = 0, M = 500;
        double l = 1;
        double low = 0, high = 1;
        double w1 = low + 0.618 * l;
        double w2 = high - 0.618 * l;
        double x = (range(0) + range(1)) / 2;
        double f1, f2, x1, x2;
        x1 = w1 * (range(1)-range(0)) + range(0);
        x2 = w2 * (range(1)-range(0)) + range(0);
        f1 = obj_func(x1);
        f2 = obj_func(x2);

        while (it < M && l > epsilon) {
            if(f1 > f2){

                x = x2;
                high = w1;
                l = high - low;
                w1 = w2;
                w2 = high - 0.618 * l;
                x1 = x2;
                x2 = w2 * (range(1)-range(0)) + range(0);
                f1 = f2;
                f2 = obj_func(x2);
                ++it;
                
            }
            else{

                x = x1;
                low = w2;
                l = high - low;
                w2 = w1;
                w1 = low + 0.618 * l;
                x2 = x1;
                x1 = w1 * (range(1)-range(0)) + range(0);
                f2 = f1;
                f1 = obj_func(x1);
                ++it;

            }

        }
        return x;
    }

    enum class BM {
        // Enum class which users can use to select the bracketing method
        // B_PHASE is for Bounding Phase
        // E_SEARCH is for Exhaustive Search
        B_PHASE,
        E_SEARCH
    };

    enum class UDM {
        // Enum class which users can use to select the unidirectional search method
        // N_RAP is for Newton Raphson
        // G_search is for Golden Section Search
        N_RAP,
        G_SEARCH
    };

    double SVOptimize (std::function<double (double)> obj_func,
            double x, BM b_meth = BM::B_PHASE,
            UDM u_search = UDM::N_RAP) {
        // This function does the single variable optimzation
        // Input parameters are :
        // obj_func - The std::function variable containing our objective function
        // x - Our initial point of type double
        // bracketing - Tells us which bracketing method to use. Is of type std::string
        // unidir - Tells us which uni-directional search to use. Is of type std::string
        // This function returns the optimum point(ans) of type double

        Eigen::Vector2d range;

        switch (b_meth) {

            case BM::B_PHASE : range = BoundingPhase(obj_func, x);
                break;
            default : range = BoundingPhase(obj_func, x);

        }

        double ans;

        switch (u_search) {

            case UDM::N_RAP : ans = NewtonRapshon(obj_func, range);
                break;
            case UDM::G_SEARCH : ans = GoldenSection(obj_func, range);
                break;
            default : ans = GoldenSection(obj_func, range);

        }

        return ans;
    }

    Eigen::VectorXd DFP (std::function<double (Eigen::VectorXd)> obj_func,
            Eigen::VectorXd x, int M = 1000, double epsilon1 = 1e-5, double epsilon2 = 1e-5,
            BM b_meth = BM::B_PHASE,
            UDM u_search = UDM::N_RAP) {
        // This function does the multi-variable optimization using the variable metric algorithm
        // Input parameters are :
        // obj_func - The std::function variable containing our objective function
        // x - The initial point from where we begin, of type Eigen::VectorXd
        // M - Number of interations, type int
        // epsilon - termination parameter, type double
        // b_meth - Selects the Bracketing method to be used, type enum class BM
        // u_search - Selects the Unidirectional search method, type enum class UDM
        // Output is of type Eigen::VectorXd, it is the optimised point

        int n = x.size();
        int it = 0;

        Eigen::VectorXd x_prev = x;
        Eigen::VectorXd G = Gradient(obj_func, x), G_prev = G;
        if (G.norm() < epsilon1)
            return x;

        Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n, n);
        Eigen::VectorXd S = -G;

        while (it < M) {

            std::function<double (double)> func = [obj_func, x, S] (double a)->double { return obj_func(x + a * S); };
            double alpha = SVOptimize(func, 0.0, b_meth, u_search);
            x_prev = x;
            x += alpha * S;
            ++it;

            if ((x - x_prev).norm() / x_prev.norm() < epsilon2)
                break;

            G_prev = G;
            G = Gradient(obj_func, x);
            if (G.norm() < epsilon2)
                break;

            Eigen::VectorXd delta_x = x - x_prev;
            Eigen::VectorXd delta_G = G - G_prev;
            A = A + (delta_x * delta_x.transpose())/(delta_x.transpose() * delta_G) -
                (A * delta_G * (A * delta_G).transpose())/(delta_G.transpose() * A * delta_G);

            S = -A * G;

        }

        return x;
    }

    Eigen::VectorXd Cauchy (std::function<double (Eigen::VectorXd)> obj_func,
            Eigen::VectorXd x, int M = 1000, double epsilon1 = 1e-5, double epsilon2 = 1e-5,
            BM b_meth = BM::B_PHASE,
            UDM u_search = UDM::N_RAP) {
        // This function does the multi-variable optimization using the Cauchy's algorithm
        // Input parameters are :
        // obj_func - The std::function variable containing our objective function
        // x - The initial point from where we begin, of type Eigen::VectorXd
        // M - Number of interations, type int
        // epsilon1 - termination parameter, type double
        // epsilon 2 - termination parameter, type double
        // b_meth - Selects the Bracketing method to be used, type enum class BM
        // u_search - Selects the Unidirectional search method, type enum class UDM
        // Output is of type Eigen::VectorXd, it is the optimised point

        int n = x.size();
        int it = 0;

        Eigen::VectorXd x_prev = x;
        Eigen::VectorXd G = Gradient(obj_func, x), G_prev = G;
        Eigen::VectorXd S = -G;

        while (it < M){

            if (G.norm() < epsilon1)
                break;

            std::function<double (double)> func = [obj_func, x, S] (double a)->double { return obj_func(x + a * S); };
            double alpha = SVOptimize(func, 0.0, b_meth, u_search);
            x_prev = x;
            x += alpha * S;
            ++it;

            if ((x - x_prev).norm() / x_prev.norm() < epsilon1)
                break;

            G_prev = G;
            G = Gradient(obj_func, x);

            if ((G * G_prev.transpose()).norm() < epsilon2)
                break;

        }

        return x;
    }

    Eigen::VectorXd Newton (std::function<double (Eigen::VectorXd)> obj_func,
                            Eigen::VectorXd x, int M = 1000, double epsilon1 = 1e-5, double epsilon2 = 1e-5,
                            BM b_meth = BM::B_PHASE,
                            UDM u_search = UDM::N_RAP) {
        // This function does the multi-variable optimization using the Newton's algorithm
        // Input parameters are :
        // obj_func - The std::function variable containing our objective function
        // x - The initial point from where we begin, of type Eigen::VectorXd
        // M - Number of interations, type int
        // epsilon1 - termination parameter, type double
        // epsilon 2 - termination parameter, type double
        // b_meth - Selects the Bracketing method to be used, type enum class BM
        // u_search - Selects the Unidirectional search method, type enum class UDM
        // Output is of type Eigen::VectorXd, it is the optimised point

        int n = x.size();
        int it = 0;

        Eigen::VectorXd x_prev = x;
        Eigen::VectorXd G = Gradient(obj_func, x), G_prev = G;
        Eigen::MatrixXd H = Hessian(obj_func, x);
        Eigen::VectorXd S = - H.inverse() * G;
        
        while (it < M){

            if (G.norm() < epsilon1)
                break;

            std::function<double (double)> func = [obj_func, x, S] (double a)->double { return obj_func(x + a * S); };
            double alpha = SVOptimize(func, 0.0, b_meth, u_search);
            x_prev = x;
            x += alpha * S;
            ++it;

            if ((x - x_prev).norm() / x_prev.norm() < epsilon1)
                break;

            G_prev = G;
            G = Gradient(obj_func, x);
            H = Hessian(obj_func, x);

        }

        return x;
    }

    Eigen::VectorXd Marquardt (std::function<double (Eigen::VectorXd)> obj_func,
                            Eigen::VectorXd x, int M = 1000, double epsilon = 1e-5, double lambda = 1e4) {
        // This function does the multi-variable optimization using the Marquardt's algorithm
        // Input parameters are :
        // obj_func - The std::function variable containing our objective function
        // x - The initial point from where we begin, of type Eigen::VectorXd
        // M - Number of interations, type int
        // epsilon - termination parameter, type double
        // Output is of type Eigen::VectorXd, it is the optimised point

        int n = x.size();
        int it = 0;

        Eigen::VectorXd x_prev = x;
        Eigen::VectorXd G = Gradient(obj_func, x);
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
        Eigen::MatrixXd H = Hessian(obj_func, x);
        Eigen::VectorXd S;
        
        while (it < M){

            if (G.norm() < epsilon)
                break;

            S = - (H + lambda * I).inverse() * G;
            x_prev = x;
            x += S;
            ++it;

            if (obj_func(x_prev) > obj_func(x))
                lambda /= 2;
            else
                lambda *=2;

            G = Gradient(obj_func, x);
            H = Hessian(obj_func, x);

        }

        return x;
    }

    Eigen::VectorXd ConjugateGradient (std::function<double (Eigen::VectorXd)> obj_func,
                                        Eigen::VectorXd x, int M = 1000, double epsilon1 = 1e-5, double epsilon2 = 1e-5,
                                        double epsilon3 = 1e-5, BM b_meth = BM::B_PHASE, UDM u_search = UDM::N_RAP) {
        // This function does the multi-variable optimization using the Conjugate Gradient algorithm
        // Input parameters are :
        // obj_func - The std::function variable containing our objective function
        // x - The initial point from where we begin, of type Eigen::VectorXd
        // M - Number of interations, type int
        // epsilon1 - termination parameter, type double
        // epsilon2 - termination parameter, type double
        // epsilon3 - termination parameter, type double
        // b_meth - Selects the Bracketing method to be used, type enum class BM
        // u_search - Selects the Unidirectional search method, type enum class UDM
        // Output is of type Eigen::VectorXd, it is the optimised point

        int n = x.size();
        int it = 0;

        Eigen::VectorXd x_prev = x;
        Eigen::VectorXd G = Gradient(obj_func, x), G_prev = G;
        Eigen::VectorXd S = - G, S_prev = S;
        std::function<double (double)> func = [obj_func, x, S] (double a)->double { return obj_func(x + a * S); };
        double alpha = SVOptimize(func, 0.0, b_meth, u_search);
        x += alpha * S;
        G = Gradient(obj_func, x);

        while (it < M){

            S = - G + (G.squaredNorm()/G_prev.squaredNorm()) * S_prev;

            std::function<double (double)> func = [obj_func, x, S] (double a)->double { return obj_func(x + a * S); };
            double alpha = SVOptimize(func, 0.0, b_meth, u_search);
            x_prev = x;
            S_prev = S;
            x += alpha * S;
            ++it;

            if ((x - x_prev).norm() / x_prev.norm() < epsilon1)
                break;

            G_prev = G;
            G = Gradient(obj_func, x);

            if (G.norm() < epsilon2)
                break;

        }

        return x;
    }

}
