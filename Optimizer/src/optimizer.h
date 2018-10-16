#ifndef OPTIMIZER_H
#define OPTIMIZER_H

/* Includes */
#include <Eigen/Dense>
#include "grad.h"
#include "util.h"
#include <iostream>
#include <queue>

#include <vector>

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

    Eigen::Vector2d Exhaustive (std::function<double (double)> obj_func, double x) {
        // The Exhaustive Search method
        // Input is a std::function(of the objective function) and the variable of type double
        // This variable is the point at which we start the algorithm
        // Output is a Eigen::Vector2d

        double delta = 0.5;
        int k = 0, M = 500;
        double x2 = x + delta, x1 = x - delta;
        double fx1 = obj_func(x1), fx = obj_func(x), fx2 = obj_func(x2);
        Eigen::Vector2d res;

        if (fx1 <= fx && fx <= fx2)
            delta = -delta;

        while(k < M){

            if ( fx1 >= fx && fx <= fx2){
                if (x1 < x2)
                    res << x1, x2;
                else
                    res << x2, x1;

                break;
            }

            x1 = x;
            x = x2;
            x2 = x2 + delta;
            fx1 = fx;
            fx = fx2;
            fx2 = obj_func(x2);

        }

        return res;
    }



    double Bisection(std::function<double (double)> obj_func, Eigen::Vector2d range, double tolerance, int max_iter)
    {
    //Bisection method -
    //Input: std::function - function object;
    //range - range over which algorithm will be evaluated
    //tolerance - sets maximum deviation from root f(x) = 0
    //max_iter - maximum iteration before operation cancelation
    //Output: value which differs from a root of f(x)=0 by less than tolerance value

        if(range(0) == range(1))
        {
            std::cout << "Bisection fail: endpoint values cannot be equal to each other." << std::endl;
            return 0.0;//TODO: return error instead of value
        }
        if(Derivative(obj_func, range(0))(0) > 0 ||
        Derivative(obj_func, range(1))(0) < 0)
        {
            //try to flip values
            Eigen::Vector2d buff = range;
            range(0) = buff(1);
            range(1) = buff(0);
            //repeat the test
                if(Derivative(obj_func, range(0))(0) > 0 ||
            Derivative(obj_func, range(1))(0) < 0)
            {
                std::cout << "Bisection fail: f'(range(0)) > 0 and/or f'(range(1)) < 0 does not meet requirments also flipping procedure failed. Aborting..." << std::endl;
                return 0.0;
            }
            //if test has been passed proceed with function
        }

        //Set variables
        int n = 0;
        double c = 0;

        //Check iterator to avoid infinite loop
        while(n <= max_iter)
        {
            //Calcluate mid point of search space
            c = (range(0) + range(1)) / 2;
            //Calculate derivative of a function
            Eigen::Vector2d dfx = Derivative(obj_func,c);

            //Solution statement
            if(std::abs(dfx(0)) <= tolerance)
            {
                return c;
            }
            //step counter up
            n++;
            //Set new, smaller interval
            if(dfx(0) < 0)
            {
                range(0) = c;
            }
            else
            {
                range(1) = c;
            }
        }
        std::cout << "Bisection fail: Max number of iteration exceeded." << std::endl;
        return 0.0;
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

            if (f1 > f2) {

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
            else {

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

    double IntervalHalving (std::function<double (double)> obj_func, Eigen::Vector2d range) {
        // The Interval Halving method
        // Input is a std::function(of the objective function) and an Eigen::Vector2d
        // This vector has the range over which the algorithm is evaluated
        // Output is a the Optimum point, type double

        double epsilon = 1e-5;
        int it = 0, M = 500;
        double x = (range(0) + range(1)) / 2;
        double f1, f2, f, x1, x2, l;
        l = range(1) - range(0);

        while (it < M && l > epsilon) {

            x1 = range(0) + l / 4;
            x2 = range(1) - l/4;
            f1 = obj_func(x1);
            f2 = obj_func(x2);
            f = obj_func(x);

            if(f1 < f){
                range(1) = x;
                x = x1;
            }
            else if(f2 < f){
                range(0) = x;
                x = x2;
            }
            else{
                range(0) = x1;
                range(1) = x2;
            }

            l = range(1) -  range(0);
            ++it;

        }

        return x;
    }

    double Fibonacci (std::function<double (double)> obj_func, Eigen::Vector2d range) {
        // The Fibonacci Search method
        // Input is a std::function(of the objective function) and an Eigen::Vector2d
        // This vector has the range over which the algorithm is evaluated
        // Output is a the Optimum point, type double
        // Be careful not to set higher values of M as double might overflow (IMPORTANT)

        int k = 2, M = 20;
        double x = (range(0) + range(1)) / 2;
        std::vector<int> fib = getFibonacci(M);
        double l = range(1)-range(0), l_star, x1, x2, f1, f2;
        l_star = ((double)fib[M-k] / (double)fib[M]) * l;

        x1 = range(0) + l_star;
        x2 = range(1) - l_star;
        f1 = obj_func(x1);
        f2 = obj_func(x2);

        while(k!=M){

            ++k;
            l_star = ((double)fib[M-k] / (double)fib[M]) * l;

            if (f1 < f2) {

                range(1) = x2;
                x = x1;
                x2 = x1;
                x1 = range(0) + l_star;
                f2 = f1;
                f1 = obj_func(x1);

            }
            else {

                range(0) = x1;
                x = x2;
                x1 = x2;
                x2 = range(1) - l_star;
                f1 = f2;
                f2 = obj_func(x2);

            }

        }

        return x;
    }

    double Secant (std::function<double (double)> obj_func, Eigen::Vector2d range) {
        // The Secant method
        // Input is a std::function(of the objective function) and an Eigen::Vector2d
        // This vector has the range over which the algorithm is evaluated
        // Output is a the Optimum point, type double

        double epsilon = 1e-5;
        int it = 0, M = 500;
        double x;
        Eigen::Vector2d f1, f2, f;
        f2 = Derivative(obj_func, range(1));
        f1 = Derivative(obj_func, range(0));

        while(it < M) {

            x = range(1) - f2(0) * (range(1) - range(0))/(f2(0) - f1(0));
            f = Derivative(obj_func, x);
            if(std::abs(f(0)) < epsilon) {
                break;
            }
            else if(f(0) < 0) {
                range(0) = x;
                f1 = f;
            }
            else {
                range(1) = x;
                f2 = f;
            }
            ++it;
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
        // I_HALVE is for Interval Halving
        N_RAP,
        G_SEARCH,
        I_HALVE
    };

    enum class MVO {
        // Enum class which users can use to select the multivariable search method
        // V_METRIC is for variable metric or DFP
        // CAUCHY is for Cauchy's Algorithm
        // Newton is for Newton's Algorithm
        // Marquardt is for Marquardt' Algorithm
        // Simplex is for Simplex Search
        // C_GRADIENT if for Conjugate Gradient
        V_METRIC,
        CAUCHY,
        NEWTON,
        MARQUARDT,
        SIMPLEX,
        C_GRADIENT
    };

    double SVOptimize (std::function<double (double)> obj_func,
            double x, BM b_meth = BM::B_PHASE,
            UDM u_search = UDM::N_RAP) {
        // This function does the single variable optimzation
        // Input parameters are :
        // obj_func - The std::function variable containing our objective function
        // x - Our initial point of type double
        // b_meth - Tells us which bracketing method to use.
        // u_search - Tells us which uni-directional search to use.
        // This function returns the optimum point(ans) of type double

        Eigen::Vector2d range;

        switch (b_meth) {

            case BM::B_PHASE : range = BoundingPhase(obj_func, x);
                break;
            case BM::E_SEARCH : range = Exhaustive(obj_func, x);
                break;
            default : range = BoundingPhase(obj_func, x);

        }

        double ans;

        switch (u_search) {

            case UDM::N_RAP : ans = NewtonRapshon(obj_func, range);
                break;
            case UDM::G_SEARCH : ans = GoldenSection(obj_func, range);
                break;
            case UDM::I_HALVE : ans = IntervalHalving(obj_func, range);
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
                                        double epsilon3 = 1e-3, BM b_meth = BM::B_PHASE, UDM u_search = UDM::N_RAP) {
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

		//TODO Check best value for epsilon3

        int it = 0;
		std::queue<Eigen::VectorXd> hist;
		hist.push(x);

        Eigen::VectorXd G = Gradient(obj_func, x), G_prev = G;
        Eigen::VectorXd S = - G, S_prev = S;

        while (it < M) {
            std::function<double (double)> func = [obj_func, x, S] (double a)->double { return obj_func(x + a * S); };
            double alpha = SVOptimize(func, 0.0, b_meth, u_search);
			hist.push(x);
			if (hist.size() > 4)
				hist.pop();
            x += alpha * S;
            ++it;

            if ((x - hist.back()).norm() / hist.back().norm() < epsilon1)
                break;

            G_prev = G;
            G = Gradient(obj_func, x);

            if (G.norm() < epsilon2)
                break;

            S_prev = S;
            S = - G + (G.squaredNorm()/G_prev.squaredNorm()) * S_prev;
			if (acos((double)(S.transpose() * S_prev)/(S.norm() * S_prev.norm())) < epsilon3 && hist.size() == 4) {
				x = hist.front();
				while (!hist.empty())
					hist.pop();

				G = Gradient(obj_func, x), G_prev = G;
				S = -G, S_prev = S;
			}
        }

        return x;
    }

    // NOT WORKING
    Eigen::VectorXd Simplex (std::function<double(Eigen::VectorXd)> obj_func, Eigen::VectorXd var, int M = 1000, double gamma = 2, double beta = 0.5, double epsilon = 1e-5) {
        // This function does the multi-variable optimization using the Simplex Search (Nedler-Mead) algorithm
        // Input parameters are :
        // obj_func - multivariable function that this algorithm is searching optimum for
        // var - point of initial search
        // gamma - Expansion coefficient needs be > 1
        // beta - shrink coefficient needs to be in (0,1) range
        // epsilon - termination criteria variable
        // M - variable used to terminate algorithm after number of iterations
        // Returns: Eigen::VectorXd of n-variables as passed in var

        Eigen::VectorXd null_vec(1);
        null_vec(0) = 0;

        //param check
        if(gamma <= 1)
        {
            std::cout << "SimplexSearch: Gamma parameter should be > 1, aborting...";
            return null_vec;
        }
        if(beta > 1 || beta < 0)
        {
            std::cout << "SimplexSearch: Beta parameter should be in (0,1) range. Aborting...";
            return null_vec;
        }
        //Creating initial simplex: see theory pdf p 114 for details
        //Get number of dimensions/variables 
        int n = var.size();
        double delta = 0;
        //calculate delta value see 114 for details appendix 3 
        if(n == 3)
        {
            delta = 0.25;
        }
        else
        {
            delta = (sqrt(n+1)-2.0)/((double)(n-3));
        }
        //Create initial simplexes, consisting of n+1 n-dimensional points 
        Eigen::MatrixXd points(n+1,n);
        std::vector<std::pair<double,int>> func_vals;
        //Add base point to simplex vector

        points.row(0) = var;

        double C = 2.4;
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n-1; j++)
            {
                if(i == j)
                {
                    points(i,j) = points(0,j) + C;
                }
                else
                {
                    points(i,j) = points(0,j) + C*delta;
                }
                //Calculate function value at currently processed points
                func_vals.push_back(std::make_pair(obj_func(points.row(i)), i));
            }
        }

        int iter = 0;
        while(true)
        {
            //STEP 2
            //Find worst point, best point, next to worst point
            //Sort func_val sector and assign special values
            std::sort(func_vals.begin(), func_vals.end()-1);
            //worst point (highest)
            Eigen::VectorXd xh(n);
            xh.col(0) = points.row(func_vals.back().second);
            //best point (lowest)
            Eigen::VectorXd xl(n);
            xl.col(0) = points.row(func_vals[0].second);
            //next to worst point
            Eigen::VectorXd xg(n);
            xg.col(0) = points.row(func_vals[func_vals.size()-2].second);
            //central point 1 row with n variabvles
            Eigen::VectorXd xc(n);
            //i < N because we exclude worst point from this sum(the worst point is also the last point in this vectpr)
         
            for(int i = 0; i < n; i++)
            {
                for(int j = 0; j < n-1; j++)
                {
                    xc(j) += points(func_vals[i].second, n-1);
                }
            }
            //multiply every variable of xc by 1/n
            xc *= 1/n;
            //STEP 3
            //Check for expansion, reflection and contraction in two types
            Eigen::VectorXd x_new(n);
            //Calculate reflected point
            Eigen::VectorXd xr(n);
            xr = 2*xc-xh;
            x_new = xr;
            //define function values to avoid unnecessary function calls
            double fxr = obj_func(xr);
            double fxl = obj_func(xl);
            double fxg = obj_func(xg);
            double fxh = obj_func(xh);
            if(fxr < fxl)
            {
                x_new = (1+gamma)*xc-gamma*xh;
            }
            else if(fxr >= fxh)
            {
                x_new = (1-beta)*xc + beta*xh;
            }
            else if(fxg < fxr && fxr < fxh)
            {
                x_new = (1+beta)*xc - beta*xh;
            }
            xh = x_new;
            //STEP 4
            //Check termination
            double termination_sum = 0;
            for(int i = 0; i < n-1; i++)
            {
                termination_sum += pow(obj_func(points.row(i)) - obj_func(xc), 2) / (n+1);
            }
            if(pow(termination_sum,0.5) <= epsilon)
            {
                //return optimal point
                return x_new;
            }
            iter++;
            if(iter >= M)
            {  
                std::cout << "SimplexSearch: exceeded number of iterations, terminating" << std::endl;
                return x_new;
            }
        }

        return null_vec;
    }

    Eigen::VectorXd MVOptimize (std::function<double (Eigen::VectorXd)> obj_func,
            Eigen::VectorXd x, MVO m_opt = MVO::MARQUARDT, BM b_meth = BM::B_PHASE,
            UDM u_search = UDM::N_RAP) {
        // This function does the Multi-variable optimzation
        // Input parameters are :
        // obj_func - The std::function variable containing our objective function
        // x - Our initial point of type Eigen::VectorXd
        // m_opt - Tells us which multi-variable search algorithm to use
        // b_meth - Tells us which bracketing method to use.
        // u_search - Tells us which uni-directional search to use.
        // This function returns the optimum point(ans) of type Eigen::VectorXd

        Eigen::VectorXd ans;

        switch (m_opt) {

            case MVO::V_METRIC : ans = DFP(obj_func, x);
                break;
            case MVO::CAUCHY : ans = Cauchy(obj_func, x);
                break;
            case MVO::NEWTON : ans = Newton(obj_func, x);
                break;
            case MVO::MARQUARDT : ans = Marquardt(obj_func, x);
                break;
            case MVO::SIMPLEX : ans = Simplex(obj_func, x);
                break;
            case MVO::C_GRADIENT : ans = ConjugateGradient(obj_func, x);
                break;
            default : ans = Marquardt(obj_func, x);

        }

        return ans;
    }

    // TODO
    Eigen::VectorXd PenaltyConstrained(std::function<double(Eigen::VectorXd)> obj_func, Eigen::VectorXd x, std::vector<std::function<double (Eigen::VectorXd)>> ineq_const,
                                       std::vector<std::function<double (Eigen::VectorXd)>> eq_const){

        return x;

    }

    // TODO
    Eigen::VectorXd MultiplierConstrained(std::function<double(Eigen::VectorXd)> obj_func, Eigen::VectorXd x, std::vector<std::function<double (Eigen::VectorXd)>> ineq_const,
                                        std::vector<std::function<double (Eigen::VectorXd)>> eq_const){

        return x;

    }
    

}

#endif   /* OPTIMIZER_H */
