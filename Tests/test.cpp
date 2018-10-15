#include <iostream>
#include <Eigen/Dense>
#include "../Optimizer/optimizer"
using namespace std;
using namespace Eigen;
using namespace Optimizer;

double func (double x) {
    // This function will return the square after adding 10 to the given number.
    return pow(x + 10, 2);
}

double SumSquares (VectorXd x) {
    return x.squaredNorm();
}

double Matyas (Vector2d x) {
    return 0.26 * (pow(x(0),2) + pow(x(1),2)) - 0.48 * x(0) * x(1);
}

double Dixon (VectorXd x) {
	double val = 0;
	for (int i = 0; i < x.size() - 1; ++i)
		val += 100 * pow(x(i + 1) - pow(x(i), 2), 2) + pow(x(i) - 1, 2);
	return val;
}

double Himmelblau(Vector2d x)
{
    return pow((pow(x(0),2) + x(1) -11),2) + pow((x(0) + pow(x(1),2) - 7),2);
}

int main () {

    cout << "Using Function: (x + 10)^2 for single variable algorithms testing." << endl;
    double ipt = 5.4;

    cout << "Test Bounding Phase:" << endl;
    Vector2d range = BoundingPhase(func, ipt);
    cout << "Range from bounding Phase for initial point :" << ipt << endl;
    cout << range << endl;

    cout << "Test Exhaustive Search:" << endl;
    range = Exhaustive(func, ipt);
    cout << "Range from Exhaustive Search for initial point :" << ipt << endl;
    cout << range << endl;

    cout << "Derivatives at " << ipt << endl;
    cout << Derivative(func, ipt) << endl;

    cout << "Finding optimal point using above range for Newton Rapshon Method." << endl;
    cout << "Optimal Point is: ";
    cout << NewtonRapshon (func, range) << endl;

    cout << "Finding optimal point using above range for Secant Method." << endl;
    cout << "Optimal Point is: ";
    cout << Secant (func, range) << endl;

    cout << "Finding optimal point using above range for Golden Section Search Method." << endl;
    cout << "Optimal Point is: ";
    cout << GoldenSection (func, range) << endl;

    cout << "Finding optimal point using above range for Interval Halving Method." << endl;
    cout << "Optimal Point is: ";
    cout << IntervalHalving (func, range) << endl;

    cout << "Finding optimal point using above range for Fibonnaci Search Method." << endl;
    cout << "Optimal Point is: ";
    cout << Fibonacci (func, range) << endl;

    cout << "Testing SVOptimize on (x + 10)^2 with initial point 5.4." << endl;
    cout << "The Optimal Point obtained is: ";
    cout << SVOptimize(func, 5.4) << endl;

    cout << "Testing Bisection method on (x + 10)^2" << endl;
    cout << "The Optimal Point obtained is: ";
    cout << Bisection(func, Eigen::Vector2d(-20,1) ,0.001, 10000) << endl;

    cout << "Testing Gradient function using SumSquares with 3 variables" << endl;
    cout << Gradient(SumSquares, Vector3d(3, 3, 4)) << endl;

    cout << "Testing Hessian function using SumSquares with 3 variables" << endl;
    cout << Hessian(SumSquares, Vector3d(3, 3, 4)) << endl;

    Vector3d x(4, -3, 7);

    cout << "Testing DFP on sum squared function with 3 variables and initial point (4, -3, 7)." << endl;
    cout << "The Optimal Point obtained is: " << endl;
    cout << DFP(SumSquares, x) << endl;

    cout << "Testing Newton's method on sum squared function with 3 variables and initial point (4, -3, 7)." << endl;
    cout << "The Optimal Point obtained is: " << endl;
    cout << Newton(SumSquares, x) << endl;

    cout << "Testing Cauchy's method on sum squared function with 3 variables and initial point (4, -3, 7)." << endl;
    cout << "The Optimal Point obtained is: " << endl;
    cout << Cauchy(SumSquares, x) << endl;

    cout << "Testing Marquardt's method on sum squared function with 3 variables and initial point (4, -3, 7)." << endl;
    cout << "The Optimal Point obtained is: " << endl;
    cout << Marquardt(SumSquares, x) << endl;

    cout << "Testing Conjugate Gradient method on sum squared function with 3 variables and initial point (4, -3, 7)." << endl;
    cout << "The Optimal Point obtained is: " << endl;
    cout << ConjugateGradient(SumSquares, x) << endl;

    cout << "Testing Conjugate Gradient method on Matyas function with 2 variables and initial point (-3, 7)." << endl;
    cout << "The Optimal Point obtained is: " << endl;
    cout << ConjugateGradient(Matyas, Vector2d(-3, 7)) << endl;

    cout << "Testing Conjugate Gradient method on Dixon function with initial point (-3, 8, 0)." << endl;
    cout << "The Optimal Point obtained is: " << endl;
    cout << ConjugateGradient(Dixon, Vector3d(-3, 8, 0), 30000) << endl;

    cout << "Testing Simplex method Himmelblau function with two variables" << endl;
    cout << "The Optimal Point obtained is: " << endl;
    VectorXd vec(2);
    vec(0) = 1;
    vec(1) = 2;
    cout << SimplexSearch(Himmelblau,vec,10000, 2, 0.5) << endl;

    cout << "Testing Data Save Methods" << endl;
    Optimum point1;
    point1.path = Hessian(SumSquares, Vector3d(3, 3, 4));
    point1.SaveAsTxt();
    point1.SaveAsCsv();
}
