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
    cout << "The Optimal Point obtained is: ";
    cout << NewtonRapshon (func, range) << endl;;

    cout << "Finding optimal point using above range for Golden Section Search Method." << endl;
    cout << "The Optimal Point obtained is: ";
    cout << GoldenSection (func, range) << endl;;

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

    cout << "Testing Data Save Methods" << endl;
    Optimum point1;
    point1.path = Hessian(SumSquares, Vector3d(3, 3, 4));
    point1.SaveAsTxt();
    point1.SaveAsCsv();
}
