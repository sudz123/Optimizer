#include <iostream>
#include <Eigen/Dense>
#include "../Optimizer/optimizer"
using namespace std;
using namespace Eigen;

double func (double x) {
    return pow(x + 10, 2);
}

double func2 (VectorXd x) {
    return x.squaredNorm();
}

int main () {
    cout << "Using Function: (x + 10)^2 for single variable algorithms testing." << endl;
    cout << "Test Bounding Phase:" << endl;
    double ipt = 5.4;
    Vector2d range = BoundingPhase(func, ipt);
    cout << "Range from bounding Phase for initial point :" << ipt << endl;
    cout << range << endl;

    cout << "Derivatives at " << ipt << endl;
    cout << Derivative(func, ipt) << endl;

    cout << "Finding optimal point using above range for Newton Rapshon Method." << endl;
    cout << "Optimal Point is: ";
    cout << NewtonRapshon (func, range) << endl;;

    cout << "Testing SVOptimize on (x + 10)^2 with initial point 5.4." << endl;
    cout << "The Optimal Point obtained is: ";
    cout << SVOptimize(func, 5.4) << endl;

    cout << "Testing Gradient function using func2 with 3 variables" << endl;
    cout << Gradient(func2, Vector3d(3, 3, 4)) << endl;

    cout << "Testing DFP on sum squared function with 3 variables and initial point (4, -3, 7)." << endl;
    cout << "The Optimal Point obtained is: " << endl;
    Vector3d x(4, -3, 7);
    cout << DFP(func2, x) << endl;

}
