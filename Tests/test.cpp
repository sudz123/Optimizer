#include <iostream>
#include <Eigen/Dense>
#include "../Optimizer/optimizer"
using namespace std;
using namespace Eigen;

double func (double x) {
    return pow(x + 10, 2);
}

int main () {
    cout << "Using Function: (x + 10)^2 for single variable algorithms testing." << endl;
    cout << "Test Bounding Phase:" << endl;
    double ipt = 5.4;
    Vector2d range = boundingPhase(func, ipt);
    cout << "Range from bounding Phase for initial point :" << ipt << endl;
    cout << range << endl;

    cout << "Derivatives at " << ipt << endl;
    cout << derivative(func, ipt) << endl;

    cout << "Finding optimal point using above range for Newton Rapshon Method." << endl;
    cout << "Optimal Point is: ";
    cout << newtonRapshon (func, range) << endl;;
}
