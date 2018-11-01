### Overview

Optimizer provides C++ solutions to optimization problems.  The goal is to build an Optimization library for C++ which can provide solutions for single variable and multi-variable optimization, constrained and unconstrained problems. The users should be able to select the algorithm to be used and also access any intermediate data which is obtained whilst running the algorithm.

### A simple first program

```
#include <iostream>
#include <Eigen/Dense>
#include "../Optimizer/optimizer"
using namespace std;
using namespace Eigen;
using namespace Optimizer;

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
}
```
