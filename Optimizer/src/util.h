// Macro for Eigen
#ifndef EIGEN_MATRIX_H
#include <Eigen/Dense>
#endif

// Macro for Optimizer
#ifndef OPTIMIZE_UTIL_H
#define OPTIMIZE_UTIL_H
#endif

namespace Optimizer {

	class Optimum{

    public:
        Eigen::VectorXd optima;
        Eigen::MatrixXd path;
        Eigen::VectorXd grad;
        std::string term_cond;
        int it;
        int func_eval;

        void SaveAsTxt(){

        }

        void SaveAsCsv(){

        }

    };
}