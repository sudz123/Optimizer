// Macro for Eigen
#ifndef EIGEN_MATRIX_H
#include <Eigen/Dense>
#endif

// Macro for Optimizer
#ifndef OPTIMIZE_UTIL_H
#define OPTIMIZE_UTIL_H
#endif

#include <fstream>
#include <vector>

namespace Optimizer {

	class Optimum
    {

    public:
        Eigen::VectorXd optima;
        Eigen::MatrixXd path;
        Eigen::VectorXd grad;
        std::string term_cond;
        int it;
        int func_eval;
     
        void SaveAsTxt()
        {
               //NOTE: REMOVE PATH BEFORE PULL REQUEST
            //create new file
            std::ofstream file_stream("data.txt"); 
            if(file_stream.is_open() == true)
            {
                //go through entire matrix
                for(int i = 0; i < path.cols(); i++)
                {
                    for(int j = 0; j < path.rows(); j++)
                    {
                        //print out whole row
                        file_stream << path(i,j) << " ";
                    }
                    //after every row add new line
                    file_stream << std::endl;
                }
            }
        }

        void SaveAsCsv()
        {
            //create new file
            std::ofstream file_stream("data.csv");
            if(file_stream.is_open() == true)
            {
                //go through entire matrix
               for(int i = 0; i < path.cols(); i++)
                {
                    for(int j = 0; j < path.rows(); j++)
                    {
                        file_stream << path(i,j);
                        //add comma after each coefficient in matrix
                        //all but not last
                        if(j < path.rows()-1)
                        {
                            file_stream << ',';
                        }
                    }
                    //add carrige return after every row
                    file_stream << '\r' << '\n';
                }

            }
        }

    };

    std::vector<int> getFibonacci(int n){

        int f0 = 1, f1 = 1;
        std::vector<int> fib = {f0, f1};

        for(int i=2; i <= n+1; ++i){
            f1 = f1 + f0;
            f0 = f1 - f0;
            fib.push_back(f1);
        }
 
        return fib;
    }
}