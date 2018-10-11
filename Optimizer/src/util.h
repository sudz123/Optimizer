// Macro for Eigen
#ifndef EIGEN_MATRIX_H
#include <Eigen/Dense>
#endif

// Macro for Optimizer
#ifndef OPTIMIZE_UTIL_H
#define OPTIMIZE_UTIL_H
#endif

//Libs for file operations
#ifndef _IOSTREAM_
#include <iostream>
#endif

#include <fstream>

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
}