#ifndef TYPES_H
#define TYPES_H
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <complex>

using cdouble = std::complex<double>;
using Vec     = Eigen::VectorXd;
using Mat     = Eigen::MatrixXd;

using cVec    = Eigen::VectorXcd;
using cMat    = Eigen::MatrixXcd;
using cDMat   = Eigen::Matrix<cdouble, Eigen::Dynamic, Eigen::Dynamic>;

using SpMat   = Eigen::SparseMatrix<double>;
using cSpMat  = Eigen::SparseMatrix<std::complex<double>>;

using Triplet = Eigen::Triplet<std::complex<double>>;


typedef struct 
{
    int N;

    float omega;
    float gamma_p;
    float gamma_m;

} Params;
 

#endif 
