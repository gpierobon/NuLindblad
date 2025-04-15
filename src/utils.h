#ifndef _UTILS_H__
#define _UTILS_H__

#include <fstream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <types.h>

cVec Mat2Vec(const cMat& mat)
{
    cVec vec = Eigen::Map<const cVec>(mat.data(), mat.size());
    return vec;
}

cMat Vec2Mat(const cVec& vec, int rows, int cols)
{
    assert(vec.size() == rows * cols);
    return Eigen::Map<const cMat>(vec.data(), rows, cols);
}


double spectral_radius(const cSpMat& A)
{
    cVec x = cVec::Random(A.cols()); 
    x.normalize();
    
    for (int i = 0; i < 20; ++i)
    {
        x = A * x; 
        x.normalize();
    }

    cdouble lambda = x.dot(A*x);
    return std::abs(lambda); 
}


cMat expm(const cMat& A, int terms = 20) {
    cMat result = cMat::Identity(A.rows(), A.cols());
    cMat term   = result;

    for (int k = 1; k <= terms; ++k) {
        term = (term * A) / static_cast<double>(k);
        result += term;
    }

    return result;
}

cSpMat expm(const cSpMat& A, int terms = 20)
{
    int rows = A.rows();
    int cols = A.cols();

    cSpMat result(rows, cols);
    result.setIdentity();

    cSpMat term = result;  

    for (int k = 1; k <= terms; ++k) {
        term = term * A / static_cast<double>(k); 
        result += term;  
    }

    return result;
}

void cache_L(const cSpMat& mat, const std::string& fname)
{
    std::ofstream file(fname);
    for (int k = 0; k < mat.outerSize(); ++k)
    {
        for (cSpMat::InnerIterator it(mat, k); it; ++it)
        {
            file << it.row() << "," << it.col() << "," << it.value().real() << "," 
                 << it.value().imag() << "\n";
        }
    }
    file.close();
}

void cache_Jz(int&N, double&t, double&Jz)
{
    std::ostringstream fname;
    fname << "output/Jz_N_" << N << ".txt"; 
    std::ofstream file(fname.str(), std::ios::app);

    if (file.is_open())
    {
        file << t << " " << Jz << "\n";
        file.close();
    }
    else
        std::cerr << "Error: could not open file for writing.\n";
}

#endif
