#ifndef _UTILS_H__
#define _UTILS_H__

#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
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

void getSz(int N, SpMat& Sz)
{
    int   dim = N+1;
    double S  = 0.5*N;
    
    for (int j = 0; j < dim; ++j)
        Sz.insert(j,j) = S-j; 
    Sz.makeCompressed();
}

void getExp(int N, double omega, cSpMat& expm)
{
    int   dim = N+1;
    double S  = 0.5*N;
    
    std::vector<Triplet> triplets;
    
    for (int j = 0; j < dim; ++j)
    {
        cdouble phase = std::exp(cdouble(0.0, omega * (S-j)));
        triplets.emplace_back(j, j, phase);
    }
    expm.setFromTriplets(triplets.begin(), triplets.end());
    triplets.clear();
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


std::string formatDuration(std::chrono::duration<double> dur) {
    using namespace std::chrono;

    double seconds = dur.count();
    std::ostringstream out;
    out << std::fixed << std::setprecision(1);

    if (seconds < 60.0) {
        out << seconds << "s";
    } else if (seconds < 3600.0) {
        int mins = static_cast<int>(seconds) / 60;
        double rem = seconds - mins * 60;
        out << mins << "m " << rem << "s";
    } else {
        int hrs  = static_cast<int>(seconds) / 3600;
        int mins = (static_cast<int>(seconds) % 3600) / 60;
        double rem = seconds - hrs * 3600 - mins * 60;
        out << hrs << "h " << mins << "m " << rem << "s";
    }

    return out.str();
}


void printStatus(int step, int num_steps, double t, double jz, 
                 std::chrono::duration<double> dur, double slope)
{
    double percent = 100.0 * step / num_steps;
    std::string time_str = formatDuration(dur);

    std::cout << std::fixed    <<  std::setprecision(4)
              << std::setw(8)  << "Meas:"     << " " << std::setw(6) << step << "/" << num_steps
              << std::setw(6)  << "t:"        << " " << std::setw(8) << t
              << std::setw(8)  << "<Jz>:"    << " " << std::setw(8) << jz
              << std::setw(8) << "n:"    << " " << std::setw(8)
              << (std::isnan(slope) ? "  ---" : std::to_string(slope).substr(0, 7))
              << std::setw(4)  << " "         << std::setw(6) << std::setprecision(1) << percent << "%"
              << std::setw(14) << "Elapsed:"  << " " << std::setw(12) << time_str << std::endl;
}


#endif
