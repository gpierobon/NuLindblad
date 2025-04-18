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

void printSparseMatrix(const cSpMat& mat) {
    int rows = mat.rows();
    int cols = mat.cols();

    std::cout << "\n\n";
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            if (mat.coeff(i, j) != cdouble(0.0, 0.0))
                std::cout << std::setprecision(3) << mat.coeff(i, j).real() << " ";
            else
                std::cout << " 0.000 ";
        }
        std::cout << std::endl;
    }
}

void getSz(int N, SpMat& Sz)
{
    int   dim = N+1;
    double S  = 0.5*N;
    
    for (int j = 0; j < dim; ++j)
        Sz.insert(j,j) = S-j; 

    Sz.makeCompressed();
}

void getSpm(int N, cSpMat& Sp, cSpMat& Sm)
{
    int    dim = N+1;
    double S   = 0.5*N;

    std::vector<double> m_vals(dim);
    for (int i = 0; i < dim; ++i)
        m_vals[i] = S - i;
    
    for (int j = 0; j < dim - 1; ++j)
    {
        double val = std::sqrt((S - m_vals[j + 1]) * (S + m_vals[j + 1] + 1));

        if (val != 0.0)
        {
            Sp.insert(j, j + 1) = cdouble(val, 0.0);
            Sm.insert(j + 1, j) = cdouble(val, 0.0);
        }
    }

    Sp.makeCompressed();
    Sm.makeCompressed();
    //printSparseMatrix(Sp);
    //printSparseMatrix(Sm);
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

void cache_Jz(Params* pars, double&t, double&Jz)
{
    std::ofstream file(pars->outf, std::ios::app);

    if (file.is_open())
    {
        file << t << " " << Jz << "\n";
        file.close();
    }
    else
        std::cerr << "Error: could not open file for writing.\n";
}

std::vector<size_t> generateLogBins(double t0, double dt, size_t num_steps, size_t n_bins) {
    std::vector<size_t> log_bins;

    double log_t_max = std::log10(t0 + num_steps * dt);
    for (size_t k = 0; k < n_bins; ++k) {
        double log_t = std::log10(t0 + dt);  // Adjusted to the actual start time
        double t_bin = std::pow(10.0, log_t + (log_t_max - log_t) * k / (n_bins - 1));
        size_t i_bin = static_cast<size_t>((t_bin - t0) / dt);
        log_bins.push_back(i_bin);
    }

    return log_bins;
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

    std::cout << std::fixed    <<  std::setprecision(5)
              << std::setw(8)  << "Meas:"     << " " << std::setw(6) << step << "/" << num_steps
              << std::setw(8)  << "t:"        << " " << std::setw(6) << t
              << std::setw(8)  << "<Jz>:"     << " " << std::setw(10) << jz
              << std::setw(6)  << "n:"        << " " << std::setw(8)
              << (std::isnan(slope) ? "  ---" : std::to_string(slope).substr(0, 7))
              << std::setw(4)  << " "         << std::setw(6) << std::setprecision(1) << percent << "%"
              << std::setw(12) << "Elapsed:"  << " " << std::setw(12) << time_str << std::endl;
}

void printStatusx(int step, int num_steps, double t, double jx,  
                  std::chrono::duration<double> dur)
{
    double percent = 100.0 * step / num_steps;
    std::string time_str = formatDuration(dur);

    std::cout << std::fixed    <<  std::setprecision(5)
              << std::setw(8)  << "Meas:"     << " " << std::setw(6) << step << "/" << num_steps
              << std::setw(8)  << "t:"        << " " << std::setw(6) << t
              << std::setw(8)  << "<Jx>:"     << " " << std::setw(10) << jx
              << std::setw(4)  << " "         << std::setw(6) << std::setprecision(1) << percent << "%"
              << std::setw(12) << "Elapsed:"  << " " << std::setw(12) << time_str << std::endl;
}


#endif
