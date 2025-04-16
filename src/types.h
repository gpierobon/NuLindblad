#ifndef TYPES_H
#define TYPES_H

#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std::complex_literals;

using Clock   = std::chrono::high_resolution_clock;

using cdouble = std::complex<double>;
using Vec     = Eigen::VectorXd;
using Mat     = Eigen::MatrixXd;

using cVec    = Eigen::VectorXcd;
using cMat    = Eigen::MatrixXcd;
using cDMat   = Eigen::Matrix<cdouble, Eigen::Dynamic, Eigen::Dynamic>;

using SpMat   = Eigen::SparseMatrix<double>;
using cSpMat  = Eigen::SparseMatrix<std::complex<double>>;

using Triplet = Eigen::Triplet<std::complex<double>>;

enum class IntegratorType
{
    Euler,
    RK4, 
    BackwardEuler,
    CrankNicolson 
};

typedef struct 
{
    int N;
    int kry;
    int outs;

    double hthr;
    double omega;
    double gamma_p;
    double gamma_m;
    double h;
    double t_i;
    double t_f;

    IntegratorType integrator;

} Params;
 

Params defaults()
{
    Params pars;
    pars.N     = 10;
    pars.kry   = 15;
    pars.outs  = 20;

    pars.hthr  = 0.5f;
    pars.omega = 1.0f;
    pars.gamma_p = 0.1;
    pars.gamma_m = 0.09;

    pars.h   = 0.01;
    pars.t_i = 1e-10;
    pars.t_f = 1.0f;

    pars.integrator = IntegratorType::RK4;

    return pars;
}

IntegratorType parseIntegrator(const std::string& name)
{
    if (name == "Euler") return IntegratorType::Euler;
    if (name == "RK4")   return IntegratorType::RK4;
    if (name == "BE")    return IntegratorType::BackwardEuler;
    if (name == "CN")    return IntegratorType::CrankNicolson;

    throw std::invalid_argument("Unknown integrator: " + name);
}


int parseArgs(int argc, char* argv[], Params* pars)
{
    for (int i = 1; i < argc; ++i)
    {
        if      (!strcmp(argv[i], "--N") && i+1 < argc)       { pars->N       = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--kr") && i+1 < argc)      { pars->kry     = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--out") && i+1 < argc)     { pars->outs    = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--thr") && i+1 < argc)     { pars->hthr    = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--h") && i+1 < argc)       { pars->h       = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--ti") && i+1 < argc)      { pars->t_i     = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--tf") && i+1 < argc)      { pars->t_f     = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--omega") && i+1 < argc)   { pars->omega   = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--gamma_p") && i+1 < argc) { pars->gamma_p = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--gamma_m") && i+1 < argc) { pars->gamma_m = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--integrator") && i + 1 < argc)
        {
            std::string integrator_str = argv[++i];
            try{
                pars->integrator = parseIntegrator(integrator_str);
            } catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << "\n";
                return 1;
            }
        }
    }

    return 0;
}

void printParams(Params* pars)
{
    std::cout << "Eigen is using " << Eigen::nbThreads() << " threads." << std::endl;
    std::cout << "                                "  << std::endl;
    std::cout << "--------------------------------"  << std::endl;
    std::cout << "N           = " << pars->N         << std::endl; 
    std::cout << "state       = Dicke"               << std::endl; 
    std::cout << "omega       = " << pars->omega     << std::endl; 
    std::cout << "dt          = " << pars->h         << std::endl; 
    std::cout << "ti          = " << pars->t_i       << std::endl; 
    std::cout << "tf          = " << pars->t_f       << std::endl; 
    std::cout << "# meas      = " << pars->outs      << std::endl; 
    std::cout << "--------------------------------"  << std::endl;
}


Params init(int argc, char* argv[])
{
    Params params = defaults();
    parseArgs(argc, argv, &params); 
    printParams(&params);
    return params;
}


#endif 
