#ifndef _STATE_H__
#define _STATE_H__
#include <Eigen/Dense>
#include <types.h>


double logFactorial(int n)
{
    double log_fact = 0.0;
    for (int i = 1; i <= n; ++i)
        log_fact += std::log(i);
    return log_fact;
}


bool isMatrixFinite(const cMat& mat)
{
    return mat.allFinite(); 
}


cMat ProductState(int N)
{
    int dim = N + 1;
    double S = N / 2.0;

    cMat X0 = cMat::Zero(dim, dim);

    std::vector<double> log_fact_table(N + 1);
    log_fact_table[0] = 0.0;
    for (int i = 1; i <= N; ++i)
        log_fact_table[i] = log_fact_table[i - 1] + std::log(i);

    double logNfact = log_fact_table[N];
    double log_pref = -N * std::log(2.0); 

    for (int i = 0; i < dim; ++i)
    {
        double m_i = S - i;
        int idx_i1 = static_cast<int>(S + m_i);
        int idx_i2 = static_cast<int>(S - m_i);
        double log_fac_i = logNfact - log_fact_table[idx_i1] - log_fact_table[idx_i2];

        for (int j = 0; j < dim; ++j)
        {
            double m_j = S - j;
            int idx_j1 = static_cast<int>(S + m_j);
            int idx_j2 = static_cast<int>(S - m_j);
            double log_fac_j = logNfact - log_fact_table[idx_j1] - log_fact_table[idx_j2];

            double log_val = log_pref + 0.5 * (log_fac_i + log_fac_j);
            X0(i, j) = std::exp(log_val);  
        }
    }

    if (!isMatrixFinite(X0))
    {
        std::cerr << "\nProduct state contains non-finite entries, shutting down!\n" << std::endl;
        exit(1);
    }
    return X0;
}


cMat DickeState(int N)
{
    if (N % 2 != 0)
    {
        std::cout << "Can't use odd N!" << std::endl;
        exit(0);
    }

    cMat X0 = cMat::Zero(N + 1, N + 1);
    X0(N / 2, N / 2) = cdouble(1.0, 0.0);
    return X0;
}

void writeState(const cMat& vstate, const std::string& fname)
{
    std::ofstream outFile(fname);
    if (!outFile.is_open())
    {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }

    size_t rows = vstate.rows(), cols = vstate.cols();
    outFile << "# Matrix dimensions: " << rows << " x " << cols << "\n\n";

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            const cdouble& element = vstate(i, j);
            outFile << element.real();
            if (j < cols - 1) outFile << "\t";
        }
        outFile << "\n";
    }

    outFile.close();
    std::cout << "Initial state matrix saved to " << fname << std::endl;
}


cMat initState(Params* pars, bool save = false) 
{
    int       N    = pars->N;
    StateType type = pars->state;
    
    cMat state;

    switch (type)
    {
        case StateType::Product:
            std::cout << "Initialising Product for " << N << " spins ..." << std::endl;
            state = ProductState(N);
            break;
        case StateType::Dicke:
            std::cout << "Initialising Dicke for "   << N << " spins ..." << std::endl;
            state = DickeState(N);
            break;
    }

    if (save) { writeState(state, "InitialState.txt"); }
    return state;
}

#endif
