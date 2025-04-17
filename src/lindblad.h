#ifndef _LINDBLAD_H__
#define _LINDBLAD_H__

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <chrono>

#include "types.h"

cSpMat kron_sparse(const cSpMat& A, const cSpMat& B)
{
    std::vector<Triplet> triplets;
    triplets.reserve(A.nonZeros() * B.nonZeros());

    int m = A.rows(), n = A.cols();
    int p = B.rows(), q = B.cols();

    for (int i = 0; i < A.outerSize(); ++i)
    {
        for (cSpMat::InnerIterator itA(A, i); itA; ++itA)
        {
            int ai = itA.row(), aj = itA.col();
            std::complex<double> a_val = itA.value();
            for (int k = 0; k < B.outerSize(); ++k)
            {
                for (cSpMat::InnerIterator itB(B, k); itB; ++itB)
                {
                    int bi = itB.row(), bj = itB.col();
                    std::complex<double> b_val = itB.value();
                    triplets.emplace_back(ai * p + bi, aj * q + bj, a_val * b_val);
                }
            }
        }
    }

    cSpMat result(m * p, n * q);
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

cSpMat diag_sparse(const Eigen::VectorXcd& diag)
{
    std::vector<Triplet> triplets;
    for (int i = 0; i < diag.size(); ++i) {
        if (diag[i] != std::complex<double>(0.0, 0.0)) {
            triplets.emplace_back(i, i, diag[i]);
        }
    }
    cSpMat D(diag.size(), diag.size());
    D.setFromTriplets(triplets.begin(), triplets.end());
    return D;
}


cSpMat Lindblad_sparse_pm(Params* pars) {
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "\n";
    std::cout << "Computing the Lindbladian operator (sparse) ... ";
    fflush(stdout);

    int N = pars->N;
    double gamma_p = pars->gamma_p;
    double gamma_m = pars->gamma_m;

    double S = N / 2.0;
    int dim = N + 1;

    Vec m_vals(dim);
    for (int i = 0; i < dim; ++i) {
        m_vals(i) = S - i;
    }

    cSpMat I(dim, dim);
    I.setIdentity();

    cSpMat Sz = diag_sparse(m_vals.cast<cdouble>());
    cSpMat Sz2 = Sz * Sz;

    std::vector<Triplet> sp_triplets, spm_triplets;
    for (int i = 0; i < dim - 1; ++i) {
        cdouble val = std::sqrt((S - m_vals(i + 1)) * (S + m_vals(i + 1) + 1));
        sp_triplets.emplace_back(i, i + 1, val);

        cdouble val2 = (S - m_vals(i + 1)) * (S + m_vals(i + 1) + 1);
        spm_triplets.emplace_back(i, i, val2);
    }

    cSpMat Sp(dim, dim), Spm(dim, dim);
    Sp.setFromTriplets(sp_triplets.begin(), sp_triplets.end());
    Spm.setFromTriplets(spm_triplets.begin(), spm_triplets.end());

    cSpMat Sm = Sp.transpose();

    cVec diag_spm(dim);
    for (int i = 0; i < dim; ++i)
        diag_spm(i) = Spm.coeff(i, i);
    cVec diag_smp = diag_spm.reverse();
    cSpMat Smp = diag_sparse(diag_smp.cast<cdouble>());

    cSpMat L1 = kron_sparse(Sp.transpose(), Sm)
               - 0.5 * kron_sparse(I, Spm)
               - 0.5 * kron_sparse(Spm.transpose(), I);

    cSpMat L2 = kron_sparse(Sm.transpose(), Sp)
               - 0.5 * kron_sparse(I, Smp)
               - 0.5 * kron_sparse(Smp.transpose(), I);

    cSpMat L = gamma_p * L1 + gamma_m * L2;

    auto t2 = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "took " << 1e-3 * dur.count() << " seconds" << std::endl;

    L.makeCompressed();

    return L;
}

#endif
