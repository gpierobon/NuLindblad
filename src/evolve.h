#ifndef _KRYLOV_H__
#define _KRYLOV_H__

#include <types.h>

cVec odefun(double t, const cVec& X_vec, const cSpMat& lind)
{
    return lind * X_vec;
}

void Euler_step(cVec& X, const cSpMat& lind, double h)
{
    X += h * (lind * X);
}

void RK4_step(cVec& X, const cSpMat& lind, double t, double h)
{
    cVec k1 = odefun(t, X, lind);
    cVec k2 = odefun(t + 0.5 * h, X + 0.5 * h * k1, lind);
    cVec k3 = odefun(t + 0.5 * h, X + 0.5 * h * k2, lind);
    cVec k4 = odefun(t + h, X + h * k3, lind);

    X += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

void arnoldi(const cSpMat& L, const cVec& v0, cMat& V, cMat& H, int m) {
    int n = v0.size();
    V = cMat::Zero(n, m); // Krylov basis
    H = cMat::Zero(m, m); // Hessenberg matrix

    // Normalize the initial vector
    cVec v = v0 / v0.norm();
    V.col(0) = v;

    for (int j = 0; j < m - 1; ++j) {
        // Apply the operator L to the current vector V.col(j)
        cVec w = L * V.col(j);

        // Orthogonalize the vector against the previous ones
        for (int i = 0; i <= j; ++i) {
            H(i, j) = V.col(i).adjoint() * w;
            w -= H(i, j) * V.col(i);
        }

        // Update the Hessenberg matrix
        H(j + 1, j) = w.norm();
        if (H(j + 1, j).real() < 1e-12) break; // Break if the vector is too small

        // Normalize and store the new Krylov basis vector
        V.col(j + 1) = w / H(j + 1, j);
    }
}

cVec krylov_step(const cSpMat& L, const cVec& rho0, double t, int m) {
    cMat V, H;
    arnoldi(L, rho0, V, H, m);  // Build Krylov subspace

    cMat expH = expm(t * H);

    // Vector of ones, first element corresponds to the norm of the initial state
    cVec e1 = cVec::Zero(m);
    e1(0) = rho0.norm();

    // Compute the result by multiplying the Krylov basis with exp(H * t) * e1
    cVec result = V.leftCols(m) * (expH * e1);
    return result;
}

void backward_euler_step(cVec& X, const cSpMat& L, double h)
{
    const int dim = L.rows();
    cSpMat I(dim, dim);
    I.setIdentity(); 

    cSpMat A = I - h * L;

    Eigen::SparseLU<cSpMat> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Decomposition failed in Backward Euler step");
    }

    X = solver.solve(X);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed in Backward Euler step");
    }
}

void crank_nicolson_step(cVec& X, const cSpMat& L, double h)
{
    const int dim = L.rows();
    cSpMat I(dim, dim);
    I.setIdentity(); 

    cSpMat A = I - (h / 2.0) * L;
    cSpMat B = I + (h / 2.0) * L;
    cVec rhs = B * X;

    Eigen::SparseLU<cSpMat> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Decomposition failed in Crank-Nicolson step");
    }

    X = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed in Crank-Nicolson step");
    }
    
}

void evolve(IntegratorType type, cVec &X, const cSpMat& L, double h, double t)
{
    switch (type)
    {
        case IntegratorType::Euler:
            Euler_step(X, L, h);
            break;

        case IntegratorType::RK4:
            RK4_step(X, L, t, h);
            break;

        case IntegratorType::BackwardEuler:
            backward_euler_step(X, L, h);
            break;

        case IntegratorType::CrankNicolson:
            crank_nicolson_step(X, L, h);
            break;
    }

}

double getJz(const cMat& mat, const cMat& Sz, const cMat& exp1, const cMat& exp2)
{
    cMat temp = mat * exp1 * Sz * exp2;
    return temp.trace().real();
}

double getJzSlope(double t1, double jz1, double t0, double jz0)
{
    if (t1 <= t0 || jz1 == 0.0 || jz0 == 0.0){
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::log(jz1 / jz0) / std::log(t1 / t0);
}

#endif
