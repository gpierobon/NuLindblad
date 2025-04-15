#ifndef _KRYLOV_H__
#define _KRYLOV_H__

#include <types.h>

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

cVec krylov_expm(const cSpMat& L, const cVec& rho0, double t, int m) {
    cMat V, H;
    arnoldi(L, rho0, V, H, m);  // Build Krylov subspace

    // Compute the matrix exponential of the small Hessenberg matrix H * t
    cMat expH = expm(t * H);  // Eigen's matrix exponential (for small dense matrices)

    // Vector of ones, first element corresponds to the norm of the initial state
    cVec e1 = cVec::Zero(m);
    e1(0) = rho0.norm();

    // Compute the result by multiplying the Krylov basis with exp(H * t) * e1
    cVec result = V.leftCols(m) * (expH * e1);
    return result;
}

void crank_nicolson_step(cVec& X, const cSpMat& L, double h)
{
    const int dim = L.rows();
    cSpMat I(dim, dim);
    I.setIdentity(); 

    // Construct A = (I - h/2 * L)
    cSpMat A = I - (h / 2.0) * L;

    // Construct B = (I + h/2 * L)
    cSpMat B = I + (h / 2.0) * L;

    // Compute RHS: B * X
    cVec rhs = B * X;

    // Solve A * X_new = rhs
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

void backward_euler_step(cVec& X, const cSpMat& L, double h)
{
    const int dim = L.rows();
    cSpMat I(dim, dim);
    I.setIdentity(); 

    // Construct (I - h * L)
    cSpMat A = I - h * L;

    // Solve (I - h * L) * X_new = X
    Eigen::SparseLU<cSpMat> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Decomposition failed in Backward Euler step");
    }

    // Solve for X_new
    X = solver.solve(X);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed in Backward Euler step");
    }
}

Vec odefun(double t, const Vec& X_vec, const SpMat& lind)
{
    return lind * X_vec;
}

void RK4(Vec& X, const SpMat& lind, double t, double h)
{
    Vec k1 = odefun(t, X, lind);
    Vec k2 = odefun(t + 0.5 * h, X + 0.5 * h * k1, lind);
    Vec k3 = odefun(t + 0.5 * h, X + 0.5 * h * k2, lind);
    Vec k4 = odefun(t + h, X + h * k3, lind);
    X += (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
}

double getJz(const cMat& mat, const cMat& Sz, const cMat& exp1, const cMat& exp2)
{
    cMat temp = mat * exp1 * Sz * exp2;
    return temp.trace().real();
}

#endif
