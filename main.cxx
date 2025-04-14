#include <iostream>
#include <complex>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

#include <types.h>
#include <utils.h>
#include <state.h>
#include <lindblad.h>

#define EIGEN_USE_OPENMP

using namespace std::complex_literals;

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

void compute_eigenvalues(const cSpMat& L)
{
    cDMat L_dense = Eigen::Matrix<cdouble, Eigen::Dynamic, Eigen::Dynamic>(L);
    
    Eigen::SelfAdjointEigenSolver<cDMat> solver(L_dense);
    if (solver.info() == Eigen::Success)
    {
        std::cout << "Eigenvalues: " << solver.eigenvalues().transpose() << std::endl;
    }
    else
    {
        std::cerr << "Eigenvalue computation failed!" << std::endl;
    }
}


double getJz_cached(const cMat& mat, const cMat& Sz, const cMat& exp1, const cMat& exp2)
{
    cMat temp = mat * exp1 * Sz * exp2;
    return temp.trace().real();
}

int main(int argc, char* argv[]) {

    //Params pars = init(argc, argv);

    int N = 2;
    for (int i = 1; i < argc; ++i)
    {
        if(!strcmp(argv[i], "--N") && i+1 < argc)
            N = atoi(argv[++i]); 
    }

    double t_i = 0;
    double t_f = 1;
    double h = 0.01;
    double omega = 1.0;

    cMat mstate = DickeState(N);
    cVec vstate = Mat2Vec(mstate);

    int mdim = sqrt(vstate.size());

    Params pars;
    pars.N = N;
    pars.omega = 1.0;
    pars.gamma_p = 0.1;
    pars.gamma_m = 0.09;

    cSpMat L = Lindblad_sparse_pm(&pars);
    L.makeCompressed();

    //compute_eigenvalues(L);

    double rad  = spectral_radius(L);
    double thr = 0.1;
    h = thr/rad;

    size_t new_steps = static_cast<size_t>((t_f-t_i)/h);
    std::cout << "Spectral radius: " << rad << std::endl; 
    std::cout << "Adjusted time step: " << h << std::endl;
    std::cout << new_steps << " steps needed!" << std::endl;

    // Final time vectors
    std::vector<cVec> vsol;
    std::vector<double> jz_vals;

    int dim = N + 1;
    double S = 0.5 * N;

    Vec m_vals(dim);
    for (int i = 0; i < dim; ++i)
        m_vals[i] = S - i;

    Mat Sz(dim, dim);
    Sz.setZero();
    for (int j = 0; j < dim; ++j)
        Sz(j, j) = m_vals[j];

    cMat exp_plus  = expm( 1i * omega * Sz, 100);
    cMat exp_minus = expm(-1i * omega * Sz, 100);

    double t0 = t_i;
    for (size_t i = 0; i < new_steps; ++i)
    {
        double t = t0 + i * h;
        backward_euler_step(vstate, L, h);
        //crank_nicolson_step(vstate, L, h);
        
        if (i % (new_steps / 20) == 0)
        {
            std::cout << "Step " << i << "/" << new_steps << ": ";
            std::cout << "t = " << t << " "; 
            vsol.push_back(vstate);
            cMat mstate = Vec2Mat(vstate,mdim,mdim);
            double jz = getJz_cached(mstate, Sz, exp_plus, exp_minus);
            jz_vals.push_back(jz);
            std::cout << "Jz = " << jz << std::endl;
        }
    }

    return 0;
}








