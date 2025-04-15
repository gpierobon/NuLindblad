#include <iostream>
#include <complex>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <types.h>
#include <utils.h>
#include <state.h>
#include <lindblad.h>
#include <evolve.h>

#define EIGEN_USE_OPENMP


int main(int argc, char* argv[])
{

    Params pars = init(argc, argv);

    bool cacheL = false;
    int  dim    = pars.N + 1;
    double S    = 0.5 * pars.N;
    double t0   = pars.t_i;

    cMat mstate = DickeState(pars.N);
    cVec vstate = Mat2Vec(mstate);
    int mdim    = sqrt(vstate.size());

    cSpMat L = Lindblad_sparse_pm(&pars);
    if (cacheL) {cache_L(L, "Lindblad.txt"); }
    L.makeCompressed();

    double rad     = spectral_radius(L);
    double thr     = 0.1;
    double h_check = thr/rad;

    std::cout << "Maximum time step: " << h_check << std::endl;
    
    // Precompute Sz, exp_+, exp_- for Jz trace
    Vec m_vals(dim);
    for (int i = 0; i < dim; ++i)
        m_vals[i] = S - i;

    fflush(stdout);
    Mat Sz(dim, dim);
    Sz.setZero();
    for (int j = 0; j < dim; ++j)
        Sz(j, j) = m_vals[j];

    std::cout << "Pre-computing exps for time evolution (sparse) ... ";
    std::vector<Triplet> triplets;
    for (int j = 0; j < dim; ++j) {
        std::complex<double> phase = std::exp(std::complex<double>(0.0, pars.omega * m_vals[j]));
        triplets.emplace_back(j, j, phase);
    }

    cSpMat exp_plus(dim, dim);
    exp_plus.setFromTriplets(triplets.begin(), triplets.end());

    triplets.clear();
    for (int j = 0; j < dim; ++j) {
        std::complex<double> phase = std::exp(std::complex<double>(0.0, -pars.omega * m_vals[j]));
        triplets.emplace_back(j, j, phase);
    }
    cSpMat exp_minus(dim, dim);
    exp_minus.setFromTriplets(triplets.begin(), triplets.end());

    std::cout << " done!" << std::endl;

    std::vector<cVec> vsol;
    std::vector<double> jz_vals;
    
    double dt = std::min(pars.h, h_check);
    size_t num_steps = static_cast<size_t>((pars.t_f-pars.t_i)/dt);
    std::cout << "\n";
    std::cout << "Running a loop with " << num_steps << " steps ... \n" << std::endl;

    for (size_t i = 0; i < num_steps; ++i)
    {
        double t = t0 + i * dt;
        cVec ev_state = krylov_expm(L, vstate, dt, pars.kry);
        //backward_euler_step(vstate, L, dt);
        //crank_nicolson_step(vstate, L, dt);

        // Intermediate measurements
        if (i % (num_steps / pars.outs) == 0)
        {
            std::cout << "Step " << i << "/" << num_steps-1 << ": ";
            std::cout << "t = " << t << " "; 
            vsol.push_back(ev_state);
            cMat mstate = Vec2Mat(ev_state, mdim, mdim);
            double jz = getJz(mstate, Sz, exp_plus, exp_minus);
            jz_vals.push_back(jz);
            std::cout << "Jz = " << jz << std::endl;
        }

        vstate = ev_state;
        
        // Last measurement
        if (i == num_steps-1)
        {
            std::cout << "Step " << i << "/" << num_steps-1 << ": ";
            std::cout << "t = " << t << " "; 
            vsol.push_back(ev_state);
            cMat mstate = Vec2Mat(ev_state, mdim, mdim);
            double jz = getJz(mstate, Sz, exp_plus, exp_minus);
            jz_vals.push_back(jz);
            std::cout << "Jz = " << jz << std::endl;
        }
    }

    return 0;
}








