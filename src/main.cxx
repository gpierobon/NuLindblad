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

    Clock::time_point start = Clock::now();
    
    bool cacheL = false;
    int    dim  = pars.N + 1;
    double t0   = pars.t_i;

    cMat mstate = DickeState(pars.N);
    cVec vstate = Mat2Vec(mstate);
    int mdim    = sqrt(vstate.size());

    cSpMat L = Lindblad_sparse_pm(&pars);
    if (cacheL) {cache_L(L, "Lindblad.txt"); }
    L.makeCompressed();

    // Precompute Sz, exp_+, exp_- for Jz trace
    std::cout << "Pre-computing exps for time evolution (sparse) ... ";
    SpMat  Sz       (dim, dim);
    cSpMat exp_plus (dim, dim);
    cSpMat exp_minus(dim, dim);
    
    getSz (pars.N, Sz);
    getExp(pars.N,  pars.omega, exp_plus);
    getExp(pars.N, -pars.omega, exp_minus);
    std::cout << " done!" << std::endl;
    
    double rad     = spectral_radius(L);
    double thr     = pars.hthr;
    double h_check = thr/rad;

    std::cout << "Time step set to: dt = " << h_check << std::endl;
    
    double dt = std::min(pars.h, h_check);
    size_t num_steps = static_cast<size_t>((pars.t_f-pars.t_i)/dt);
    std::cout << "\n";
    std::cout << "Running a loop with " << num_steps << " steps ... \n" << std::endl;

    // For the slope
    double t_prev = -1.0;
    double jz_prev = -1.0;

    for (size_t i = 0; i < num_steps; ++i)
    {
        double t = t0 + i * dt;

        evolve(pars.integrator, vstate, L, dt, t);
        //cVec ev_state = krylov_expm(L, vstate, dt, pars.kry); // Krylov case
                                                            
        // Intermediate measurements
        if (i % (num_steps / pars.outs) == 0)
        {
            //cMat mstate = Vec2Mat(ev_state, mdim, mdim);
            cMat mstate = Vec2Mat(vstate, mdim, mdim);
            double jz = getJz(mstate, Sz, exp_plus, exp_minus);
            cache_Jz(pars.N, t, jz);

            double slope = std::numeric_limits<double>::quiet_NaN();
            if (t_prev > 0.0 && std::abs(jz_prev) > 1e-12)
                slope = getJzSlope(t, -jz, t_prev, -jz_prev);
            
            Clock::time_point curr = Clock::now();
            auto dur = curr-start;
            printStatus(i, num_steps, t, jz, dur, slope);

            t_prev = t;
            jz_prev = jz;
        }

        //vstate = ev_state; // Krylov case
        
        // Last measurement
        if (i == num_steps-1)
        {
            //cMat mstate = Vec2Mat(ev_state, mdim, mdim);
            cMat mstate = Vec2Mat(vstate, mdim, mdim);
            double jz = getJz(mstate, Sz, exp_plus, exp_minus);
            cache_Jz(pars.N, t, jz);
            
            double slope = std::numeric_limits<double>::quiet_NaN();
            if (t_prev > 0.0 && std::abs(jz_prev) > 1e-12)
                slope = getJzSlope(t, -jz, t_prev, -jz_prev);

            Clock::time_point curr = Clock::now();
            auto dur = curr-start;
            printStatus(i, num_steps, t, jz, dur, slope);
        }
    }

    Clock::time_point end = Clock::now();
    auto wall = end-start;
    std::string time_str = formatDuration(wall);
    std::cout << "\nFinished: it took " << time_str << std::endl;  
    
    return 0;
}








