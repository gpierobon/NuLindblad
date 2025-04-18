#include <iostream>
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
    
    int    dim  = pars.N + 1;
    double t0   = pars.t_i;

    cMat mstate = ProductState(pars.N);
    cVec vstate = Mat2Vec(mstate);
    int mdim    = sqrt(vstate.size());

    cSpMat L = Lindblad_sparse_pm(&pars);

    std::cout << "Pre-computing exps for time evolution (sparse) ... ";
    cSpMat Splus    (dim, dim);
    cSpMat Sminus   (dim, dim);
    cSpMat exp_plus (dim, dim);
    cSpMat exp_minus(dim, dim);
    
    getSpm(pars.N,  Splus, Sminus);
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

    for (size_t i = 0; i < num_steps; ++i)
    {
        double t = t0 + i * dt;

        evolve(pars.integrator, vstate, L, dt, t);
                                                            
        if (i % (num_steps / pars.outs) == 0)
        {
            cMat mstate = Vec2Mat(vstate, mdim, mdim);
            double jx = getJx(mstate, Splus, Sminus, exp_plus, exp_minus);
            //cache_Jx(&pars, t, jz); // Implement!
            
            Clock::time_point curr = Clock::now();
            auto dur = curr-start;
            printStatusx(i, num_steps, t, jx, dur);
        }
    }

    Clock::time_point end = Clock::now();
    auto wall = end-start;
    std::string time_str = formatDuration(wall);
    std::cout << "\nFinished: it took " << time_str << std::endl;  
    
    return 0;
}








