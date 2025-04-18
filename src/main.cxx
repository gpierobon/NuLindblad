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

    cMat mstate = initState(&pars);
    cVec vstate = Mat2Vec(mstate);
    int mdim    = sqrt(vstate.size());

    cSpMat L = Lindblad_sparse_pm(&pars);

    std::cout << "Pre-computing exps for time evolution (sparse) ... ";
    SpMat  Sz       (dim, dim);
    cSpMat Splus    (dim, dim);
    cSpMat Sminus   (dim, dim);
    cSpMat exp_plus (dim, dim);
    cSpMat exp_minus(dim, dim);
    
    getSz (pars.N,  Sz);
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

    size_t n_bins = pars.outs; 
    std::vector<size_t> log_bins = generateLogBins(t0, dt, num_steps, n_bins);

    double t_prev = -1.0;
    double jz_prev = -1.0;
    double slope;
    double slope_prev = std::numeric_limits<double>::quiet_NaN();
    const double slope_tol = 0.025;
    bool check_pert   = true;
    bool reached_pert = false;

    for (size_t i = 0; i < num_steps; ++i)
    {
        double t = t0 + i * dt;

        evolve(pars.integrator, vstate, L, dt, t);
                                                            
        if (std::binary_search(log_bins.begin(), log_bins.end(), i))
        {
            cMat mstate = Vec2Mat(vstate, mdim, mdim);
            double jz = getJz(mstate, Sz, exp_plus, exp_minus);
            cache_Jz(&pars, t, jz);

            slope = std::numeric_limits<double>::quiet_NaN();
            if (t_prev > 0.0 && std::abs(jz_prev) > 1e-12)
                slope = getJzSlope(t, -jz, t_prev, -jz_prev);

            if (check_pert)
            {
                if (!std::isnan(slope_prev) && !std::isnan(slope))
                {
                    if ( (slope_prev < (1.0 - slope_tol) && slope >= (1.0 - slope_tol)) ||
                         (slope_prev > (1.0 + slope_tol) && slope <= (1.0 + slope_tol)) )
                    {
                        check_pert = false;
                        reached_pert = true;
                        std::cout << "\n ---------------------------------------------------";
                        std::cout << "\n Reached perturbative solution at t = " 
                                  << std::setprecision(5) << t << std::endl;
                        std::cout << " ---------------------------------------------------\n\n";
                    }
                }
            }

            if (reached_pert && slope < pars.sthr)
            {
                std::cout << "\nSlope reached " << std::setprecision(3) 
                          << pars.sthr << ", closing the loop!\n";
                break;
            }
            
            Clock::time_point curr = Clock::now();
            auto dur = curr-start;
            printStatus(i, num_steps, t, jz, dur, slope);

            t_prev = t;
            jz_prev = jz;
            slope_prev = slope;
        }
    }

    Clock::time_point end = Clock::now();
    auto wall = end-start;
    std::string time_str = formatDuration(wall);
    std::cout << "\nFinished: it took " << time_str << std::endl;  
    
    return 0;
}








