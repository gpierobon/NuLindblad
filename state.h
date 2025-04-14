#ifndef _STATE_H__
#define _STATE_H__
#include <Eigen/Dense>
#include <types.h>

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

//cMat initState(Params* pars)
//{
//    int N  = pars->N;
//    bool p = pars->pstate;
//    cMat state;
//
//    if (p){
//        std::cout << "Initialising product state (to implement) ..." << std::endl;
//        // state = ProductState(N);  // Still needs implementation
//    }
//    else{
//        std::cout << "Initialising Dicke for " << N << " spins ..." << std::endl;
//        state = DickeState(N);
//    }
//    return state;
//}

#endif
