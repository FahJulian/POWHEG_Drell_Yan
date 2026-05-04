#pragma once 

#include <complex>

namespace powheg_dy
{
    class Process;
    
    struct _ZHelicityCouplings
    {
        double minus = 0.0;     // helicity -1
        double plus  = 0.0;     // helicity +1

        double operator[](int hel) const
        {
            return (hel == -1) ? minus : plus;
        }
    };

    _ZHelicityCouplings _zCouplings(const Process& process, bool upType, double charge);
    std::complex<double> _photonPropagator(double q2);
    std::complex<double> _zPropagator(const Process& process, double q2);
}