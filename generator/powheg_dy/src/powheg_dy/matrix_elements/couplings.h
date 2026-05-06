#pragma once 

#include "powheg_dy/base.h"
#include "powheg_dy/config.h"

namespace powheg_dy
{
    struct _ZHelicityCouplings
    {
        double minus = 0.0;     // helicity -1
        double plus  = 0.0;     // helicity +1

        double operator[](int hel) const
        {
            return (hel == -1) ? minus : plus;
        }
    };

    _ZHelicityCouplings _zCouplings(const Config& config, bool upType, double charge);
    std::complex<double> _photonPropagator(double q2);
    std::complex<double> _zPropagator(const Config& config, double q2);
}