#pragma once 

#include "powheg_dy/config/config.h"

namespace powheg_dy
{
    struct ZHelicityCouplings
    {
        double minus = 0.0;     // helicity -1
        double plus  = 0.0;     // helicity +1

        double operator[](int hel) const
        {
            return (hel == -1) ? minus : plus;
        }
    };

    struct DrellYanConfig : public Config 
    {
        ////////// mandatory process params ///////////
        double M_MIN;
        double M_MAX;

        //////////// independent constants ////////////
        double M_Z = 91.1876;
        double M_W = 80.398;
        double GAMMA_Z = 2.4952;
        double GAMMA_W = 2.141;

        double ALPHA_EW = 1.0 / 128.89;

        ///////////// dependent constants /////////////
        double M_Z_SQ;
        double M_W_SQ;

        double S_W_SQ;
        double C_W_SQ;
        double S_W;
        double C_W;
        
        double E_SQ;

        std::complex<double> photonPropagator(const double q2) const;
        std::complex<double> zPropagator(const double q2) const;
        ZHelicityCouplings zCouplings(const bool upType, const double charge) const;
    };

} // namespace powheg_dy
