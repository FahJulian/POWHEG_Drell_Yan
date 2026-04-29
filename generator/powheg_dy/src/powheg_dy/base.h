#pragma once

#include <math.h>
#include <stdexcept>

#define ASSERT(cond, msg) if (!(cond)) throw std::runtime_error(msg);

namespace powheg_dy
{
    namespace Math
    {
        static constexpr double PI = 3.14159265358979323846;

    } // namespace Math

    namespace Physics
    {
        static constexpr double ALPHA = 1.0 / 137.035999084;
        static constexpr double NC = 3.0;
        static constexpr double GEV2_TO_PB = 0.389379338e9;
        static constexpr double S_W_SQ  = 0.23126;
        static constexpr double C_W_SQ  = 1.0 - S_W_SQ;
        static constexpr double M_Z     = 91.1876;   
        static constexpr double GAMMA_Z = 2.4952;    
        static constexpr double KAPPA   = 1.0 / (4.0 * S_W_SQ * C_W_SQ);
        static constexpr double LAMBDA_SQ_QCD = 0.2*0.2;
        static constexpr double C_F = 4.0 / 3.0;

        // charged lepton axial and vector couplings
        static constexpr double A_L = -0.5;
        static constexpr double V_L = -0.5 + 2.0 * S_W_SQ;

        inline double alphaSOneLoop(double qSq, int nF)
        {
            ASSERT (qSq > LAMBDA_SQ_QCD, "Invalid Scale");

            double beta0 = 11.0 - 2.0 / 3.0 * nF;
            return 4.0 * Math::PI / beta0 / std::log(qSq / LAMBDA_SQ_QCD);
        }

    } // namespace Physics

} // namespace powheg_dy
