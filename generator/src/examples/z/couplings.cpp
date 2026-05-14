#include "couplings.h"

#include "powheg_dy/config/config.h"

namespace powheg_dy
{
    ZHelicityCouplings zCouplings(const Config& config, bool upType, double charge)
    {
        int sign = upType ? 1 : -1;

        const double v = 0.5 * sign - 2.0 * charge * config.S_W_SQ;
        const double a = 0.5 * sign;

        ZHelicityCouplings out;
        out.minus = v + a;      // Fortran Zcoup(-1)
        out.plus  = v - a;      // Fortran Zcoup(+1)
        return out;
    }

    std::complex<double> photonPropagator(double q2)
    {
        return 1.0 / q2;
    }

    std::complex<double> zPropagator(const Config& config, double q2)
    {
        const double mZ  = config.M_Z;
        const double gZ  = config.GAMMA_Z;
        const double mZ2 = mZ * mZ;

        // no running width
        return 1.0 / std::complex<double>(q2 - mZ2, mZ * gZ);
    }

} // namespace powheg_dy
