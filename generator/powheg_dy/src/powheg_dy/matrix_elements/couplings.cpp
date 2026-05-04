#include "couplings.h"

#include "powheg_dy/process.h"

namespace powheg_dy
{
    _ZHelicityCouplings _zCouplings(const Process& process, bool upType, double charge)
    {
        int sign = upType ? 1 : -1;

        const double v = 0.5 * sign - 2.0 * charge * process.S_W_SQ();
        const double a = 0.5 * sign;

        _ZHelicityCouplings out;
        out.minus = v + a;      // Fortran Zcoup(-1)
        out.plus  = v - a;      // Fortran Zcoup(+1)
        return out;
    }

    std::complex<double> _photonPropagator(double q2)
    {
        return 1.0 / q2;
    }

    std::complex<double> _zPropagator(const Process& process, double q2)
    {
        const double mZ  = process.M_Z();
        const double gZ  = process.GAMMA_Z();
        const double mZ2 = mZ * mZ;

        // no running width
        return 1.0 / std::complex<double>(q2 - mZ2, mZ * gZ);
    }

} // namespace powheg_dy
