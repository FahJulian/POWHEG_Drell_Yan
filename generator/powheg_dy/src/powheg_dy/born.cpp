#include "matrix_elements.h"

#include "powheg_dy/math.h"
#include "powheg_dy/process.h"

namespace powheg_dy
{
    namespace 
    {
        struct _NeutralCurrentCouplings
        {
            double symmetric = 0.0;        // multiplies (1 + cos^2 theta)
            double forwardBackward = 0.0;  // multiplies 2 cos theta
        };

        _NeutralCurrentCouplings _neutralCurrentCouplings(const Process& process, int flavour, double sHat)
        {
            bool upType = flavour % 2 == 0;

            // Quark charge and axial/vector Z couplings.
            double qQ = upType ? 2.0 / 3.0 : -1.0 / 3.0;
            double aQ = upType ? 0.5 : -0.5;
            double vQ = upType ? 0.5 - (4.0 / 3.0) * process.S_W_SQ()
                              : -0.5 + (2.0 / 3.0) * process.S_W_SQ();

            double dz = (sHat - process.M_Z() * process.M_Z())
                      * (sHat - process.M_Z() * process.M_Z())
                      + process.M_Z() * process.M_Z() * process.GAMMA_Z() * process.GAMMA_Z();

            double reChi = process.KAPPA() * sHat * (sHat - process.M_Z() * process.M_Z()) / dz;
            double absChiSq = process.KAPPA() * process.KAPPA() * sHat * sHat / dz;

            _NeutralCurrentCouplings couplings;

            couplings.symmetric = qQ * qQ
                - 2.0 * qQ * process.V_L() * vQ * reChi
                + (process.V_L() * process.V_L() + process.A_L() * process.A_L())
                * (vQ * vQ + aQ * aQ) * absChiSq;

            couplings.forwardBackward = -2.0 * qQ * process.A_L() * aQ * reChi
                + 4.0 * process.V_L() * process.A_L() * vQ * aQ * absChiSq;

            return couplings;
        }

    } // namespace 

    double MatrixElements::bornAngularFactor(const Process& process, int flavour, double sHat, double cosTheta)
    {
        auto couplings = _neutralCurrentCouplings(process, flavour, sHat);
        return couplings.symmetric * (1.0 + cosTheta * cosTheta)
                + 2.0 * couplings.forwardBackward * cosTheta;
    }

} // namespace powheg_dy
