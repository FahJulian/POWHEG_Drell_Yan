#include "phase_space.h"

#include "powheg_dy/math.h"
#include "powheg_dy/process.h"

namespace powheg_dy
{
    PhaseSpacePoint PhaseSpaceSampler::samplePoint(double rands[3]) const
    {
        PhaseSpacePoint point;

        const double zMass = m_process.zMass();
        const double zWidth = m_process.zWidth();

        const double zMass2 = zMass * zMass;
        const double zMassWidth = zMass * zWidth;

        const double m2Low = m_process.mMin() * m_process.mMin();
        const double m2High = m_process.mMax() * m_process.mMax();

        // Breit-Wigner-like sampling of m^2 around the Z peak
        const double zlow = atan((m2Low - zMass2) / zMassWidth);
        const double zhigh = atan((m2High - zMass2) / zMassWidth);

        const double z = zlow + (zhigh - zlow) * rands[0];

        const double cosZ = cos(z);

        const double m2 = zMassWidth * tan(z) + zMass2;

        point.mBoson = sqrt(m2);
        point.sHat = m2;

        // Jacobian dm^2 / dr
        const double jacobianM2 =
            (zhigh - zlow) * zMassWidth / (cosZ * cosZ);

        // Sample boson rapidity uniformly from the kinematically allowed range
        const double yBosonMax = log(m_process.sqrtS() / point.mBoson);
        point.yBoson = (2.0 * rands[1] - 1.0) * yBosonMax;

        // Sample cos(theta) from p(c) = 3(1+c^2)/8
        point.cosTh =
            2.0 * sinh(asinh(4.0 * rands[2] - 2.0) / 3.0);

        // Calculate x1 and x2
        point.x1 = point.mBoson / m_process.sqrtS()
                * exp(point.yBoson);

        point.x2 = point.mBoson / m_process.sqrtS()
                * exp(-point.yBoson);

        const double jacobianY =
            2.0 * log(m_process.sqrtS() / point.mBoson);

        const double jacobianCosTh =
            8.0 / 3.0 / (1.0 + point.cosTh * point.cosTh);

        const double integratedPhi = 2.0 * PI;

        // See explanation below:
        // If your weight formula is still written for dM, use jacobianM.
        const double jacobianM = jacobianM2 / (2.0 * point.mBoson);

        point.invSamplingFact =
            jacobianM * jacobianY * jacobianCosTh * integratedPhi;

        return point;
    }

} // namespace powheg_dy