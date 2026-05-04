#include "born_phase_space.h"

#include "powheg_dy/math/math.h"
#include "powheg_dy/process.h"

namespace powheg_dy
{
    namespace 
    {
        static constexpr double __ALLOWED_MISMATCH = 1.0e-10;

    } // namespace

    BornPhSpPt BornPhaseSpace::samplePoint(double rands[3]) const
    {
        BornPhSpPt point;

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

        point.mB = sqrt(m2);
        point.sHat = m2;

        // Jacobian dm^2 / dr
        const double jacobianM2 =
            (zhigh - zlow) * zMassWidth / (cosZ * cosZ);

        // Sample boson rapidity uniformly from the kinematically allowed range
        const double yBosonMax = log(m_process.sqrtS() / point.mB);
        point.yB = (2.0 * rands[1] - 1.0) * yBosonMax;

        // Sample cos(theta) from p(c) = 3(1+c^2)/8
        point.cosTh =
            2.0 * sinh(asinh(4.0 * rands[2] - 2.0) / 3.0);

        // Calculate x1Bar and x2Bar
        point.x1Bar = point.mB / m_process.sqrtS()
                * exp(point.yB);

        point.x2Bar = point.mB / m_process.sqrtS()
                * exp(-point.yB);

        const double jacobianY =
            2.0 * log(m_process.sqrtS() / point.mB);

        const double jacobianCosTh =
            8.0 / 3.0 / (1.0 + point.cosTh * point.cosTh);

        const double integratedPhi = 2.0 * PI;

        // See explanation below:
        // If your weight formula is still written for dM, use jacobianM.
        const double jacobianM = jacobianM2 / (2.0 * point.mB);

        point.jacobian =
            jacobianM * jacobianY * jacobianCosTh * integratedPhi;

        return point;
    }

    void BornPhaseSpace::reconstructMomenta(BornPhSpPt& point) const
    {
        const double cosThLeg1 = point.channel.id1 > 0 ? point.cosTh : -point.cosTh;
        
        point.pB = {
            point.mB * cosh(point.yB),
            0.0,
            0.0,
            point.mB * sinh(point.yB)
        };

        point.p1Bar = {
            0.5 * point.x1Bar * m_process.sqrtS(),
            0.0,
            0.0,
            0.5 * point.x1Bar * m_process.sqrtS()
        };

        point.p2Bar = {
            0.5 * point.x2Bar * m_process.sqrtS(),
            0.0,
            0.0,
            -0.5 * point.x2Bar * m_process.sqrtS()
        };

        const double p = point.mB / 2.0;
        const double sinTh = sqrt(1.0 - cosThLeg1 * cosThLeg1);

        const FourVector p1Rest = {
            p,
            p * sinTh * cos(point.phi),
            p * sinTh * sin(point.phi),
            p * cosThLeg1
        };

        const FourVector p2Rest = { p1Rest.e, -p1Rest.getThreeVec() };

        point.pLMinus = p1Rest.boost(point.pB.getBeta());
        point.pLPlus = p2Rest.boost(point.pB.getBeta());

        const FourVector totalIn = point.p1Bar + point.p2Bar;
        const FourVector totalOut = point.pLMinus + point.pLPlus;
        double mismatch = (totalIn - totalOut) * (totalIn - totalOut) / point.sHat;

        assert(abs(mismatch) < __ALLOWED_MISMATCH);
    }

} // namespace powheg_dy