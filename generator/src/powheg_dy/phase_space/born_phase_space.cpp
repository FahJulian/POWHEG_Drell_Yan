#include "born_phase_space.h"

#include "powheg_dy/process.h"

namespace powheg_dy
{
namespace 
{
    static constexpr double ALLOWED_MISMATCH = 1.0e-10;

} // anonymous namespace

    BornPhSpPt BornPhaseSpace::samplePoint(double rands[3]) const
    {
        BornPhSpPt born;

        const double zMass = m_config.M_Z;
        const double zWidth = m_config.GAMMA_Z;

        const double zMass2 = zMass * zMass;
        const double zMassWidth = zMass * zWidth;

        const double m2Low = m_config.M_MIN * m_config.M_MIN;
        const double m2High = m_config.M_MAX * m_config.M_MAX;

        // Breit-Wigner-like sampling of m^2 around the Z peak
        const double zlow = std::atan((m2Low - zMass2) / zMassWidth);
        const double zhigh = std::atan((m2High - zMass2) / zMassWidth);

        const double z = zlow + (zhigh - zlow) * rands[0];

        const double cosZ = std::cos(z);

        const double m2 = zMassWidth * std::tan(z) + zMass2;

        born.mBoson = std::sqrt(m2);
        born.sHat = m2;

        // Jacobian dm^2 / dr
        const double jacobianM2 =
            (zhigh - zlow) * zMassWidth / (cosZ * cosZ);

        // Sample boson rapidity uniformly from the kinematically allowed range
        const double yBosonMax = std::log(m_config.SQRT_S / born.mBoson);
        born.yB = (2.0 * rands[1] - 1.0) * yBosonMax;

        // Sample cos(theta) from p(c) = 3(1+c^2)/8
        born.cosTh =
            2.0 * std::sinh(std::asinh(4.0 * rands[2] - 2.0) / 3.0);

        // Calculate x1Bar and x2Bar
        born.x1Bar = born.mBoson / m_config.SQRT_S
                * std::exp(born.yB);

        born.x2Bar = born.mBoson / m_config.SQRT_S
                * std::exp(-born.yB);

        const double jacobianY =
            2.0 * std::log(m_config.SQRT_S / born.mBoson);

        const double jacobianCosTh =
            8.0 / 3.0 / (1.0 + born.cosTh * born.cosTh);

        const double integratedPhi = 2.0 * PI;

        born.jacobianOld =
            jacobianM2 * jacobianY * jacobianCosTh * integratedPhi;
        
        born.jacobian = 
            1.0 / (32.0 * PI * PI)              // dPhi = 1/32\pi^2 d\cos\theta d\phi
            * 1.0 / m_config.S                  // Change of variables from (x1, x2) to (sHat, yB)
            * jacobianM2 * jacobianY * jacobianCosTh * integratedPhi;   // Sampling jacobians
        
            born.jacobianOld * 1.0 / (32.0 * PI * PI);

        return born;
    }

    void BornPhaseSpace::reconstructMomenta(BornPhSpPt& born) const
    {
        const double cosThLeg1 = born.channel.id1 > 0 ? born.cosTh : -born.cosTh;
        
        born.pBoson = {
            born.mBoson * std::cosh(born.yB),
            0.0,
            0.0,
            born.mBoson * std::sinh(born.yB)
        };

        born.p1Bar = {
            0.5 * born.x1Bar * m_config.SQRT_S,
            0.0,
            0.0,
            0.5 * born.x1Bar * m_config.SQRT_S
        };

        born.p2Bar = {
            0.5 * born.x2Bar * m_config.SQRT_S,
            0.0,
            0.0,
            -0.5 * born.x2Bar * m_config.SQRT_S
        };

        const double p = born.mBoson / 2.0;
        const double sinTh = std::sqrt(1.0 - cosThLeg1 * cosThLeg1);

        const FourVector p1Rest = {
            p,
            p * sinTh * std::cos(born.phi),
            p * sinTh * std::sin(born.phi),
            p * cosThLeg1
        };

        const FourVector p2Rest = { p1Rest.e, -p1Rest.getThreeVec() };

        born.pLMinus = p1Rest.boost(born.pBoson.getBeta());
        born.pLPlus = p2Rest.boost(born.pBoson.getBeta());

        const FourVector totalIn = born.p1Bar + born.p2Bar;
        const FourVector totalOut = born.pLMinus + born.pLPlus;
        double mismatch = dot(totalIn - totalOut, totalIn - totalOut) / born.sHat;

        powheg_assert(abs(mismatch) < ALLOWED_MISMATCH);
    }

} // namespace powheg_dy