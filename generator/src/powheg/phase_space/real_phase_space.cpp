#include "real_phase_space.h"

#include "powheg/process.h"

namespace powheg
{
namespace
{   
    static constexpr double ALLOWED_REL_MISMATCH = 1.0e-10;

    struct IncomingFractions
    {
        double x1;
        double x2;
    };

    IncomingFractions computeRealFractions(
        const BornPhSpPt& born,
        const RadiationVariables& rad
    )
    {
        const double a = 2.0 - rad.xi * (1.0 - rad.y);
        const double b = 2.0 - rad.xi * (1.0 + rad.y);

        const double sqrt1minusXi = std::sqrt(1.0 - rad.xi);

        // eq. (5.11) of the FKS inverse construction.
        const double x1 = born.x1Bar / sqrt1minusXi * std::sqrt(a / b);
        const double x2 = born.x2Bar / sqrt1minusXi * std::sqrt(b / a);

        return { x1, x2 };
    }

    double radJacobian(
        const RealPhSpPt& real,
        const RadiationVariables& rad
    ) 
    {
        // eq. (5.10) in the paper
        return real.sHatReal
            / std::pow(4.0 * PI, 3)
            * rad.xi
            / (1.0 - rad.xi);
    }

} // anonymous namespace 

    RealPhSpPt ISRRealPhaseSpace::reconstruct(
        const BornPhSpPt& born, 
        const RadiationVariables& rad
    ) const
    {
        RealPhSpPt real;

        const auto [x1, x2] = computeRealFractions(born, rad);
        real.x1 = x1;
        real.x2 = x2;

        real.p1In = { 0.5 * x1 * m_config.SQRT_S, 0.0, 0.0,  0.5 * x1 * m_config.SQRT_S };
        real.p2In = { 0.5 * x2 * m_config.SQRT_S, 0.0, 0.0, -0.5 * x2 * m_config.SQRT_S };
        
        const FourVector totalIn = real.p1In + real.p2In;
        real.sHatReal = totalIn.square();
        
        real.radJacobian = radJacobian(real, rad);

        // eq. (5.1) in the paper
        const double eRadCM = 0.5 * std::sqrt(real.sHatReal) * rad.xi;
        const double sinTh = std::sqrt(1.0 - rad.y * rad.y);

        const FourVector pRadiatedCM = {
            eRadCM,
            eRadCM * sinTh * sin(rad.phi),
            eRadCM * sinTh * cos(rad.phi),
            eRadCM * rad.y
        };
        
        const FourVector pRadiated = pRadiatedCM.boost(totalIn.getBeta());
        const FourVector pBoson = totalIn - pRadiated;

        // construct a longitudinal boost B_L such that B_L k_Boson has zero longitudinal momentum
        const double realRapidity = pBoson.rapidity();
        const ThreeVector longBoost = { 0.0, 0.0, -tanh(realRapidity) };

        // construct a transverse boost B_T such that B_T B_L k_Boson has zero transverse momentum
        const FourVector pBosonBoosted = pBoson.boost(longBoost);
        const ThreeVector transvBoost = { pBosonBoosted.x / pBosonBoosted.e, 
                pBosonBoosted.y / pBosonBoosted.e, 0.0 };

        FourVector totalOut = { };
        for (const FourVector& pBorn : born.pOut)
        {
            // eq. (5.16) in the paper
            const FourVector pReal = pBorn.boost(longBoost).boost(transvBoost).boost(-longBoost);
            
            totalOut += pReal;
            real.pOut.push_back(pReal);
        }

        totalOut += pRadiated;
        real.pOut.push_back(pRadiated);

        FourVector mismatch = totalIn - totalOut;
        const double tolerance = ALLOWED_REL_MISMATCH * std::sqrt(real.sHatReal);

        powheg_assert(std::abs(mismatch.e) < tolerance, "Energy not conserved.");
        powheg_assert(std::abs(mismatch.x) < tolerance, "Momentum not conserved.");
        powheg_assert(std::abs(mismatch.y) < tolerance, "Momentum not conserved.");
        powheg_assert(std::abs(mismatch.z) < tolerance, "Momentum not conserved.");

        return real;
    }

    double ISRRealPhaseSpace::xiMax(
        const BornPhSpPt& born,
        double y
    ) const
    {
        const double xb1sq = born.x1Bar * born.x1Bar;
        const double xb2sq = born.x2Bar * born.x2Bar;

        const double onePlusY  = 1.0 + y;
        const double oneMinusY = 1.0 - y;

        const double denom1 =
            std::sqrt(
                std::pow(1.0 + xb1sq, 2) * oneMinusY * oneMinusY
                + 16.0 * y * xb1sq
            )
            + oneMinusY * (1.0 - xb1sq);

        const double denom2 =
            std::sqrt(
                std::pow(1.0 + xb2sq, 2) * onePlusY * onePlusY
                - 16.0 * y * xb2sq
            )
            + onePlusY * (1.0 - xb2sq);

        const double bound1 = 2.0 * onePlusY * xb1sq / denom1;
        const double bound2 = 2.0 * oneMinusY * xb2sq / denom2;

        // eq. (5.13) in the paper
        return 1.0 - std::max(bound1, bound2);
    }

} // namespace powheg
