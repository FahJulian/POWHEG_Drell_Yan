#include "real_phase_space.h"

#include "powheg_dy/process.h"

namespace powheg_dy
{
namespace
{   
    static constexpr double ALLOWED_MISMATCH = 1.0e-6;

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

    double kt2(
        const BornPhSpPt& born,
        const RadiationVariables& rad
    ) 
    {
        // eq. (7.233) in the paper
        return born.sHat / (4.0 * (1.0 - rad.xi))
                * rad.xi * rad.xi
                * (1.0 - rad.y * rad.y);
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

    void assertRealKinematics(
        const BornPhSpPt& born,
        const RealPhSpPt& real
    ) 
    {
        const FourVector totalIn = real.p1In + real.p2In;
        const FourVector totalOut = real.pLMinus + real.pLPlus + real.pRadiated;
        const FourVector dilepton = real.pLMinus + real.pLPlus;
        const double kt2FromMomentum =  real.pRadiated.pX * real.pRadiated.pX
            + real.pRadiated.pY * real.pRadiated.pY;

        const double totalMomentumMismatch = (totalIn - totalOut).square() / born.sHat;
        const double gluonMassMismatch = real.pRadiated.square() / born.sHat;
        const double bosonMassMismatch = (dilepton.square() - born.sHat) / born.sHat;
        const double ktMismatch = (kt2FromMomentum - real.kt2) / born.sHat;
        const double bosonMomentumMismatch = (dilepton - real.pBoson).square() / born.sHat;

        powheg_assert(abs(totalMomentumMismatch) < ALLOWED_MISMATCH);
        powheg_assert(abs(gluonMassMismatch) < ALLOWED_MISMATCH);
        powheg_assert(abs(bosonMassMismatch) < ALLOWED_MISMATCH);
        powheg_assert(abs(bosonMomentumMismatch) < ALLOWED_MISMATCH);
        powheg_assert(abs(ktMismatch) < ALLOWED_MISMATCH);
    }

} // anonymous namespace 

    RealPhSpPt FKSRealPhaseSpace::reconstruct(
        const BornPhSpPt& born, 
        const RadiationVariables& rad
    ) const
    {
        RealPhSpPt real;
        real.underlyingBorn = born;
        real.rad = rad;

        const auto [x1, x2] = computeRealFractions(born, rad);
        real.x1 = x1;
        real.x2 = x2;

        const double sqrtS = m_config.SQRT_S;
        real.p1In = { 0.5 * x1 * sqrtS, 0.0, 0.0,  0.5 * x1 * sqrtS };
        real.p2In = { 0.5 * x2 * sqrtS, 0.0, 0.0, -0.5 * x2 * sqrtS };
        
        const FourVector realIncoming = real.p1In + real.p2In;
        real.sHatReal = realIncoming.square();
        
        real.kt2 = kt2(born, rad);
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

        real.pRadiated = pRadiatedCM.boost(realIncoming.getBeta());
        real.pBoson = realIncoming - real.pRadiated;

        // construct a longitudinal boost B_L such that B_L k_Boson has zero longitudinal momentum
        const double realRapidity = real.pBoson.rapidity();
        const ThreeVector longBoost = { 0.0, 0.0, -tanh(realRapidity) };

        // construct a transverse boost B_T such that B_T B_L k_Boson has zero transverse momentum
        const FourVector pBosonBoosted = real.pBoson.boost(longBoost);
        const ThreeVector transvBoost = { pBosonBoosted.x / pBosonBoosted.e, 
                pBosonBoosted.y / pBosonBoosted.e, 0.0 };

        // eq. (5.16) in the paper
        real.pLMinus = born.pOut[0].boost(longBoost).boost(transvBoost).boost(-longBoost);  
        real.pLPlus  = born.pOut[1].boost(longBoost).boost(transvBoost).boost(-longBoost);  
        
        assertRealKinematics(born, real);

        return real;
    }

    double FKSRealPhaseSpace::xiMax(
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

} // namespace powheg_dy
