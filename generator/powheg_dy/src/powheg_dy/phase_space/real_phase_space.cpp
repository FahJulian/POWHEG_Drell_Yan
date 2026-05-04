#include "real_phase_space.h"

#include "powheg_dy/math/math.h"
#include "powheg_dy/process.h"

namespace powheg_dy
{
    namespace
    {   
        static constexpr double __ALLOWED_MISMATCH = 1.0e-9;

        struct __IncomingFractions
        {
            double x1;
            double x2;
        };

        __IncomingFractions __computeRealFractions(
            const BornPhSpPt& born,
            const RadiationVariables& rad
        )
        {
            const double a = 2.0 - rad.xi * (1.0 - rad.y);
            const double b = 2.0 - rad.xi * (1.0 + rad.y);

            const double sqrt1minusXi = sqrt(1.0 - rad.xi);

            // eq. (5.11) of the FKS inverse construction.
            const double x1 = born.x1Bar / sqrt1minusXi * sqrt(a / b);
            const double x2 = born.x2Bar / sqrt1minusXi * sqrt(b / a);

            return { x1, x2 };
        }

        double __kt2(
            const BornPhSpPt& born,
            const RadiationVariables& rad
        ) 
        {
            // eq. (7.233) in the paper
            return born.sHat / (4.0 * (1.0 - rad.xi))
                 * rad.xi * rad.xi
                 * (1.0 - rad.y * rad.y);
        }

        double __radJacobian(
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

        void __assertRealKinematics(
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
            const double rapidityMismatch = real.pBoson.rapidity() - born.yB;

            assert(abs(totalMomentumMismatch) < __ALLOWED_MISMATCH);
            assert(abs(gluonMassMismatch) < __ALLOWED_MISMATCH);
            assert(abs(bosonMassMismatch) < __ALLOWED_MISMATCH);
            assert(abs(bosonMomentumMismatch) < __ALLOWED_MISMATCH);
            assert(abs(ktMismatch) < 1.0e-8);
            assert(abs(rapidityMismatch) < 1.0e-8);
        }

    } // namespace 

    RealPhSpPt FKSRealPhaseSpace::reconstruct(
        const BornPhSpPt& born, 
        const RadiationVariables& rad
    ) const
    {
        RealPhSpPt real;
        real.underlyingBorn = born;
        real.rad = rad;

        const auto [x1, x2] = __computeRealFractions(born, rad);
        real.x1 = x1;
        real.x2 = x2;

        const double sqrtS = m_process.sqrtS();
        real.p1In = { 0.5 * x1 * sqrtS, 0.0, 0.0,  0.5 * x1 * sqrtS };
        real.p2In = { 0.5 * x2 * sqrtS, 0.0, 0.0, -0.5 * x2 * sqrtS };
        
        const FourVector realIncoming = real.p1In + real.p2In;
        real.sHatReal = realIncoming.square();
        
        real.kt2 = __kt2(born, rad);
        real.radJacobian = __radJacobian(real, rad);

        // eq. (5.1) in the paper
        const double eRadCM = 0.5 * sqrt(real.sHatReal) * rad.xi;
        const double sinTh = sqrt(1.0 - rad.y * rad.y);

        const FourVector pRadiatedCM = {
            eRadCM,
            eRadCM * sinTh * cos(rad.phi),
            eRadCM * sinTh * sin(rad.phi),
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
        real.pLMinus = born.pLMinus.boost(longBoost).boost(transvBoost).boost(-longBoost);  
        real.pLPlus  = born.pLPlus .boost(longBoost).boost(transvBoost).boost(-longBoost);  
        
        __assertRealKinematics(born, real);

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
        return 1.0 - max(bound1, bound2);
    }

} // namespace powheg_dy
