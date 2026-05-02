#include "emission_generator.h"

#include "powheg_dy/math.h"
#include "powheg_dy/rand.h"
#include "powheg_dy/process.h"
#include "powheg_dy/matrix_elements/matrix_elements.h"

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/toms748_solve.hpp>

namespace powheg_dy
{
namespace 
{
    static constexpr double __XI_EPS     = 1.0e-12;
    static constexpr double __Y_EPS      = 1.0e-10;

    static constexpr int __MAX_TRIALS = 100000;

    constexpr double __LAMBDA2_OVEREST_FACTOR = 4.0;
    constexpr double __N_Q = 3.0;
    constexpr int __MAX_TRIALS_V = 5000;

    double __alphaSSudakov(const Process& process, double qSq, int nF)
    {
        double alphaS = process.alphaSOneLoop(qSq, nF);
        double bracket = (67.0 / 18.0 - PI*PI / 6.0) * process.C_A() - 5.0 / 9.0 * nF;

        return alphaS * (1.0 + alphaS / 2.0 / PI * bracket);
    }

    double __regionSign(RadiationRegion region)
    {
        switch (region)
        {
            case RadiationRegion::Plus:  return +1.0;
            case RadiationRegion::Minus: return -1.0;
            case RadiationRegion::None:  return  0.0;
        }

        assert(false);
        return 0.0;
    }

    Emission __makeAcceptedEmission(
        const RealPhSpPt& real,
        RadiationRegion region,
        double exact,
        double upper,
        double ratio
    )
    {
        Emission emission;

        emission.rad = real.rad;
        emission.reg = region;
        emission.kt2 = real.kt2;

        emission.exactDensity = exact;
        emission.upperDensity = upper;
        emission.acceptanceRatio = ratio;

        emission.rejected = false;

        return emission;
    }

} // namespace

    Emission EmissionGenerator::generateEmission(const BornPhSpPt& born) const 
    {
        // Generate one emission for each region and pick the harder one

        const Emission plus = _generateCandidate(born, RadiationRegion::Plus);
        const Emission minus = _generateCandidate(born, RadiationRegion::Minus);

        const double kt2Plus = plus.rejected ? 0.0 : plus.kt2;
        const double kt2Minus = minus.rejected ? 0.0 : minus.kt2;

        return kt2Plus > kt2Minus ? plus : minus;
    }

    Emission EmissionGenerator::_generateCandidate(
        const BornPhSpPt& born,
        RadiationRegion region
    ) const
    {
        double kt2Max = _globalKt2Max(born);

        for (int trial = 0; trial < __MAX_TRIALS && kt2Max > m_process.pt2Cutoff(); ++trial)
        {
            const double kt2Trial = _sampleTrialKt2(born, kt2Max);
            if (kt2Trial < m_process.pt2Cutoff())
                return Emission().reject();
            
            RadiationVariables rad = _sampleTrialRadiation(born, region, kt2Trial);

            const double xiMax = m_realPhaseSpace.xiMax(born, rad.y);
            if (rad.xi <= __XI_EPS || rad.xi >= 1.0 - __XI_EPS || rad.xi >= xiMax)
            {
                kt2Max = kt2Trial;
                continue;
            }

            RealPhSpPt real = m_realPhaseSpace.reconstruct(born, rad);

            const double exact = _exactRadiationDensity(real);
            const double upper = _upperRadiationDensity(real, kt2Trial);

            assert(upper > 0);
            const double ratio = exact / upper;

            assert(ratio <= 1.0);

            if (rand() < ratio)
            {
                return __makeAcceptedEmission(
                    real,
                    region,
                    exact,
                    upper,
                    ratio
                );
            }
            else
                kt2Max = kt2Trial;
        }
        
        return Emission().reject();
    }

    double EmissionGenerator::_globalKt2Max(const BornPhSpPt& born) const
    {
        const double rho = born.sHat / m_process.S();
        return 0.25 * m_process.S() * (1.0 - rho) * (1.0 - rho);
    }

    RadiationVariables EmissionGenerator::_sampleTrialRadiation(
        const BornPhSpPt& born,
        RadiationRegion region,
        double kt2Trial
    ) const
    {
        RadiationVariables rad;

        rad.xi = _sampleTrialXi(born, kt2Trial);
        rad.y = _computeYForRegion(born, kt2Trial, rad.xi, region);
        rad.phi = _sampleTrialPhi();

        return rad;
    }

    double EmissionGenerator::_exactRadiationDensity(const RealPhSpPt& real) const
    {
        return real.radJacobian * MatrixElements::realOverBorn(m_process, real, real.kt2, real.kt2);
    }

    double EmissionGenerator::_upperRadiationDensity(const RealPhSpPt& real, double kt2Trial) const
    {
        const double alphaS = __alphaSSudakov(m_process, kt2Trial, 5);

        return __N_Q
            * alphaS
            / real.rad.xi 
            / (1 - real.rad.y * real.rad.y);
    }

    double EmissionGenerator::_sampleTrialPhi() const
    {
        return rand() * 2.0 * PI;
    }

    double EmissionGenerator::_sampleTrialXi(const BornPhSpPt& born, double pT2) const
    {
        const double b = sqrt(pT2) / born.mB;
        const double a = sqrt(1.0 + b * b);
        const double xiPlus = 1.0 - (a - b) * (a - b);
        const double xiMinus = 1.0 - (a + b) * (a + b);
        const double rho = born.sHat / m_process.S();

        const double delta = xiPlus - xiMinus;
        const double eta0 = sqrt(delta);
        const double etaMax = sqrt((1 - rho) - xiPlus) + sqrt((1 - rho) - xiMinus);

        const double u = rand();
        const double eta = eta0 * pow(etaMax / eta0, u);

        const double bracket = (eta - delta / eta);
        return xiPlus + bracket * bracket / 4.0;
    }

    double EmissionGenerator::_computeYForRegion(const BornPhSpPt& born, double pT2, double xi, RadiationRegion region) const
    {
        const double xiFactor = (1.0 - xi) / xi / xi;
        const double absY = sqrt(1.0 - 4.0 * pT2 / born.sHat * xiFactor);

        return __regionSign(region) * absY;
    }

    // sample a k_T^2 after the procedure describes in Appendix D of ref. [22]
    double EmissionGenerator::_sampleTrialKt2(const BornPhSpPt& born, double ptMax2) const
    {
        assert(ptMax2 >= m_process.pt2Cutoff());

        for (int trial = 1; trial <= __MAX_TRIALS_V; trial++)
        {
            const double logR = log(rand());

            const double uLow  = log(m_process.pt2Cutoff());
            const double uHigh = log(ptMax2);

            const double maxIntegral = _integrateVTildeLog(uLow, uHigh, born.sHat, 5);

            if (-logR > maxIntegral)
                return -1.0;

            auto f = [&](double u) -> double
            {
                return _integrateVTildeLog(u, uHigh, born.sHat, 5) + logR;
            };

            const double fLow  = f(uLow);   // >= 0
            const double fHigh = logR;  // = -target < 0

            boost::uintmax_t maxIter = 100;
            boost::math::tools::eps_tolerance<double> tol(40);

            auto bracket = boost::math::tools::toms748_solve(
                f,
                uLow,
                uHigh,
                fLow,
                fHigh,
                tol,
                maxIter
            );

            const double uRoot = 0.5 * (bracket.first + bracket.second);
            const double trialPt2 = std::exp(uRoot);

            if (trialPt2 < m_process.pt2Cutoff())
                return -1.0;

            const double accRatio = _VExact(trialPt2, born.sHat, 5) / _VTildeLog(uRoot, born.sHat, 5);
            assert(accRatio <= 1.0);

            if (rand() < accRatio)
                return trialPt2;
            else 
                ptMax2 = trialPt2;
        }

        assert(false); // Max trials exceeded!
        return -1.0;
    }

    double EmissionGenerator::_VExact(double pt2, double sHat, int nF) const 
    {
        const double prefactor = PI * __N_Q;      // One emission per region (y pos, y neg) generated, so only half the real prefactor
        const double alphaS = __alphaSSudakov(m_process, pt2, nF);

        const double ratio = pt2 / sHat;
        const double a = sqrt(1 + ratio);
        const double b = sqrt(ratio);

        const double xPlus = (a + b) * (a + b);
        const double xMinus = (a - b) * (a - b);
        const double rho = sHat / m_process.S();

        const double sqrtPlus = sqrt(xPlus - rho);
        const double sqrtMinus = sqrt(xMinus - rho);

        const double logFactor = log((sqrtPlus + sqrtMinus) / (sqrtPlus - sqrtMinus));

        return prefactor * alphaS * logFactor;
    }

    double EmissionGenerator::_VTildeLog(double logPt2, double sHat, int nF) const 
    {
        const double beta0 = (33.0 - 2.0 * nF) / 12.0 / PI;
        const double prefactor = PI * __N_Q / 2.0 / beta0;          // Only half the real prefactor
        const double logLambda2 = log(__LAMBDA2_OVEREST_FACTOR * m_process.LAMBDA_SQ_QCD());

        const double oneMinusRho = 1.0 - sHat / m_process.S();
        const double ktmaxGlobal2 = m_process.S() / 4.0 * oneMinusRho * oneMinusRho;
        const double logQ2 = log(ktmaxGlobal2 + sHat);

        return prefactor 
            / (logPt2 - logLambda2)
            * (logQ2 - logPt2);
    }

    double EmissionGenerator::_integrateVTildeLog(double logPt2, double logKtmax2, double sHat, int nF) const 
    {
        const double beta0 = (33.0 - 2.0 * nF) / 12.0 / PI;
        const double prefactor = PI * __N_Q / 2.0 / beta0;          // Only half the real prefactor 
        const double logLambda2 = log(__LAMBDA2_OVEREST_FACTOR * m_process.LAMBDA_SQ_QCD());

        const double oneMinusRho = 1.0 - sHat / m_process.S();
        const double ktmaxGlobal2 = m_process.S() / 4.0 * oneMinusRho * oneMinusRho;
        const double logQ2 = log(ktmaxGlobal2 + sHat);

        return prefactor * (
            (logQ2 - logLambda2) 
            * log((logKtmax2 - logLambda2) / (logPt2 - logLambda2))
            - logKtmax2 + logPt2
        );
    }

} // namespace powheg_dy
