#include "emission_sampler.h"

#include "powheg_dy/math/rand.h"

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/toms748_solve.hpp>

namespace powheg_dy
{
namespace 
{
    constexpr int MAX_TRIALS_V = 5000;

    constexpr double LAMBDA_SQ_OVEREST_FACTOR = 4.0;        // like in powheg
    constexpr double N_Q = 10;

} // anonymous namespace

    double EmissionSampler::globalKt2Max(const BornPhSpPt& born) const
    {
        const double sumXbar = born.x1Bar + born.x2Bar;
        return born.sHat * (1 - born.x2Bar * born.x2Bar) * (1 - born.x1Bar * born.x1Bar)
            / (sumXbar * sumXbar);
    }

    RadiationVariables EmissionSampler::sampleTrialRadiation(
        const BornPhSpPt& born,
        double kt2Trial
    ) const
    {
        RadiationVariables rad;
        
        rad.xi = sampleTrialXi(born, kt2Trial);
        rad.y = sampleY(born, kt2Trial, rad.xi);
        rad.phi = sampleTrialPhi();
        rad.xiMax = m_realPhaseSpace->xiMax(born, rad.y);

        return rad;
    }

    double EmissionSampler::upperRadiationDensity(const RadiationVariables& rad, double kt2Trial) const
    {
        const double alphaS = m_config.alphaSCMW(kt2Trial);

        return N_Q * alphaS
            / rad.xi 
            / (1 - rad.y * rad.y);
    }

    double EmissionSampler::sampleTrialPhi() const
    {
        return rand() * 2.0 * PI;
    }

    double EmissionSampler::sampleTrialXi(const BornPhSpPt& born, double pT2) const
    {
        const double r = pT2 / born.sHat;
        const double a = std::sqrt(1.0 + r);
        const double b = std::sqrt(r);

        const double xPlus  = (a + b) * (a + b);
        const double xMinus = (a - b) * (a - b);

        // This must match vExactOverVTilde
        const double xMin = std::min(born.x1Bar, born.x2Bar) / (2.0 * a);

        powheg_assert(xMin < xMinus);

        const double delta = xPlus - xMinus;

        const double eta0   = std::sqrt(delta);
        const double etaMax = std::sqrt(xPlus - xMin) + std::sqrt(xMinus - xMin);

        const double u = rand();
        const double eta = eta0 * std::pow(etaMax / eta0, u);

        const double bracket = eta - delta / eta;

        const double x = xMinus - 0.25 * bracket * bracket;
        const double xi = 1.0 - x;

        return xi;
    }

    double EmissionSampler::sampleY(const BornPhSpPt& born, double pT2, double xi) const
    {
        const double xiFactor = (1.0 - xi) / xi / xi;
        const double absY = std::sqrt(std::abs(1.0 - 4.0 * pT2 / born.sHat * xiFactor));
        const double sign = rand() > 0.5 ? 1.0 : -1.0;

        return sign * absY;
    }

    // sample a k_T^2 after the procedure describes in Appendix D of ref. [22]
    double EmissionSampler::sampleTrialKt2(const BornPhSpPt& born, double ptMax2, double& logR) const
    {
        powheg_assert(ptMax2 >= m_config.PT_SQ_CUTOFF);

        const int nF = 5;
        const double LAMBDA2 = LAMBDA_SQ_OVEREST_FACTOR * m_config.LAMBDA_MSB_5_SQ;
        const double BETA0 = (11.0 * m_config.C_A - 4.0 * m_config.T_F * nF) / 12.0 / PI;

        // logR = 0.0;
        for (int trial = 1; trial <= MAX_TRIALS_V; trial++)
        {
            logR = std::log(rand());

            const double maxIntegral = integrateVTilde(m_config.PT_SQ_CUTOFF, ptMax2, born.sHat, LAMBDA2, N_Q, BETA0);

            if (-logR > maxIntegral)
                return -1.0;

            auto f = [&](double pt2) -> double
            {
                return integrateVTilde(pt2, ptMax2, born.sHat, LAMBDA2, N_Q, BETA0) + logR;
            };

            const double fLow  = f(m_config.PT_SQ_CUTOFF);   // >= 0
            const double fHigh = logR;  // = -target < 0

            const boost::uintmax_t MAX_ITER = 100;
            boost::uintmax_t maxIter = MAX_ITER;
            boost::math::tools::eps_tolerance<double> tol(40);

            auto bracket = boost::math::tools::toms748_solve(
                f,
                m_config.PT_SQ_CUTOFF,
                ptMax2,
                fLow,
                fHigh,
                tol,
                maxIter
            );

            powheg_assert(maxIter <= MAX_ITER);    // If maxIter > MAX_ITER, the solver failed

            const double trialPt2 = 0.5 * (bracket.first + bracket.second);

            if (trialPt2 < m_config.PT_SQ_CUTOFF)
                return -1.0;
            
            const double alphaSCorr = m_config.alphaSCMW(trialPt2) / m_config.alphaS0customLambda(trialPt2, 5, LAMBDA2);
            const double accRatio = alphaSCorr * vExactOverVTilde(born, trialPt2);
            powheg_assert(accRatio <= 1.0, "Acceptance Ratio" << accRatio << "exceeds one in emission generation");

            if (rand() < accRatio)
                return trialPt2;
            else 
                ptMax2 = trialPt2;
        }

        powheg_assert(false); // Max trials exceeded!
        return -1.0;
    }

    double EmissionSampler::vExactOverVTilde(const BornPhSpPt& born, double pt2) const
    {
        const double ratio = pt2 / born.sHat;
        
        const double a = std::sqrt(1.0 + ratio);
        const double b = std::sqrt(ratio);

        const double xPlus = (a + b) * (a + b);
        const double xMinus = (a - b) * (a - b);

        const double xMin = std::min(born.x1Bar, born.x2Bar) / 2.0 / a;

        const double rootPlus = std::sqrt(xPlus - xMin);
        const double rootMinus = std::sqrt(xMinus - xMin);

        const double prefactor = std::log((rootPlus + rootMinus) / (rootPlus - rootMinus));

        if (pt2 < born.sHat)
            return 2.0 * prefactor / std::log(2.0 / ratio);
        else 
            return 2.0 * prefactor / std::log(2.0);
    }

    double EmissionSampler::integrateVTilde(
        double pt2,
        double kt2max,
        double sBorn,
        double lambda2,
        double nQ,
        double beta0
    ) const
    {
        const double Lpt  = std::log(pt2 / lambda2);
        const double Lmax = std::log(kt2max / lambda2);
        const double Ls   = std::log(sBorn / lambda2);

        const double prefactor = M_PI * nQ / beta0;

        double integral = 0.0;
        if (pt2 < sBorn) 
        {
            if (sBorn < kt2max) 
            {
                integral = std::log(2.0 * sBorn / lambda2) * std::log(Ls / Lpt)
                    - std::log(sBorn / pt2)
                    + std::log(2.0) * std::log(Lmax / Ls);
            } 
            else 
            {
                integral = std::log(2.0 * sBorn / lambda2) * std::log(Lmax / Lpt)
                    - std::log(kt2max / pt2);
            }
        } 
        else 
        {
            integral = std::log(2.0) * std::log(Lmax / Lpt);
        }

        return prefactor * integral;
    }

} // namespace powheg_dy
