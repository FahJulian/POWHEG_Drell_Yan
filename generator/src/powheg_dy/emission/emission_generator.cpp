#include "emission_generator.h"

#include "powheg_dy/process.h"
#include "powheg_dy/alpha_s.h"
#include "powheg_dy/math/rand.h"
#include "powheg_dy/matrix_elements.h"

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/toms748_solve.hpp>

namespace powheg_dy
{
namespace 
{
    static constexpr double EPS_XI     = 1.0e-12;
    static constexpr double EPS_Y      = 1.0e-10;

    static constexpr int MAX_TRIALS = 100000;
    constexpr int MAX_TRIALS_V = 5000;

    constexpr double LAMBDA_SQ_OVEREST_FACTOR = 4.0;        // like in powheg
    constexpr double N_Q = 10;

    RealChannel chooseChannel(const RealOverBornContributions& contributions)
    {
        double u = rand(0.0, contributions.total);

        for (const auto& [channel, contribution] : contributions.channels)
        {
            if (u < contribution)
                return channel;
            
            u -= contribution;
        }

        throw std::runtime_error("Invalid RealOverBornContributions.");
        return { };
    }

} // anonymous namespace

    Emission EmissionGenerator::generateEmission(
        const BornPhSpPt& born,
        int region
    ) const
    {
        if (region != 1)
            assert(false);      // Only ISR Radiation is implemented

        return generateISREmission(born);
    }

    Emission EmissionGenerator::generateISREmission(const BornPhSpPt& born) const
    {
        double kt2Max = globalKt2Max(born);
        double amp2Born = m_process.bornAmp2(born);

        double logR = 0;
        for (int trial = 0; trial < MAX_TRIALS && kt2Max > m_config.PT_SQ_CUTOFF; ++trial)
        {
            const double kt2Trial = sampleTrialKt2(born, kt2Max, logR);
            if (kt2Trial < m_config.PT_SQ_CUTOFF)
                return Emission().reject();
            
            const RadiationVariables rad = sampleTrialRadiation(born, kt2Trial);

            if (rad.xi <= EPS_XI || rad.xi >= 1.0 - EPS_XI || rad.xi >= rad.xiMax)
            {
                kt2Max = kt2Trial;
                continue;
            }

            const double muF2 = kt2Trial;
            const double muR2 = kt2Trial;
            const RealPhSpPt real = m_realPhaseSpace->reconstruct(born, rad);

            const RealOverBornContributions contributions 
                = getRealOverBornContributions(real, born, amp2Born, muF2, muR2);

            const double accRatio = contributions.total / upperRadiationDensity(real, kt2Trial);

            if (accRatio > 1.0)
                Log::warn << "Acceptance ratio " << accRatio << " is greater than one, accepting emission." << std::endl;

            if (rand() < accRatio)
            {
                RealChannel channel = chooseChannel(contributions);

                return {
                    .rad = real.rad,
                    .channel = channel,
                    .kt2 = real.kt2,
                    .rejected = false
                };
            }
            else
                kt2Max = kt2Trial;
        }
        
        return Emission().reject();
    }

    double EmissionGenerator::globalKt2Max(const BornPhSpPt& born) const
    {
        const double sumXbar = born.x1Bar + born.x2Bar;
        return born.sHat * (1 - born.x2Bar * born.x2Bar) * (1 - born.x1Bar * born.x1Bar)
            / (sumXbar * sumXbar);
    }

    RadiationVariables EmissionGenerator::sampleTrialRadiation(
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

    double EmissionGenerator::upperRadiationDensity(const RealPhSpPt& real, double kt2Trial) const
    {
        const double alphaS = alphaSCMW(m_config, kt2Trial);

        return N_Q
            * alphaS
            / real.rad.xi 
            / (1 - real.rad.y * real.rad.y);
    }

    double EmissionGenerator::sampleTrialPhi() const
    {
        return rand() * 2.0 * PI;
    }

    double EmissionGenerator::sampleTrialXi(const BornPhSpPt& born, double pT2) const
    {
        const double r = pT2 / born.sHat;
        const double a = std::sqrt(1.0 + r);
        const double b = std::sqrt(r);

        const double xPlus  = (a + b) * (a + b);
        const double xMinus = (a - b) * (a - b);

        // This must match vExactOverVTilde
        const double xMin = std::min(born.x1Bar, born.x2Bar) / (2.0 * a);

        assert(xMin < xMinus);

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

    double EmissionGenerator::sampleY(const BornPhSpPt& born, double pT2, double xi) const
    {
        const double xiFactor = (1.0 - xi) / xi / xi;
        const double absY = std::sqrt(std::abs(1.0 - 4.0 * pT2 / born.sHat * xiFactor));
        const double sign = rand() > 0.5 ? 1.0 : -1.0;

        return sign * absY;
    }

    // sample a k_T^2 after the procedure describes in Appendix D of ref. [22]
    double EmissionGenerator::sampleTrialKt2(const BornPhSpPt& born, double ptMax2, double& logR) const
    {
        assert(ptMax2 >= m_config.PT_SQ_CUTOFF);

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

            assert(maxIter <= MAX_ITER);    // If maxIter > MAX_ITER, the solver failed

            const double trialPt2 = 0.5 * (bracket.first + bracket.second);

            if (trialPt2 < m_config.PT_SQ_CUTOFF)
                return -1.0;
            
            const double alphaSCorr = alphaSCMW(m_config, trialPt2) / alphaS0customLambda(m_config, trialPt2, 5, LAMBDA2);
            const double accRatio = alphaSCorr * vExactOverVTilde(born, trialPt2);
            assert(accRatio <= 1.0);

            if (rand() < accRatio)
                return trialPt2;
            else 
                ptMax2 = trialPt2;
        }

        assert(false); // Max trials exceeded!
        return -1.0;
    }

    double EmissionGenerator::vExactOverVTilde(const BornPhSpPt& born, double pt2) const
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

    double EmissionGenerator::integrateVTilde(
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

    RealOverBornContributions EmissionGenerator::getRealOverBornContributions(
        const RealPhSpPt& real,
        const BornPhSpPt& born,
        const double amp2Born,
        const double muF2,
        const double muR2
    ) const
    {
        RealOverBornContributions contributions;

        const double lumBorn = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar, muF2) / born.x1Bar
            * m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar, muF2) / born.x2Bar;

        for (const auto& realChannel: m_process.realChannels(born.channel))
        {
            const double lumReal = m_config.PDF->xfxQ2(realChannel.id1, real.x1, muF2) / real.x1
                * m_config.PDF->xfxQ2(realChannel.id2, real.x2, muF2) / real.x2;

            const double realAmp2 = m_process.realAmp2(real, realChannel, alphaSCMW(m_config, muR2));
            const double realOverBornPartonic = real.radJacobian * born.sHat / real.sHatReal * realAmp2 / amp2Born;
            const double realOverBorn = lumReal / lumBorn * realOverBornPartonic;

            contributions.total += realOverBorn;
            contributions.channels.push_back({ realChannel, realOverBorn });
        }

        return contributions;
    }

} // namespace powheg_dy
