#include "emission_generator.h"

#include "powheg_dy/math/math.h"
#include "powheg_dy/math/rand.h"
#include "powheg_dy/process.h"

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
    constexpr double N_Q = 7.5;
    constexpr int __MAX_TRIALS_V = 5000;

    Emission __makeAcceptedEmission(
        const RealPhSpPt& real,
        RadiationChannel channel,
        double exact,
        double upper,
        double ratio
    )
    {
        Emission emission;

        emission.rad = real.rad;
        emission.channel = channel;
        emission.kt2 = real.kt2;

        emission.exactDensity = exact;
        emission.upperDensity = upper;
        emission.acceptanceRatio = ratio;

        emission.rejected = false;

        return emission;
    }

    RadiationChannel __chooseChannel(const MatrixElements::RealOverBornContributions& contributions,
        const BornChannel& bornChannel)
    {
        // Choose radiation channel, i.e. the emitted parton based on relative weight
        const double u = rand(0.0, contributions.total());

        if (u < contributions.qqbar)
            return { bornChannel.id1, bornChannel.id2, 21 };
        else if (u < contributions.qqbar + contributions.gqbar)
            return { 21, bornChannel.id2, bornChannel.id2 };
        else
            return { bornChannel.id1, 21, bornChannel.id1 };
    }

} // namespace

    Emission EmissionGenerator::generateEmission(
        const BornPhSpPt& born,
        int region
    ) const
    {
        if (region != 1)
            assert(false);      // Only ISR Radiation is implemented

        return _generateEmissionISR(born);
    }

    Emission EmissionGenerator::_generateEmissionISR(const BornPhSpPt& born
    ) const
    {
        double kt2Max = _globalKt2Max(born);

        double logR = 0;
        for (int trial = 0; trial < __MAX_TRIALS && kt2Max > m_process.pt2Cutoff(); ++trial)
        {
            const double kt2Trial = _sampleTrialKt2(born, kt2Max, logR);
            if (kt2Trial < m_process.pt2Cutoff())
                return Emission().reject();
            
            const RadiationVariables rad = _sampleTrialRadiation(born, kt2Trial);

            const double xiMax = m_realPhaseSpace.xiMax(born, rad.y);
            if (rad.xi <= __XI_EPS || rad.xi >= 1.0 - __XI_EPS || rad.xi >= xiMax)
            {
                // kt2Max = kt2Trial;
                continue;
            }

            const RealPhSpPt real = m_realPhaseSpace.reconstruct(born, rad);

            const auto contributions = MatrixElements::realOverBornContributions(m_process, real, kt2Trial, kt2Trial, true);

            const double exact = real.radJacobian * contributions.total();
            const double upper = _upperRadiationDensity(real, kt2Trial);

            assert(upper > 0);
            const double ratio = exact / upper;

            if (ratio > 1.0)
            {
                std::cout << "WARNING: Acceptance ratio is greater than one, accepting emission." << std::endl;
                
                auto channel = __chooseChannel(contributions, born.channel);

                return __makeAcceptedEmission(
                    real,
                    channel,
                    exact,
                    upper,
                    ratio
                );
            }

            if (rand() < ratio)
            {
                auto channel = __chooseChannel(contributions, born.channel);

                return __makeAcceptedEmission(
                    real,
                    channel,
                    exact,
                    upper,
                    ratio
                );
            }
            // else
            //     kt2Max = kt2Trial;
        }
        
        return Emission().reject();
    }

    double EmissionGenerator::_globalKt2Max(const BornPhSpPt& born) const
    {
        const double sumXbar = born.x1Bar + born.x2Bar;
        return born.sHat * (1 - born.x2Bar * born.x2Bar) * (1 - born.x1Bar * born.x1Bar)
            / (sumXbar * sumXbar);
    }

    RadiationVariables EmissionGenerator::_sampleTrialRadiation(
        const BornPhSpPt& born,
        double kt2Trial
    ) const
    {
        RadiationVariables rad;

        rad.xi = _sampleTrialXi(born, kt2Trial);
        rad.y = _sampleY(born, kt2Trial, rad.xi);
        rad.phi = _sampleTrialPhi();

        return rad;
    }

    double EmissionGenerator::_upperRadiationDensity(const RealPhSpPt& real, double kt2Trial) const
    {
        const double alphaS = m_process.alphaSCMW(kt2Trial);

        return N_Q
            * alphaS
            / real.rad.xi 
            / (1 - real.rad.y * real.rad.y);
    }

    double EmissionGenerator::_sampleTrialPhi() const
    {
        return rand() * 2.0 * PI;
    }

    // double EmissionGenerator::_sampleTrialXi(const BornPhSpPt& born, double pT2) const
    // {
    //     const double b = sqrt(pT2) / born.mB;
    //     const double a = sqrt(1.0 + b * b);
    //     const double xiPlus = 1.0 - (a - b) * (a - b);
    //     const double xiMinus = 1.0 - (a + b) * (a + b);
    //     const double rho = born.sHat / m_process.S();

    //     const double delta = xiPlus - xiMinus;
    //     const double eta0 = sqrt(delta);
    //     const double etaMax = sqrt((1 - rho) - xiPlus) + sqrt((1 - rho) - xiMinus);

    //     const double u = rand();
    //     const double eta = eta0 * pow(etaMax / eta0, u);

    //     const double bracket = (eta - delta / eta);
    //     return xiPlus + bracket * bracket / 4.0;
    // }

    double EmissionGenerator::_sampleTrialXi(const BornPhSpPt& born, double pT2) const
    {
        const double r = pT2 / born.sHat;
        const double a = std::sqrt(1.0 + r);
        const double b = std::sqrt(r);

        const double xPlus  = (a + b) * (a + b);
        const double xMinus = (a - b) * (a - b);

        // This must match _VExactOverVTilde
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

    double EmissionGenerator::_sampleY(const BornPhSpPt& born, double pT2, double xi) const
    {
        const double xiFactor = (1.0 - xi) / xi / xi;
        const double absY = sqrt(1.0 - 4.0 * pT2 / born.sHat * xiFactor);
        const double sign = rand() > 0.5 ? 1.0 : -1.0;

        return sign * absY;
    }

    // sample a k_T^2 after the procedure describes in Appendix D of ref. [22]
    double EmissionGenerator::_sampleTrialKt2(const BornPhSpPt& born, double ptMax2, double& logR) const
    {
        assert(ptMax2 >= m_process.pt2Cutoff());

        const int nF = 5;
        const double LAMBDA2 = __LAMBDA2_OVEREST_FACTOR * m_process.LAMBDA_SQ_QCD();
        const double BETA0 = (11.0 * m_process.C_A() - 4.0 * m_process.T_F() * nF) / 12.0 / PI;

        // double logR = 0.0;
        for (int trial = 1; trial <= __MAX_TRIALS_V; trial++)
        {
            logR += log(rand());

            const double maxIntegral = _integrateVTilde(m_process.pt2Cutoff(), ptMax2, born.sHat, LAMBDA2, N_Q, BETA0);

            if (-logR > maxIntegral)
                return -1.0;

            auto f = [&](double pt2) -> double
            {
                return _integrateVTilde(pt2, ptMax2, born.sHat, LAMBDA2, N_Q, BETA0) + logR;
            };

            const double fLow  = f(m_process.pt2Cutoff());   // >= 0
            const double fHigh = logR;  // = -target < 0

            const boost::uintmax_t MAX_ITER = 100;
            boost::uintmax_t maxIter = MAX_ITER;
            boost::math::tools::eps_tolerance<double> tol(40);

            auto bracket = boost::math::tools::toms748_solve(
                f,
                m_process.pt2Cutoff(),
                ptMax2,
                fLow,
                fHigh,
                tol,
                maxIter
            );

            assert(maxIter <= MAX_ITER);    // If maxIter > MAX_ITER, the solver failed

            const double trialPt2 = 0.5 * (bracket.first + bracket.second);

            if (trialPt2 < m_process.pt2Cutoff())
                return -1.0;
            
            const double alphaS0 = 1.0 / BETA0 / std::log(trialPt2 / LAMBDA2);
            const double alphaSCorr = m_process.alphaSCMW(trialPt2) / alphaS0;
            const double accRatio = alphaSCorr * _VExactOverVTilde(born, trialPt2);
            assert(accRatio <= 1.0);

            if (rand() < accRatio)
                return trialPt2;
            // else 
                // ptMax2 = trialPt2;
        }

        assert(false); // Max trials exceeded!
        return -1.0;
    }

    double EmissionGenerator::_VExactOverVTilde(const BornPhSpPt& born, double pt2) const
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

    double EmissionGenerator::_integrateVTilde(
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
