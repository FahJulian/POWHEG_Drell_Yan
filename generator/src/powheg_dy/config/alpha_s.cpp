#include "alpha_s.h"

#include "powheg_dy/config/config.h"

namespace powheg_dy
{
namespace 
{
    int activeFlavours(const Config& config, double qSq)
    {
        if (qSq < config.CMW_CHARM_TRSHLD_SQ) 
            return 3;
        else if (qSq < config.CMW_BOTTOM_TRSHLD_SQ) 
            return 4;
        else 
            return 5;
    }

    double betaPrimePowheg(int nf)
    {
        return (153.0 - 19.0 * nf) / (2.0 * PI * (33.0 - 2.0 * nf));
    }

    double alphaRawPowheg(const Config& config, double qSq, int nF)
    {
        const double L = std::log(qSq / config.LAMBDA_MSB_5_SQ);
        powheg_assert(L > 0.0);

        const double b  = config.beta0(nF);
        const double bp = betaPrimePowheg(nF);

        return 1.0 / (b * L) - bp * std::log(L) / ((b * L) * (b * L));
    }

} // anonymous namespace

    double powhegAlphaS0(const Config& config, double qSq, int nF)
    {
        return 1.0 / (config.beta0(nF) * std::log(qSq / config.LAMBDA_MSB_5_SQ));
    }

    double powhegAlphaS0customLambda(const Config& config, double qSq, int nF, double lambdaSq)
    {
        return 1.0 / (config.beta0(nF) * std::log(qSq / lambdaSq));
    }

    double powhegAlphaS(const Config& config, double qSq) 
    {
        if (config.ALPHA_S_FROM_PDF)
            return config.PDF->alphasQ2(qSq);

        const double a5_mb_raw = alphaRawPowheg(config, config.CMW_BOTTOM_TRSHLD_SQ, 5);
        const double a4_mb_raw = alphaRawPowheg(config, config.CMW_BOTTOM_TRSHLD_SQ, 4);

        const double c45 = 1.0 / a5_mb_raw - 1.0 / a4_mb_raw;

        const double a4_mc_raw = alphaRawPowheg(config, config.CMW_CHARM_TRSHLD_SQ, 4);
        const double a3_mc_raw = alphaRawPowheg(config, config.CMW_CHARM_TRSHLD_SQ, 3);

        const double c35 = 1.0 / a4_mc_raw - 1.0 / a3_mc_raw + c45;

        const int nF = activeFlavours(config, qSq);

        if (nF == 5)
        {
            return alphaRawPowheg(config, qSq, 5);
        }
        else if (nF == 4)
        {
            const double a4_raw = alphaRawPowheg(config, qSq, 4);
            return 1.0 / (1.0 / a4_raw + c45);
        }
        else
        {
            const double a3_raw = alphaRawPowheg(config, qSq, 3);
            return 1.0 / (1.0 / a3_raw + c35);
        }
    }

    double powhegAlphaSCMW(const Config& config, double qSq) 
    {
        const int nF = activeFlavours(config, qSq);
        const double alphaS = powhegAlphaS(config, qSq);

        const double K = (67.0 / 18.0 - PI * PI / 6.0) * config.C_A - 5.0 / 9.0 * nF;

        return alphaS * (1.0 + alphaS / (2.0 * PI) * K);
    }

} // namespace powheg_dy