#include "config.h"

#include "powheg_dy/base.h"

namespace powheg_dy
{
namespace 
{
    double lambda5FromAlphaSNNLO(const Config& config, double alphaSMZ)
    {
        const int nF = 5;

        const double b0 = config.beta0(nF);
        const double b1 = (153.0 - 19.0 * nF) / (24.0 * PI * PI);
        const double b2 = (2857.0 / 2.0 - 5033.0 / 18.0 * nF + 325.0 / 54.0 * nF * nF)
            / std::pow(4.0 * PI, 3);

        double t = 1.0 / (b0 * alphaSMZ);

        for (int iter = 0; iter < 10000; ++iter)
        {
            const double old = t;
            const double lt = std::log(t);

            const double c = b1 / (b0 * b0);

            const double as = 1.0 / (t * b0)
                * (
                    1.0
                    - c * lt / t
                    + std::pow(c * lt / t, 2)
                    - (b1 * b1 * (lt + 1.0) - b0 * b2)
                        / std::pow(b0, 4)
                        / (t * t)
                );

            // Numerical derivative is good enough here and avoids copying POWHEG's long as1.
            const double eps = 1e-6 * std::max(1.0, std::abs(t));
            auto alphaAt = [&](double tt)
            {
                const double ltt = std::log(tt);
                const double cc = b1 / (b0 * b0);
                return
                    1.0 / (tt * b0)
                    * (
                        1.0
                        - cc * ltt / tt
                        + std::pow(cc * ltt / tt, 2)
                        - (b1 * b1 * (ltt + 1.0) - b0 * b2)
                            / std::pow(b0, 4)
                            / (tt * tt)
                    );
            };

            const double das_dt = (alphaAt(t + eps) - alphaAt(t - eps)) / (2.0 * eps);
            t += (alphaSMZ - as) / das_dt;

            if (std::abs((old - t) / old) <= 1e-8)
                break;
        }

        return config.M_Z / std::exp(0.5 * t);
    }

    double lambda5FromAlphaSNLO(const Config& config, double alphaSMZ)
    {
        const int nF = 5;

        const double b0 = config.beta0(nF);
        const double bp = (153.0 - 19.0 * nF) / (2.0 * PI * (33.0 - 2.0 * nF));

        // POWHEG initial guess
        double t = 1.0 / (b0 * alphaSMZ);

        while (true)
        {
            const double old = t;
            const double logt = std::log(t);

            const double as0 =
                1.0 / (b0 * t)
                - bp * logt / ((b0 * t) * (b0 * t));

            const double das_dt =
                -1.0 / (b0 * t * t)
                - (bp / (b0 * b0)) * (1.0 - 2.0 * logt) / (t * t * t);

            t += (alphaSMZ - as0) / das_dt;

            if (std::abs((old - t) / old) <= 1e-8)
                break;
        }

        return config.M_Z / std::exp(0.5 * t);
    }

} // anonymous namespace

    void setDependentParams(Config& config)
    {
        config.S_W_SQ = 1.0 - config.M_W * config.M_W / config.M_Z / config.M_Z;
        config.C_W_SQ = 1.0 - config.S_W_SQ;
        config.S_W = std::sqrt(config.S_W_SQ);
        config.C_W = std::sqrt(config.C_W_SQ);

        config.E_SQ = 4.0 * PI * config.ALPHA_EW;

        config.S = config.SQRT_S * config.SQRT_S;
    }

    void extractLambdaFromPdf(Config& config)
    {
        const double alphaSMZ = config.PDF->alphasQ(config.M_Z);

        // LHAPDF6 exposes the alpha_s order through info metadata.
        // Depending on the set/version this can be "AlphaS_Order".
        const int iord = config.PDF->info().get_entry_as<int>("AlphaS_OrderQCD");

        double lambdaMsb5;
        switch (iord)
        {
        case 0:  lambdaMsb5 = lambda5FromAlphaSNLO(config, alphaSMZ); break;
        case 1:  lambdaMsb5 = lambda5FromAlphaSNLO(config, alphaSMZ); break;
        case 2:  lambdaMsb5 = lambda5FromAlphaSNNLO(config, alphaSMZ); break;
        default: lambdaMsb5 = lambda5FromAlphaSNLO(config, alphaSMZ); break;
        }

        Log::info << "LAMBDA_QCD extracted from LHAPDF: " << lambdaMsb5 << std::endl; 
        
        // std::cout << "LAMBDA5 is " << lambdaMsb5 << std::endl;

        config.LAMBDA_MSB_5_SQ = lambdaMsb5 * lambdaMsb5;
    }

} // namespace powheg_dy
