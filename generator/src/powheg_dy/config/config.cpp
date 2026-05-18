#include "config.h"

#include "powheg_dy/base.h"
#include "powheg_dy/config/alpha_s.h"

namespace powheg
{
namespace 
{
    double M_Z_FOR_LAMBDA_QCD_EXTRACTION = 91.1876;

    double lambda5FromAlphaSNNLO(int nF, double cA, double tF, double mZ, double alphaSMZ)
    {
        const double b0 = (11.0 * cA - 4.0 * tF * nF) / (12.0 * PI);
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

            // Numerical derivative to avoid copying POWHEG's long as1.
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

        return mZ / std::exp(0.5 * t);
    }

    double lambda5FromAlphaSNLO(int nF, double cA, double tF, double mZ, double alphaSMZ)
    {
        const double b0 = (11.0 * cA - 4.0 * tF * nF) / (12.0 * PI);
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

        return mZ / std::exp(0.5 * t);
    }

} // anonymous namespace

    void Config::initPdf()
    {
        std::stringstream buffer;

        auto* oldCout = std::cout.rdbuf(buffer.rdbuf());
        auto* oldCerr = std::cerr.rdbuf(buffer.rdbuf());
        PDF = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(PDF_NAME, 0));

        std::cout.rdbuf(oldCout);
        std::cerr.rdbuf(oldCerr);

        std::string line;
        while (std::getline(buffer, line, '\n'))
            Log::info << line << Log::endl;
    }

    void Config::extractLambdaFromPdf()
    {
        const double alphaSMZ = PDF->alphasQ(M_Z_FOR_LAMBDA_QCD_EXTRACTION);

        // LHAPDF6 exposes the alpha_s order through info metadata.
        // Depending on the set/version this can be "AlphaS_Order".
        const int iord = PDF->info().get_entry_as<int>("AlphaS_OrderQCD");

        double lambdaMsb5;
        switch (iord)
        {
        case 0:  lambdaMsb5 = lambda5FromAlphaSNLO (5, C_A, T_F, M_Z_FOR_LAMBDA_QCD_EXTRACTION, alphaSMZ); break;
        case 1:  lambdaMsb5 = lambda5FromAlphaSNLO (5, C_A, T_F, M_Z_FOR_LAMBDA_QCD_EXTRACTION, alphaSMZ); break;
        case 2:  lambdaMsb5 = lambda5FromAlphaSNNLO(5, C_A, T_F, M_Z_FOR_LAMBDA_QCD_EXTRACTION, alphaSMZ); break;
        default: lambdaMsb5 = lambda5FromAlphaSNLO (5, C_A, T_F, M_Z_FOR_LAMBDA_QCD_EXTRACTION, alphaSMZ); break;
        }

        Log::info << "LAMBDA_QCD extracted from LHAPDF: " << lambdaMsb5 << Log::endl; 
        
        LAMBDA_MSB_5_SQ = lambdaMsb5 * lambdaMsb5;
    }

    double Config::beta0(int nF) const
    {
        return (11.0 * C_A - 4.0 * T_F * nF) / (12.0 * PI);
    }

    double Config::alphaS0(double qSq, int nF) const
    {
        return powhegAlphaS0(*this, qSq, nF);
    }

    double Config::alphaS0customLambda(double qSq, int nF, double lambdaSq) const
    {
        return powhegAlphaS0customLambda(*this, qSq, nF, lambdaSq);
    }

    double Config::alphaS(double qSq) const
    {
        return powhegAlphaS(*this, qSq);
    }

    double Config::alphaSCMW(double qSq) const
    {
        return powhegAlphaSCMW(*this, qSq);
    }

} // namespace powheg_dy
