#include "matrix_elements.h"

#include "couplings.h"

#include "powheg_dy/math/bra_ket.h"

namespace powheg_dy
{
    namespace 
    {
        double _bornAmp2(
            const Config& config,
            int flavour,
            const FourVector& pQ,
            const FourVector& pQbar,
            const FourVector& pLMinus,
            const FourVector& pLPlus
        )
        {
            const WeylSpinors psiQ    = ket(pQ,      FermionType::FERMION);
            const WeylSpinors psiQbar = bra(pQbar,   FermionType::ANTIFERMION);
            const WeylSpinors psiLm   = bra(pLMinus, FermionType::FERMION);
            const WeylSpinors psiLp   = ket(pLPlus,  FermionType::ANTIFERMION);

            const bool upType = flavour % 2 == 0;

            const double qQ = (upType) ? 2.0 / 3.0 : -1.0 / 3.0;
            const double qL = -1.0;

            const auto zQ = _zCouplings(config, upType, qQ);
            const auto zL = _zCouplings(config, false, qL);

            const double q2 = (pLMinus + pLPlus).square();

            const std::complex<double> propGamma = _photonPropagator(q2);
            const std::complex<double> propZ     = _zPropagator(config, q2);

            double amp2 = 0.0;
            for (int helL = -1; helL <= +1; helL += 2)
            {
                const FourVectorC jL = braGammaKet(psiLm, psiLp, helL);

                for (int helQ = -1; helQ <= +1; helQ += 2)
                {
                    const FourVectorC jQ = braGammaKet(psiQbar, psiQ, helQ);
                    std::complex<double> contraction = dot(jL, jQ);

                    const std::complex<double> ewFactor = qQ * qL * propGamma
                        + propZ * zQ[helQ] * zL[helL] / std::pow(2.0 * config.S_W * config.C_W, 2);

                    amp2 += std::norm(contraction * ewFactor);
                }
            }

            return amp2 * config.E_SQ * config.E_SQ / (4.0 * config.N_C);
        }

    } // namespace 

    double MatrixElements::born(const Config& config, const BornPhSpPt& born)
    {
        // arrange momenta in the same convention as POWHEG's Fortran:
        // argument 1 is the quark, argument 2 is the antiquark.
        const bool leg1IsQuark = (born.channel.id1 > 0);

        const FourVector pQ    = leg1IsQuark ? born.p1Bar  : born.p2Bar;
        const FourVector pQbar = leg1IsQuark ? born.p2Bar  : born.p1Bar;

        return _bornAmp2(
            config,
            born.channel.flavour,
            pQ,
            pQbar,
            born.pLMinus,
            born.pLPlus
        );
    }

} // namespace powheg_dy
