#include "drell_yan.h"

#include "powheg/math/bra_ket.h"

namespace powheg
{
    double DrellYanProcess::bornAmp2(
        const int flavour,
        const FourVector& pQ,
        const FourVector& pQbar,
        const FourVector& pLMinus,
        const FourVector& pLPlus
    ) const
    {
        const WeylSpinors psiQ    = ket(pQ,      FermionType::FERMION);
        const WeylSpinors psiQbar = bra(pQbar,   FermionType::ANTIFERMION);
        const WeylSpinors psiLm   = bra(pLMinus, FermionType::FERMION);
        const WeylSpinors psiLp   = ket(pLPlus,  FermionType::ANTIFERMION);

        const bool upType = flavour % 2 == 0;

        const double qQ = (upType) ? 2.0 / 3.0 : -1.0 / 3.0;
        const double qL = -1.0;

        const auto zQ = m_config.zCouplings(upType, qQ);
        const auto zL = m_config.zCouplings(false, qL);

        const double q2 = (pLMinus + pLPlus).square();

        const std::complex<double> propGamma = m_config.photonPropagator(q2);
        const std::complex<double> propZ     = m_config.zPropagator(q2);

        double amp2 = 0.0;
        for (int helL = -1; helL <= +1; helL += 2)
        {
            const FourVectorC jL = braGammaKet(psiLm, psiLp, helL);

            for (int helQ = -1; helQ <= +1; helQ += 2)
            {
                const FourVectorC jQ = braGammaKet(psiQbar, psiQ, helQ);
                std::complex<double> contraction = dot(jL, jQ);

                const std::complex<double> ewFactor = qQ * qL * propGamma
                    + propZ * zQ[helQ] * zL[helL] / std::pow(2.0 * m_config.S_W * m_config.C_W, 2);

                amp2 += std::norm(contraction * ewFactor);
            }
        }

        return amp2 * m_config.E_SQ * m_config.E_SQ / (4.0 * m_config.N_C);
    }

    double DrellYanProcess::bornAmp2(const BornPhSpPt& born, const BornChannel& bornChannel) const
    {
        const bool leg1IsQuark = (bornChannel.id1 > 0);

        const FourVector pQ    = leg1IsQuark ? born.p1Bar : born.p2Bar;
        const FourVector pQbar = leg1IsQuark ? born.p2Bar : born.p1Bar;

        return bornAmp2(
            std::abs(bornChannel.id1),
            pQ,
            pQbar,
            born.pOut[0],
            born.pOut[1]
        );
    }

} // namespace powheg
