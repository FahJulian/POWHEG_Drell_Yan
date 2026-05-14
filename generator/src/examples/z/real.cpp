#include "drell_yan.h"

#include "powheg_dy/config/config.h"
#include "powheg_dy/math/bra_ket.h"

namespace powheg_dy
{
    double DrellYanProcess::realAmp2_crossed(
        const int flavour,
        const FourVector& qSlot,
        const FourVector& qbarSlot,
        const FourVector& gluonSlot,
        const FourVector& pLMinus,
        const FourVector& pLPlus,
        const double aS,
        const double colorFactor
    ) const
    {
        const FourVector qSpinMom     = (qSlot.e     < 0.0) ? -qSlot     : qSlot;
        const FourVector qbarSpinMom  = (qbarSlot.e  < 0.0) ? -qbarSlot  : qbarSlot;
        const FourVector gluonPolMom  = (gluonSlot.e < 0.0) ? -gluonSlot : gluonSlot;

        const WeylSpinors psiQ    = ket(qSpinMom,    FermionType::FERMION);
        const WeylSpinors psiQbar = bra(qbarSpinMom, FermionType::ANTIFERMION);

        const WeylSpinors psiLm   = bra(pLMinus, FermionType::FERMION);
        const WeylSpinors psiLp   = ket(pLPlus,  FermionType::ANTIFERMION);
        
        const bool upType = flavour % 2 == 0;

        const double qQ = (upType) ? 2.0 / 3.0 : -1.0 / 3.0;
        const double qL = -1.0;

        const auto zQ = m_config.zCouplings(upType, qQ);
        const auto zL = m_config.zCouplings(false, qL);

        const double sw = std::sqrt(m_config.S_W_SQ);
        const double cw = std::sqrt(m_config.C_W_SQ);

        const double q2 = (pLMinus + pLPlus).square();

        const std::complex<double> propGamma = m_config.photonPropagator(q2);
        const std::complex<double> propZ     = m_config.zPropagator(q2);

        double amp2 = 0.0;
        for (int pol = 1; pol <= 2; pol++)
        {
            const FourVector epsG = polVec(gluonPolMom, 0.0, pol);

            for (int helL = -1; helL <= 1; helL += 2)
            {
                const FourVectorC jL = braGammaKet(psiLm, psiLp, helL);

                for (int helQ = -1; helQ <= 1; helQ += 2)
                {
                    const FourVectorC jQ = braGammaKetGluon(psiQbar, psiQ, helQ, -qbarSlot, qSlot, gluonSlot, epsG);
                    const std::complex<double> contraction = dot(jL, jQ);

                    const std::complex<double> ewFactor = qQ * qL * propGamma
                        + propZ * zQ[helQ] * zL[helL] / std::pow(2.0 * sw * cw, 2);

                    amp2 += std::norm(contraction * ewFactor);
                }   
            }
        }

        const double gs2 = 4.0 * PI * aS;

        return amp2 * m_config.E_SQ * m_config.E_SQ * gs2 * colorFactor / 4.0 / m_config.N_C;
    }

    double DrellYanProcess::realAmp2qqbar(const RealPhSpPt& real, const BornChannel& bornChannel, const double alphaS) const
    {
        FourVector qSlot;
        FourVector qbarSlot;
        FourVector gluonSlot;

        if (bornChannel.id1 > 0)
        {
            qSlot     = real.p1In;
            qbarSlot  = real.p2In;
            gluonSlot = real.pOut[2];
        }
        else
        {
            qSlot     = real.p2In;
            qbarSlot  = real.p1In;
            gluonSlot = real.pOut[2];
        }

        return realAmp2_crossed(
            std::abs(bornChannel.id1),
            qSlot,
            qbarSlot,
            gluonSlot,
            real.pOut[0],
            real.pOut[1],
            alphaS,
            m_config.C_F
        );
    }

    double DrellYanProcess::realAmp2gluonLeg1(const RealPhSpPt& real, const BornChannel& bornChannel, const double alphaS) const
    {
        FourVector qSlot;
        FourVector qbarSlot;
        FourVector gluonSlot;

        if (bornChannel.id1 > 0)
        {
            qSlot     = -real.pOut[2];
            qbarSlot  =  real.p2In;
            gluonSlot = -real.p1In;
        }
        else
        {
            qSlot     =  real.p2In;
            qbarSlot  = -real.pOut[2];
            gluonSlot = -real.p1In;
        }

        return realAmp2_crossed(
            std::abs(bornChannel.id1),
            qSlot,
            qbarSlot,
            gluonSlot,
            real.pOut[0],
            real.pOut[1],
            alphaS,
            m_config.T_F
        );
    }

    double DrellYanProcess::realAmp2gluonLeg2(const RealPhSpPt& real, const BornChannel& bornChannel, const double alphaS) const
    {
        FourVector qSlot;
        FourVector qbarSlot;
        FourVector gluonSlot;

        if (bornChannel.id1 > 0)
        {
            qSlot     =  real.p1In;
            qbarSlot  = -real.pOut[2];
            gluonSlot = -real.p2In;
        }
        else
        {
            qSlot     = -real.pOut[2];
            qbarSlot  =  real.p1In;
            gluonSlot = -real.p2In;
        }

        return realAmp2_crossed(
            std::abs(bornChannel.id1),
            qSlot,
            qbarSlot,
            gluonSlot,
            real.pOut[0],
            real.pOut[1],
            alphaS,
            m_config.T_F
        );
    }
    
    double DrellYanProcess::realAmp2(
        const RealPhSpPt& real, 
        const BornChannel& bornChannel,
        const RealChannel& realChannel, 
        const double alphaS) const 
    {
        if (std::abs(realChannel.id1) <= 5 && realChannel.id1 == -realChannel.id2 && realChannel.outIDs[2] == 21)
            return realAmp2qqbar(real, bornChannel, alphaS);
        else if (realChannel.id1 == 21 && std::abs(realChannel.id2) <= 5)
            return realAmp2gluonLeg1(real, bornChannel, alphaS);
        else if (realChannel.id2 == 21 && std::abs(realChannel.id1) <= 5)
            return realAmp2gluonLeg2(real, bornChannel, alphaS);
        else
            throw std::runtime_error("Invalid real channel");

        return 0.0;
    }

} // namespace powheg_dy
