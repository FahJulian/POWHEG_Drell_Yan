#include "drell_yan.h"
#include "couplings.h"

#include "powheg_dy/config.h"
#include "powheg_dy/math/bra_ket.h"

namespace powheg_dy
{
namespace 
{
    constexpr int NF = 5;
    constexpr double EPS = 1.0e-14;
    
} // anonymous namespace

    double DrellYanProcess::realAmp2_crossed(
        int flavour,
        const FourVector& qSlot,
        const FourVector& qbarSlot,
        const FourVector& gluonSlot,
        const FourVector& pLMinus,
        const FourVector& pLPlus,
        double aS,
        double colorFactor
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

        const auto zQ = zCouplings(m_config, upType, qQ);
        const auto zL = zCouplings(m_config, false, qL);

        const double sw = std::sqrt(m_config.S_W_SQ);
        const double cw = std::sqrt(m_config.C_W_SQ);

        const double q2 = (pLMinus + pLPlus).square();

        const std::complex<double> propGamma = photonPropagator(q2);
        const std::complex<double> propZ     = zPropagator(m_config, q2);

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

    RealOverBornContributionsOld DrellYanProcess::realOverBornContributions(
        const RealPhSpPt& real,
        double muF2,
        double muR2, 
        bool useCMWALphaS
    ) const
    {
        RealOverBornContributionsOld out;

        const BornPhSpPt& born = real.underlyingBorn;
        const int id1 = born.channel.id1;
        const int id2 = born.channel.id2;

        const double fBorn1 = m_config.PDF->xfxQ2(id1, born.x1Bar, muF2) / born.x1Bar;
        const double fBorn2 = m_config.PDF->xfxQ2(id2, born.x2Bar, muF2) / born.x2Bar;
        const double bornLuminosity = fBorn1 * fBorn2;

        const double fReal1Q  = m_config.PDF->xfxQ2(id1, real.x1, muF2) / real.x1; 
        const double fReal2QB = m_config.PDF->xfxQ2(id2, real.x2, muF2) / real.x2; 

        const double fReal1G  = m_config.PDF->xfxQ2(21, real.x1, muF2) / real.x1; 
        const double fReal2G  = m_config.PDF->xfxQ2(21, real.x2, muF2) / real.x2; 

        const double lumQQbar = fReal1Q * fReal2QB;
        const double lumGQbar = fReal1G * fReal2QB;
        const double lumQG    = fReal1Q * fReal2G;

        out.qqbar = lumQQbar / bornLuminosity * realOverBornQQbar(real, muR2, useCMWALphaS);
        out.gluonLeg1 = lumGQbar / bornLuminosity * realOverBornGQbar(real, muR2, useCMWALphaS);
        out.gluonLeg2 = lumQG / bornLuminosity * realOverBornQG(real, muR2, useCMWALphaS);

        return out;
    }

    double DrellYanProcess::realOverBornQQbar(const RealPhSpPt& real, double muR2, bool useCMWALphaS) const
    {
        const auto& born = real.underlyingBorn;

        // arrange momenta in the same convention as POWHEG's Fortran:
        // argument 1 is the quark, argument 2 is the antiquark.
        const bool leg1IsQuark = (born.channel.id1 > 0);

        const FourVector pQReal    = leg1IsQuark ? real.p1In  : real.p2In;
        const FourVector pQbarReal = leg1IsQuark ? real.p2In  : real.p1In;

        const FourVector pQBorn    = leg1IsQuark ? born.p1Bar : born.p2Bar;
        const FourVector pQbarBorn = leg1IsQuark ? born.p2Bar : born.p1Bar;

        // TODO: Change to the corrected alphaS
        const double aS = useCMWALphaS ? m_config.alphaSCMW(muR2) : m_config.alphaS(muR2);

        const double realamp = realAmp2_crossed(
            born.channel.flavour,
            pQReal,
            pQbarReal,
            real.pRadiated,
            real.pLMinus,
            real.pLPlus,
            aS,
            m_config.C_F
        );

        const double bornamp = bornAmp2(
            born.channel.flavour,
            pQBorn,
            pQbarBorn,
            born.pLMinus,
            born.pLPlus
        );

        return born.sHat / real.sHatReal * realamp / bornamp;
    }

    double DrellYanProcess::realOverBornGQbar(const RealPhSpPt& real, double muR2, bool useCMWALphaS) const
    {
        const BornPhSpPt& born = real.underlyingBorn;

        const double aS = useCMWALphaS ? m_config.alphaSCMW(muR2) : m_config.alphaS(muR2);

        FourVector qSlot;
        FourVector qbarSlot;
        FourVector gluonSlot;

        const bool leg1IsQuark = (born.channel.id1 > 0);

        const FourVector pQBorn    = leg1IsQuark ? born.p1Bar : born.p2Bar;
        const FourVector pQbarBorn = leg1IsQuark ? born.p2Bar : born.p1Bar;

        if (born.channel.id1 > 0)
        {
            qSlot     = -real.pRadiated;
            qbarSlot  =  real.p2In;
            gluonSlot = -real.p1In;
        }
        else
        {
            qSlot     =  real.p2In;
            qbarSlot  = -real.pRadiated;
            gluonSlot = -real.p1In;
        }

        const double realamp = realAmp2_crossed(
            born.channel.flavour,
            qSlot,
            qbarSlot,
            gluonSlot,
            real.pLMinus,
            real.pLPlus,
            aS,
            m_config.T_F
        );

        const double bornamp = bornAmp2(
            born.channel.flavour,
            pQBorn,
            pQbarBorn,
            born.pLMinus,
            born.pLPlus
        );

        return born.sHat / real.sHatReal * realamp / bornamp;
    }

    double DrellYanProcess::realOverBornQG(const RealPhSpPt& real, double muR2, bool useCMWALphaS) const
    {
        const BornPhSpPt& born = real.underlyingBorn;

        const double aS = useCMWALphaS ? m_config.alphaSCMW(muR2) : m_config.alphaS(muR2);

        FourVector qSlot;
        FourVector qbarSlot;
        FourVector gluonSlot;

        const bool leg1IsQuark = (born.channel.id1 > 0);
        
        const FourVector pQBorn    = leg1IsQuark ? born.p1Bar : born.p2Bar;
        const FourVector pQbarBorn = leg1IsQuark ? born.p2Bar : born.p1Bar;

        if (born.channel.id1 > 0)
        {
            qSlot     =  real.p1In;
            qbarSlot  = -real.pRadiated;
            gluonSlot = -real.p2In;
        }
        else
        {
            qSlot     = -real.pRadiated;
            qbarSlot  =  real.p1In;
            gluonSlot = -real.p2In;
        }

        const double realamp = realAmp2_crossed(
            born.channel.flavour,
            qSlot,
            qbarSlot,
            gluonSlot,
            real.pLMinus,
            real.pLPlus,
            aS,
            m_config.T_F
        );

        const double bornamp = bornAmp2(
            born.channel.flavour,
            pQBorn,
            pQbarBorn,
            born.pLMinus,
            born.pLPlus
        );

        return born.sHat / real.sHatReal * realamp / bornamp;
    }

    double DrellYanProcess::realAmp2qqbar(const RealPhSpPt& real, const double alphaS) const
    {
        FourVector qSlot;
        FourVector qbarSlot;
        FourVector gluonSlot;

        if (real.underlyingBorn.channel.id1 > 0)
        {
            qSlot     = real.p1In;
            qbarSlot  = real.p2In;
            gluonSlot = real.pRadiated;
        }
        else
        {
            qSlot     = real.p2In;
            qbarSlot  = real.p1In;
            gluonSlot = real.pRadiated;
        }

        return realAmp2_crossed(
            real.underlyingBorn.channel.flavour,
            qSlot,
            qbarSlot,
            gluonSlot,
            real.pLMinus,
            real.pLPlus,
            alphaS,
            m_config.C_F
        );
    }

    double DrellYanProcess::realAmp2gluonLeg1(const RealPhSpPt& real, const double alphaS) const
    {
        FourVector qSlot;
        FourVector qbarSlot;
        FourVector gluonSlot;

        if (real.underlyingBorn.channel.id1 > 0)
        {
            qSlot     = -real.pRadiated;
            qbarSlot  =  real.p2In;
            gluonSlot = -real.p1In;
        }
        else
        {
            qSlot     =  real.p2In;
            qbarSlot  = -real.pRadiated;
            gluonSlot = -real.p1In;
        }

        return realAmp2_crossed(
            real.underlyingBorn.channel.flavour,
            qSlot,
            qbarSlot,
            gluonSlot,
            real.pLMinus,
            real.pLPlus,
            alphaS,
            m_config.T_F
        );
    }

    double DrellYanProcess::realAmp2gluonLeg2(const RealPhSpPt& real, const double alphaS) const
    {
        FourVector qSlot;
        FourVector qbarSlot;
        FourVector gluonSlot;

        if (real.underlyingBorn.channel.id1 > 0)
        {
            qSlot     =  real.p1In;
            qbarSlot  = -real.pRadiated;
            gluonSlot = -real.p2In;
        }
        else
        {
            qSlot     = -real.pRadiated;
            qbarSlot  =  real.p1In;
            gluonSlot = -real.p2In;
        }

        return realAmp2_crossed(
            real.underlyingBorn.channel.flavour,
            qSlot,
            qbarSlot,
            gluonSlot,
            real.pLMinus,
            real.pLPlus,
            alphaS,
            m_config.T_F
        );
    }
    
    double DrellYanProcess::realAmp2(
        const RealPhSpPt& real, 
        const RealChannel& channel, 
        const double alphaS) const 
    {
        if (std::abs(channel.id1) <= 5 && channel.id1 == -channel.id2 && channel.outIDs[2] == 21)
            return realAmp2qqbar(real, alphaS);
        else if (channel.id1 == 21 && std::abs(channel.id2) <= 5)
            return realAmp2gluonLeg1(real, alphaS);
        else if (channel.id2 == 21 && std::abs(channel.id1) <= 5)
            return realAmp2gluonLeg2(real, alphaS);
        else
            throw std::runtime_error("Invalid real channel");

        return 0.0;
    }

} // namespace powheg_dy
