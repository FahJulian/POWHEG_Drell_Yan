#include "matrix_elements.h"

#include "couplings.h"

#include "powheg_dy/config.h"
#include "powheg_dy/alpha_s.h"
#include "powheg_dy/math/bra_ket.h"

#include <iostream> 

namespace powheg_dy
{
    namespace 
    {
        constexpr int __NF = 5;
        constexpr double __EPS = 1.0e-14;

        double _realAmp2_crossed(
            const Config& config,
            int flavour,
            const FourVector& qSlot,
            const FourVector& qbarSlot,
            const FourVector& gluonSlot,
            const FourVector& pLMinus,
            const FourVector& pLPlus,
            double alphaS,
            double colorFactor
        )
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

            const auto zQ = _zCouplings(config, upType, qQ);
            const auto zL = _zCouplings(config, false, qL);

            const double sw = std::sqrt(config.S_W_SQ);
            const double cw = std::sqrt(config.C_W_SQ);

            const double q2 = (pLMinus + pLPlus).square();

            const std::complex<double> propGamma = _photonPropagator(q2);
            const std::complex<double> propZ     = _zPropagator(config, q2);

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

            const double e2  = 4.0 * PI * config.ALPHA_EW;
            const double gs2 = 4.0 * PI * alphaS;

            return amp2 * e2 * e2 * gs2 * colorFactor / 4.0 / config.N_C;
        }

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

            const double sw = std::sqrt(config.S_W_SQ);
            const double cw = std::sqrt(config.C_W_SQ);

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
                        + propZ * zQ[helQ] * zL[helL] / std::pow(2.0 * sw * cw, 2);

                    amp2 += std::norm(contraction * ewFactor);
                }
            }

            const double e2 = 4.0 * PI * config.ALPHA_EW;

            return amp2 * e2 * e2 / (4.0 * config.N_C);
        }

    } // namespace

    MatrixElements::RealOverBornContributions MatrixElements::realOverBornContributions(
        const Config& config,
        const RealPhSpPt& real,
        double muF2,
        double muR2, 
        bool useCMWALphaS
    )
    {
        RealOverBornContributions out;

        const BornPhSpPt& born = real.underlyingBorn;
        const int id1 = born.channel.id1;
        const int id2 = born.channel.id2;

        const double fBorn1 = config.PDF->xfxQ2(id1, born.x1Bar, muF2) / born.x1Bar;
        const double fBorn2 = config.PDF->xfxQ2(id2, born.x2Bar, muF2) / born.x2Bar;
        const double bornLuminosity = fBorn1 * fBorn2;

        const double fReal1Q  = config.PDF->xfxQ2(id1, real.x1, muF2) / real.x1; 
        const double fReal2QB = config.PDF->xfxQ2(id2, real.x2, muF2) / real.x2; 

        const double fReal1G  = config.PDF->xfxQ2(21, real.x1, muF2) / real.x1; 
        const double fReal2G  = config.PDF->xfxQ2(21, real.x2, muF2) / real.x2; 

        const double lumQQbar = fReal1Q * fReal2QB;
        const double lumGQbar = fReal1G * fReal2QB;
        const double lumQG    = fReal1Q * fReal2G;

        out.qqbar = lumQQbar / bornLuminosity * realOverBornQQbar(config, real, muR2, useCMWALphaS);
        out.gqbar = lumGQbar / bornLuminosity * realOverBornGQbar(config, real, muR2, useCMWALphaS);
        out.qg = lumQG / bornLuminosity * realOverBornQG(config, real, muR2, useCMWALphaS);

        return out;
    }

    double MatrixElements::realOverBornQQbar(const Config& config, const RealPhSpPt& real, double muR2, bool useCMWALphaS)
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
        const double alphaS = useCMWALphaS ? alphaSCMW(config, muR2) : alphaSPowheg(config, muR2);

        const double realAmp2 = _realAmp2_crossed(
            config,
            born.channel.flavour,
            pQReal,
            pQbarReal,
            real.pRadiated,
            real.pLMinus,
            real.pLPlus,
            alphaS,
            config.C_F
        );

        const double bornAmp2 = _bornAmp2(
            config,
            born.channel.flavour,
            pQBorn,
            pQbarBorn,
            born.pLMinus,
            born.pLPlus
        );

        return born.sHat / real.sHatReal * realAmp2 / bornAmp2;
    }

    double MatrixElements::realOverBornGQbar(const Config& config, const RealPhSpPt& real, double muR2, bool useCMWALphaS)
    {
        const BornPhSpPt& born = real.underlyingBorn;

        const double alphaS = useCMWALphaS ? alphaSCMW(config, muR2) : alphaSPowheg(config, muR2);

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

        const double realAmp2 = _realAmp2_crossed(
            config,
            born.channel.flavour,
            qSlot,
            qbarSlot,
            gluonSlot,
            real.pLMinus,
            real.pLPlus,
            alphaS,
            config.T_F
        );

        const double bornAmp2 = _bornAmp2(
            config,
            born.channel.flavour,
            pQBorn,
            pQbarBorn,
            born.pLMinus,
            born.pLPlus
        );

        return born.sHat / real.sHatReal * realAmp2 / bornAmp2;
    }

    double MatrixElements::realOverBornQG(const Config& config, const RealPhSpPt& real, double muR2, bool useCMWALphaS)
    {
        const BornPhSpPt& born = real.underlyingBorn;

        const double alphaS = useCMWALphaS ? alphaSCMW(config, muR2) : alphaSPowheg(config, muR2);

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

        const double realAmp2 = _realAmp2_crossed(
            config,
            born.channel.flavour,
            qSlot,
            qbarSlot,
            gluonSlot,
            real.pLMinus,
            real.pLPlus,
            alphaS,
            config.T_F
        );

        const double bornAmp2 = _bornAmp2(
            config,
            born.channel.flavour,
            pQBorn,
            pQbarBorn,
            born.pLMinus,
            born.pLPlus
        );

        return born.sHat / real.sHatReal * realAmp2 / bornAmp2;
    }

} // namespace powheg_dy
