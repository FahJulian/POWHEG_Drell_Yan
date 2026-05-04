#include "matrix_elements.h"
#include "couplings.h"

#include "powheg_dy/math.h"
#include "powheg_dy/process.h"
#include "powheg_dy/bra_ket.h"

namespace powheg_dy
{
    namespace 
    {
        struct _NeutralCurrentCouplings
        {
            double symmetric = 0.0;        // multiplies (1 + cos^2 theta)
            double forwardBackward = 0.0;  // multiplies 2 cos theta
        };

        _NeutralCurrentCouplings _neutralCurrentCouplings(const Process& process, int flavour, double sHat)
        {
            bool upType = flavour % 2 == 0;

            // Quark charge and axial/vector Z couplings.
            double qQ = upType ? 2.0 / 3.0 : -1.0 / 3.0;
            double aQ = upType ? 0.5 : -0.5;
            double vQ = upType ? 0.5 - (4.0 / 3.0) * process.S_W_SQ()
                              : -0.5 + (2.0 / 3.0) * process.S_W_SQ();

            double dz = (sHat - process.M_Z() * process.M_Z())
                      * (sHat - process.M_Z() * process.M_Z())
                      + process.M_Z() * process.M_Z() * process.GAMMA_Z() * process.GAMMA_Z();

            double reChi = process.KAPPA() * sHat * (sHat - process.M_Z() * process.M_Z()) / dz;
            double absChiSq = process.KAPPA() * process.KAPPA() * sHat * sHat / dz;

            _NeutralCurrentCouplings couplings;

            couplings.symmetric = qQ * qQ
                - 2.0 * qQ * process.V_L() * vQ * reChi
                + (process.V_L() * process.V_L() + process.A_L() * process.A_L())
                * (vQ * vQ + aQ * aQ) * absChiSq;

            couplings.forwardBackward = -2.0 * qQ * process.A_L() * aQ * reChi
                + 4.0 * process.V_L() * process.A_L() * vQ * aQ * absChiSq;

            return couplings;
        }

        double _bornAmp2(
            const Process& process,
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

            const auto zQ = _zCouplings(process, upType, qQ);
            const auto zL = _zCouplings(process, false, qL);

            const double sw = std::sqrt(process.S_W_SQ());
            const double cw = std::sqrt(process.C_W_SQ());

            const double q2 = (pLMinus + pLPlus).square();

            const std::complex<double> propGamma = _photonPropagator(q2);
            const std::complex<double> propZ     = _zPropagator(process, q2);

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

            const double e2 = 4.0 * PI * process.ALPHA();

            return amp2 * e2 * e2 / (4.0 * process.NC());
        }

    } // namespace 

    double MatrixElements::bornAngularFactorOld(const Process& process, int flavour, double sHat, double cosTheta)
    {
        auto couplings = _neutralCurrentCouplings(process, flavour, sHat);
        return couplings.symmetric * (1.0 + cosTheta * cosTheta)
                + 2.0 * couplings.forwardBackward * cosTheta;
    }

    double MatrixElements::born(const Process& process, const BornPhSpPt& born, double muR2)
    {
        // arrange momenta in the same convention as POWHEG's Fortran:
        // argument 1 is the quark, argument 2 is the antiquark.
        const bool leg1IsQuark = (born.channel.id1 > 0);

        const FourVector pQ    = leg1IsQuark ? born.p1Bar  : born.p2Bar;
        const FourVector pQbar = leg1IsQuark ? born.p2Bar  : born.p1Bar;

        return _bornAmp2(
            process,
            born.channel.flavour,
            pQ,
            pQbar,
            born.pLMinus,
            born.pLPlus
        );
    }

} // namespace powheg_dy
