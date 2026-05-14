#include "bra_ket.h"

namespace powheg_dy
{
namespace 
{
    static constexpr double TINY = 1.0e-20;
    static constexpr std::complex<double> I = { 0.0, 1.0 };

    std::array<std::complex<double>, 2> slashKet(
        const FourVectorC& p, 
        const std::array<std::complex<double>, 2>& ket, 
        int helicity)
    {
        powheg_assert(helicity == 1 || helicity == -1);
        const double h = static_cast<double>(helicity);
        
        const std::complex<double> pSlash11 = p.x0 - h * p.x3;
        const std::complex<double> pSlash12 = -h * (p.x1 - I * p.x2);
        const std::complex<double> pSlash21 = -h * (p.x1 + I * p.x2);
        const std::complex<double> pSlash22 = p.x0 + h * p.x3;

        return {
            pSlash11 * ket[0] + pSlash12 * ket[1],
            pSlash21 * ket[0] + pSlash22 * ket[1]
        };
    }

    // std::array<std::complex<double>, 2> slashBra(
    //     const FourVectorC& p, 
    //     const std::array<std::complex<double>, 2>& bra, 
    //     int helicity)
    // {
    //     powheg_assert(helicity == 1 || helicity == -1);
    //     const double h = static_cast<double>(helicity);
        
    //     const std::complex<double> pSlash11 = p.x0 - h * p.x3;
    //     const std::complex<double> pSlash12 = -h * (p.x1 - I * p.x2);
    //     const std::complex<double> pSlash21 = -h * (p.x1 + I * p.x2);
    //     const std::complex<double> pSlash22 = p.x0 + h * p.x3;

    //     return {
    //         pSlash11 * bra[0] + pSlash21 * bra[1],
    //         pSlash12 * bra[0] + pSlash22 * bra[1]
    //     };
    // }

    std::complex<double> braSlashKet(
        const WeylSpinors& bra,
        const WeylSpinors& ket,
        int ketHel,
        const std::vector<FourVectorC>& momenta
    )
    {
        powheg_assert(ketHel == 1 || ketHel == -1);

        const int n = static_cast<int>(momenta.size());
        const int minusOneToN = (n % 2 == 0) ? 1 : -1;
        const int braHel = -minusOneToN * ketHel;

        const auto& b = (braHel == 1) ? bra.plus : bra.minus;
        auto k        = (ketHel == 1) ? ket.plus : ket.minus;

        for (int i = n - 1; i >= 0; --i)
        {
            k = slashKet(momenta[i], k, ketHel);
            ketHel = -ketHel;
        }

        return b[0] * k[0] + b[1] * k[1];
    }

} // anonymous namespace

    FourVector polVec(const FourVector& p, double mass, int pol)
    {
        FourVector epsilon = { };
        epsilon.e = 0.0;

        if (pol == 3)
        {
            powheg_assert(p.square() > 0.0);     // Long. pol. not allowed for massless gauge boson

            const double absP = sqrt(p.getThreeVec().square());
            powheg_assert(absP > 0.0);
            
            const double normal = mass / p.e / absP;

            epsilon.x = normal * p.x;
            epsilon.y = normal * p.y;
            epsilon.z = normal * p.z;
        }
        else if (pol == 2)
        {
            const double pT = std::sqrt(p.x * p.x + p.y * p.y);

            if (pT > 0.0)
            {
                epsilon.x = -p.y / pT;
                epsilon.y =  p.x / pT;
                epsilon.z = 0.0;
            }
            else 
            {
                epsilon.x = 0.0;
                epsilon.y = 1.0;
                epsilon.z = 0.0;
            }
        }
        else if (pol == 1)
        {
            const double pT = std::sqrt(p.x * p.x + p.y * p.y);

            if (pT > 0)
            {
                const double absP = (mass > 0) ? sqrt(p.getThreeVec().square()) : p.e;
                const double normal = 1.0 / absP / pT;
                
                epsilon.x = p.x * p.z * normal;
                epsilon.y = p.y * p.z * normal;
                epsilon.z = -pT *  pT * normal;
            }
            else 
            {
                epsilon.x = p.z > 0.0 ? 1.0 : -1.0;
                epsilon.y = 0.0;
                epsilon.z = 0.0;
            }
        }
        else 
        {
            powheg_assert(false);          // Invariant pol. mode
            epsilon = { };
        }

        return epsilon;
    }

    WeylSpinors ket(const FourVector& p, FermionType type)
    {
        WeylSpinors psi = { };

        const double ePlusPz = p.z > 0.0 ? p.e + p.z : (p.x * p.x + p.y * p.y) / (p.e - p.z);

        if (type == FermionType::FERMION)
        {
            if (ePlusPz < TINY * p.e)
            {
                const double norm = std::sqrt(2.0 * p.e);

                psi.plus[0] = 0.0;
                psi.plus[1] = norm;
                psi.minus[0] = -norm;
                psi.minus[1] = 0.0;
            }
            else 
            {
                const double sqrt_ePlusPz = std::sqrt(ePlusPz);
                const double norm = 1.0 / sqrt_ePlusPz;

                psi.plus[0] = sqrt_ePlusPz;
                psi.plus[1] = { p.x * norm, p.y * norm };
                psi.minus[0] = { -p.x * norm, p.y * norm };
                psi.minus[1] = sqrt_ePlusPz;
            }
        }
        else 
        {
            if (ePlusPz < TINY * p.e)
            {
                const double norm = std::sqrt(2.0 * p.e);

                psi.plus[0] = 0.0;
                psi.plus[1] = -norm;
                psi.minus[0] = norm;
                psi.minus[1] = 0.0;
            }
            else 
            {
                const double sqrt_ePlusPz = std::sqrt(ePlusPz);
                const double norm = 1.0 / sqrt_ePlusPz;

                psi.plus[0] = -sqrt_ePlusPz;
                psi.plus[1] = { -p.x * norm, -p.y * norm };
                psi.minus[0] = { p.x * norm, -p.y * norm };
                psi.minus[1] = -sqrt_ePlusPz;
            }
        }

        return psi;
    }

    WeylSpinors bra(const FourVector& p, FermionType type)
    {
        WeylSpinors psi = { };

        const double ePlusPz = p.z > 0.0 ? p.e + p.z : (p.x * p.x + p.y * p.y) / (p.e - p.z);

        if (type == FermionType::FERMION)
        {
            if (ePlusPz < TINY * p.e)
            {
                const double norm = std::sqrt(2.0 * p.e);

                psi.plus[0] = 0.0;
                psi.plus[1] = norm;
                psi.minus[0] = -norm;
                psi.minus[1] = 0.0;
            }
            else 
            {
                const double sqrt_ePlusPz = std::sqrt(ePlusPz);
                const double norm = 1.0 / sqrt_ePlusPz;

                psi.plus[0] = sqrt_ePlusPz;
                psi.plus[1] = { p.x * norm, -p.y * norm };
                psi.minus[0] = { -p.x * norm, -p.y * norm };
                psi.minus[1] = sqrt_ePlusPz;
            }
        }
        else 
        {
            if (ePlusPz < TINY * p.e)
            {
                const double norm = std::sqrt(2.0 * p.e);

                psi.plus[0] = 0.0;
                psi.plus[1] = -norm;
                psi.minus[0] = norm;
                psi.minus[1] = 0.0;
            }
            else 
            {
                const double sqrt_ePlusPz = std::sqrt(ePlusPz);
                const double norm = 1.0 / sqrt_ePlusPz;

                psi.plus[0] = -sqrt_ePlusPz;
                psi.plus[1] = { -p.x * norm, p.y * norm };
                psi.minus[0] = { p.x * norm, p.y * norm };
                psi.minus[1] = -sqrt_ePlusPz;
            }
        }

        return psi;
    }

    FourVectorC braGammaKet(const WeylSpinors& bra, const WeylSpinors& ket, int helicity)
    {
        powheg_assert(helicity == 1 || helicity == -1);
        const std::complex<double> h = static_cast<double>(helicity);

        const auto& b = (helicity == 1) ? bra.plus  : bra.minus;
        const auto& k = (helicity == 1) ? ket.plus  : ket.minus;

        return {
            b[0] * k[0] + k[1] * b[1],

            h * (b[0] * k[1] + k[0] * b[1]),

            h * I * (k[0] * b[1] - b[0] * k[1]),

            h * (b[0] * k[0] - k[1] * b[1])
        };
    }

    FourVectorC braGammaKetGluon(
        const WeylSpinors& bra, 
        const WeylSpinors& ket, 
        int helicity,
        const FourVector& pBra,
        const FourVector& pKet,
        const FourVector& pGluon,
        const FourVector& eps
    )
    {
        const FourVector pProp1 = pGluon + pBra;
        const FourVector pProp2 = pKet - pGluon;

        FourVectorC curr1 = {
            braSlashKet(bra, ket, helicity, { eps, pProp1, { 1.0,  0.0,  0.0,  0.0 } }),
            braSlashKet(bra, ket, helicity, { eps, pProp1, { 0.0, -1.0,  0.0,  0.0 } }),
            braSlashKet(bra, ket, helicity, { eps, pProp1, { 0.0,  0.0, -1.0,  0.0 } }),
            braSlashKet(bra, ket, helicity, { eps, pProp1, { 0.0,  0.0,  0.0, -1.0 } }),
        };

        FourVectorC curr2 = {
            braSlashKet(bra, ket, helicity, { { 1.0,  0.0,  0.0,  0.0 }, pProp2, eps }),
            braSlashKet(bra, ket, helicity, { { 0.0, -1.0,  0.0,  0.0 }, pProp2, eps }),
            braSlashKet(bra, ket, helicity, { { 0.0,  0.0, -1.0,  0.0 }, pProp2, eps }),
            braSlashKet(bra, ket, helicity, { { 0.0,  0.0,  0.0, -1.0 }, pProp2, eps }),
        };

        const double denom1 = 2.0 * dot(pGluon, pBra);
        const double denom2 = -2.0 * dot(pGluon, pKet);

        return curr1 / denom1 + curr2 / denom2;
    }
    
} // namespace powheg_dy
