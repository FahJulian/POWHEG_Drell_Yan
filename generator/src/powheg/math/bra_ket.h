#pragma once

#include "powheg/base.h"
#include "powheg/math/four_vector.h"

namespace powheg
{
    struct FourVectorC
    {
        FourVectorC()
            : x0(0.0), x1(0.0), x2(0.0), x3(0.0)
        {
        }

        FourVectorC(const FourVector& v)
            : x0(v.x0), x1(v.x1), x2(v.x2), x3(v.x3)
        {
        }

        FourVectorC(std::complex<double> xx0, std::complex<double> xx1, 
            std::complex<double> xx2, std::complex<double> xx3)
            : x0(xx0), x1(xx1), x2(xx2), x3(xx3)
        {
        }
        
        std::complex<double> x0, x1, x2, x3;
    };

    struct WeylSpinors
    {
        std::array<std::complex<double>, 2> minus;
        std::array<std::complex<double>, 2> plus;
    };

    enum class FermionType
    {
        FERMION,
        ANTIFERMION
    };

    inline FourVectorC operator*(std::complex<double> a, const FourVectorC& v)
    {
        return { a * v.x0, a * v.x1, a * v.x2, a * v.x3 };
    }

    inline FourVectorC operator*(const FourVectorC& v, std::complex<double> a)
    {
        return a * v;
    }

    inline FourVectorC operator/(const FourVectorC& v, std::complex<double> a)
    {
        return (1.0 / a) * v;
    }

    inline FourVectorC operator+(const FourVectorC& v1, const FourVectorC& v2)
    {
        return { v1.x0 + v2.x0, v1.x1 + v2.x1, v1.x2 + v2.x2, v1.x3 + v2.x3 }; 
    }

    inline FourVectorC operator-(const FourVectorC& v1, const FourVectorC& v2)
    {
        return { v1.x0 - v2.x0, v1.x1 - v2.x1, v1.x2 - v2.x2, v1.x3 - v2.x3 }; 
    }

    inline std::complex<double> dot(const FourVectorC& v1, const FourVectorC& v2)
    {
        return v1.x0 * v2.x0 - v1.x1 * v2.x1 - v1.x2 * v2.x2 - v1.x3 * v2.x3;
    }

    FourVector polVec(const FourVector& p, double mass, int pol);

    WeylSpinors ket(const FourVector& p, FermionType type);
    WeylSpinors bra(const FourVector& p, FermionType type);

    FourVectorC braGammaKet(
        const WeylSpinors& bra, 
        const WeylSpinors& ket, 
        int helicity
    );

    FourVectorC braGammaKetGluon(
        const WeylSpinors& bra, 
        const WeylSpinors& ket, 
        int helicity,
        const FourVector& pBra,
        const FourVector& pKet,
        const FourVector& pGluon,
        const FourVector& eps
    );
    
} // namespace powheg
