#include "four_vector.h"

namespace powheg_dy
{
    FourVector operator+(const FourVector& v1, const FourVector& v2)
    {
        return { v1.x0 + v2.x0, v1.x1 + v2.x1, v1.x2 + v2.x2, v1.x3 + v2.x3 };
    }

    ThreeVector operator+(const ThreeVector& v1, const ThreeVector& v2)
    {
        return { v1.x1 + v2.x1, v1.x2 + v2.x2, v1.x3 + v2.x3 };
    }

    FourVector operator-(const FourVector& v1, const FourVector& v2)
    {
        return { v1.x0 - v2.x0, v1.x1 - v2.x1, v1.x2 - v2.x2, v1.x3 - v2.x3 };
    }

    ThreeVector operator-(const ThreeVector& v1, const ThreeVector& v2)
    {
        return { v1.x1 - v2.x1, v1.x2 - v2.x2, v1.x3 - v2.x3 };
    }

    double dot(const FourVector& v1, const FourVector& v2)
    {
        return v1.x0 * v2.x0 - v1.x1 * v2.x1 - v1.x2 * v2.x2 - v1.x3 * v2.x3;
    }

    double dot(const ThreeVector& v1, const ThreeVector& v2)
    {
        return v1.x1 * v2.x1 + v1.x2 * v2.x2 + v1.x3 * v2.x3;
    }

    FourVector operator*(const FourVector& v, double a)
    {
        return { a * v.x0, a * v.x1, a * v.x2, a * v.x3 }; 
    }

    ThreeVector operator*(const ThreeVector& v, double a)
    {
        return { a * v.x1, a * v.x2, a * v.x3 };
    }

    FourVector operator*(double a, const FourVector& v)
    {
        return v * a;
    }

    ThreeVector operator*(double a, const ThreeVector& v)
    {
        return v * a;
    }

    FourVector operator/(const FourVector& v, double a)
    {
        return v * (1.0 / a);
    }

    ThreeVector operator/(const ThreeVector& v, double a)
    {
        return v * (1.0 / a);
    }

    FourVector FourVector::operator-() const
    {
        return { -x0, -x1, -x2, -x3 };
    }

    ThreeVector ThreeVector::operator-() const
    {
        return { -x1, -x2, -x3 };
    }
    
    ThreeVector FourVector::getThreeVec() const
    {
        return { x1, x2, x3 };
    }

    double FourVector::square() const 
    {
        return dot(*this, *this);
    }

    double FourVector::rapidity() const
    {
        return 0.5 * std::log((e + pZ) / (e - pZ));
    }

    double ThreeVector::square() const 
    {
        return dot(*this, *this);
    }

    ThreeVector FourVector::getBeta() const
    {
        return getThreeVec() / e;
    }

    FourVector FourVector::boost(const ThreeVector& beta) const
    {
        double betaSq = beta.square();
        double gamma = 1.0 / std::sqrt(1.0 - betaSq);

        if (betaSq == 0)
            return *this;

        double dotProd = dot(beta, getThreeVec());
        double factor = (gamma - 1.0) / betaSq * dotProd + gamma * x0;

        return {
            gamma * (e + dotProd),
            x1 + factor * beta.x1,
            x2 + factor * beta.x2,
            x3 + factor * beta.x3
        };
    }

} // namespace powheg_dy
