#pragma once

#include "powheg_dy/base.h"

namespace powheg_dy
{
    struct ThreeVector
    {
        union 
        {
            double x1;
            double x;
            double pX;
        };
        
        union 
        {
            double x2;
            double y;
            double pY;
        };

        union 
        {
            double x3;
            double z;
            double pZ;
        };

        ThreeVector operator-() const;
        double square() const;
    };

    struct FourVector
    {
        union
        {
            double x0;
            double t;
            double e;
        };

        union 
        {
            double x1; 
            double x;
            double pX;
        };

        union 
        {
            double x2;
            double y;
            double pY;
        };
        
        union
        {
            double x3;
            double z;
            double pZ;
        };

        FourVector()
            : x0(0.0), x1(0.0), x2(0.0), x3(0.0)
        {
        }

        FourVector(double x0, double x1, double x2, double x3)
            : x0(x0), x1(x1), x2(x2), x3(x3)
        {
        }

        FourVector(double x0, ThreeVector x)
            : x0(x0), x1(x.x1), x2(x.x2), x3(x.x3)
        {
        }

        FourVector operator-() const;

        ThreeVector getThreeVec() const;
        double square() const;
        FourVector boost(const FourVector& delta) const;
    };

    FourVector operator+(const FourVector& v1, const FourVector& v2);
    ThreeVector operator+(const ThreeVector& v1, const ThreeVector& v2);
    FourVector operator-(const FourVector& v1, const FourVector& v2);
    ThreeVector operator-(const ThreeVector& v1, const ThreeVector& v2);
    double operator*(const FourVector& v1, const FourVector& v2);
    double operator*(const ThreeVector& v1, const ThreeVector& v2);
    FourVector operator*(const FourVector& v, double a);
    ThreeVector operator*(const ThreeVector& v, double a);
    FourVector operator*(double a, const FourVector& v);
    ThreeVector operator*(double a, const ThreeVector& v);
    FourVector operator/(const FourVector& v, double a);
    ThreeVector operator/(const ThreeVector& v, double a);
    
} // namespace powheg_dy
