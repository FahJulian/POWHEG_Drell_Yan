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

        FourVector(double xx0, double xx1, double xx2, double xx3)
            : x0(xx0), x1(xx1), x2(xx2), x3(xx3)
        {
        }

        FourVector(double xx0, ThreeVector xx)
            : x0(xx0), x1(xx.x1), x2(xx.x2), x3(xx.x3)
        {
        }

        FourVector operator-() const;

        ThreeVector getThreeVec() const;
        ThreeVector getBeta() const;
        double rapidity() const;
        double square() const;
        FourVector boost(const ThreeVector& beta) const;
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
