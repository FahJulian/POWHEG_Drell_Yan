#include "math.h"

#include <math.h>

namespace powheg_dy
{
    int abs(int x) 
    { 
        return x > 0 ? x : -x; 
    }

    double abs(double x) 
    { 
        return x > 0.0 ? x : -x; 
    }

    int max(int x, int y) 
    { 
        return x > y ? x : y; 
    }

    double max(double x, double y) 
    { 
        return x > y ? x : y; 
    }

    double log(double x)
    {
        return std::log(x);
    }

    double exp(double x)
    {
        return std::exp(x);
    }

    double sqrt(double x)
    {
        return std::sqrt(x);
    }

    double pow(double x, int n)
    {
        return std::pow(x, n);
    }
    
    double sin(double x)
    {
        return std::sin(x);
    }

    double cos(double x)
    {
        return std::cos(x);
    }

    double tan(double x)
    {
        return std::tan(x);
    }

    double sinh(double x)
    {
        return std::sinh(x);
    }

    double cosh(double x)
    {
        return std::cosh(x);
    }

    double tanh(double x)
    {
        return std::tanh(x);
    }

    double asin(double x)
    {
        return std::asin(x);
    }

    double acos(double x)
    {
        return std::acos(x);
    }

    double atan(double x)
    {
        return std::atan(x);
    }

    double asinh(double x)
    {
        return std::asinh(x);
    }

    double acosh(double x)
    {
        return std::acosh(x);
    }

    double atanh(double x)
    {
        return std::atanh(x);
    }
    
} // namespace powheg_dy
