#pragma once

namespace powheg_dy
{
    static constexpr double PI = 3.14159265358979323846;

    int abs(int x);
    double abs(double x);

    double log(double x);
    double exp(double x);

    double sqrt(double x);
    double pow(double x, int n);
    
    double sin(double x);
    double cos(double x);
    double tan(double x);

    double sinh(double x);
    double cosh(double x);
    double tanh(double x);

    double asin(double x);
    double acos(double x);
    double atan(double x);

    double asinh(double x);
    double acosh(double x);
    double atanh(double x);
    
} // namespace powheg_dy
