#pragma once

#include "powheg_dy/base.h"

namespace powheg
{
    class Config;

    double powhegAlphaS0(const Config& config, double qSq, int nF);
    double powhegAlphaS0customLambda(const Config& config, double qSq, int nF, double lambdaSq);
    double powhegAlphaS(const Config& config, double qSq);
    double powhegAlphaSCMW(const Config& config, double qSq);

} // namespace powheg_dy
