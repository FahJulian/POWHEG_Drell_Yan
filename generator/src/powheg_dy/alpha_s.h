#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/config.h"

namespace powheg_dy
{
    double alphaS0(const Config& config, double qSq, int nF);
    double alphaS0customLambda(const Config& config, double qSq, int nF, double lambdaSq);
    double alphaS(const Config& config, double qSq);
    double alphaSCMW(const Config& config, double qSq);

} // namespace powheg_dy
