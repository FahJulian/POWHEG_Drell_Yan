#pragma once
 
#include "powheg/base.h"
#include "powheg/math/four_vector.h"

namespace powheg
{
    struct BornPhSpPt
    {
        double sHat;
        double x1Bar;
        double x2Bar;

        FourVector p1Bar;
        FourVector p2Bar;
        std::vector<FourVector> pOut;

        double jacobian;
    };

} // namespace powheg
