#pragma once
 
#include "powheg_dy/base.h"
#include "powheg_dy/math/four_vector.h"

namespace powheg_dy
{
    struct BornChannel
    {
        int id1;
        int id2;

        std::vector<int> outIDs;
    };

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

} // namespace powheg_dy
