#pragma once

#include "powheg/base.h"
#include "powheg/flavour/channel.h"



/// OLD ///
#include "powheg/flavour/born_channel.h"
#include "powheg/phase_space/real_phase_space.h"
#include "powheg/phase_space/born_phase_space_point.h"

namespace powheg
{
    struct RadiationCache
    {
        double xi = 0.0;
        double y = 0.0;
        double phi = 0.0;
        
        double xiTilde = 0.0;
        double xiMax = 0.0;
        double xiMaxLeg1 = 0.0;
        double xiMaxLeg2 = 0.0;

        double xiLeg1 = 0.0;
        double xiLeg2 = 0.0;

        double zLeg1 = 0.0;
        double zLeg2 = 0.0;

        double dzduLeg1 = 0.0;
        double dzduLeg2 = 0.0;

        double jacobianRad = 0.0;
        double jacobianRadSoftOverXiSq = 0.0;
        double jacobianRadLeg1 = 0.0;
        double jacobianRadLeg2 = 0.0;

        RealPhSpPt real;
    };

    struct ChannelCache
    {
        Channel channel;

        double amp2Born = 0.0;
        double amp2Virt = 0.0;

        double f1Born = 0.0;
        double f2Born = 0.0;

        double luminosityBorn = 0.0;

        // TODO: Put real amplitude and real luminosities here
        std::vector<double> f1Coll;
        std::vector<double> f2Coll;
    };

    struct BBarCache
    {
        BornPhSpPt born;

        double muF2 = 0.0;
        double muR2 = 0.0;
        double alphaS = 0.0;

        double fluxFactorBorn = 0.0;

        RadiationCache rad;

        std::vector<ChannelCache> channels;
    };

    // TODO: Split into CollRemnIntPoint and RealIntPoint
    struct BBarCacheOld
    {
        BornPhSpPt born;
        BornChannel bornChannel;
        RealPhSpPt real;

        double xi = 0.0;
        double y = 0.0;
        double phi = 0.0;

        double xiTilde = 0.0;
        double xiLeg1 = 0.0;
        double xiLeg2 = 0.0;

        double xiMax = 0.0;
        double xiMaxLeg1 = 0.0;
        double xiMaxLeg2 = 0.0;

        double jacobian = 0.0;
        double jacobianOverXiSqSoft = 0.0;
        double jacobianLeg1 = 0.0;
        double jacobianLeg2 = 0.0;

        double dzduLeg1 = 0.0;
        double dzduLeg2 = 0.0;

        double zLeg1 = 0.0;
        double zLeg2 = 0.0;

        double f1Born = 0.0;
        double f2Born = 0.0;

        double f1RealQ = 0.0;
        double f2RealQ = 0.0;

        double f1RealG = 0.0;
        double f2RealG = 0.0;

        double alphaS = 0.0;

        double amp2Born = 0.0;
    };

} // namespace powheg