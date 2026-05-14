#include "real_minus_ct_handler.h"

namespace powheg_dy
{
    double RealMinusCTHandler::dSigmaRealMinusCT(
        const BBarIntegrationPoint& point,
        const SingularRegion& region
    ) const
    {   
        if (region.type == SingularRegionType::SOFT)
            return dSigmaSoft(point, region);
        else if (region.type == SingularRegionType::ISR_LEG1)
            return dSigmaISRLeg1(point, region);
        else if (region.type == SingularRegionType::ISR_LEG2)
            return dSigmaISRLeg2(point, region);
        else if (region.type == SingularRegionType::FSR)
            throw std::runtime_error("FSR not yet implemented.");
        else
            throw std::runtime_error("Unknown singular region.");

        return 0.0;
    }   

    double RealMinusCTHandler::dSigmaSoft(
        const BBarIntegrationPoint& point,
        const SingularRegion& region
    ) const
    {
        return 0.0;
    }

    double RealMinusCTHandler::dSigmaISRLeg1(
        const BBarIntegrationPoint& point,
        const SingularRegion& region
    ) const
    {
        return 0.0;
    }

    double RealMinusCTHandler::dSigmaISRLeg2(
        const BBarIntegrationPoint& point,
        const SingularRegion& region
    ) const
    {
        return 0.0;
    }

} // namespace powheg_dy
