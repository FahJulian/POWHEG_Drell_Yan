#include "born_virtual_handler.h"

#include "powheg_dy/process.h"

namespace powheg_dy
{
    double BornVirtualHandler::dSigmaBorn(
        const BBarIntegrationPoint& point 
    ) const
    {
        const double fluxFactor = 1.0 / (2.0 * point.born.sHat);
        const double luminosity = point.f1Born * point.f2Born;
        
        return luminosity * fluxFactor * point.born.jacobian * point.amp2Born;
    }

    double BornVirtualHandler::dSigmaVirtual(
        const BBarIntegrationPoint& point, 
        const double muR2
    ) const
    {
        const double fluxFactor = 1.0 / (2.0 * point.born.sHat);
        const double luminosity = point.f1Born * point.f2Born;
        const double virtAmp2 = m_process.virtAmp2(point.born, point.bornChannel, point.amp2Born, muR2);
        
        return luminosity * fluxFactor * point.born.jacobian * virtAmp2;
    }

} // namespace powheg_dy
