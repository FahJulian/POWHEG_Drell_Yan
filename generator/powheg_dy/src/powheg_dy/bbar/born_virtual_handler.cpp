#include "born_virtual_handler.h"

#include "powheg_dy/process.h"

namespace powheg_dy
{
    double BornVirtualHandler::dSigmaBorn(
        const BornPhSpPt& born 
    ) const
    {
        const double fluxFactor = 1.0 / (2.0 * born.sHat);
        const double luminosity = born.f1 * born.f2;
        
        return luminosity * fluxFactor * born.jacobian * born.amp2;
    }

    double BornVirtualHandler::dSigmaVirtual(
        const BornPhSpPt& born, 
        const double muR2
    ) const
    {
        const double fluxFactor = 1.0 / (2.0 * born.sHat);
        const double luminosity = born.f1 * born.f2;
        const double virtAmp2 = m_process.virtAmp2(born, muR2);
        
        return luminosity * fluxFactor * born.jacobian * virtAmp2;
    }

} // namespace powheg_dy
