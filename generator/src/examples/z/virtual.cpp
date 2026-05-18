#include "drell_yan.h"

namespace powheg
{
    double DrellYanProcess::virtAmp2(
        const BornPhSpPt& born,
        const BornChannel& bornChannel, 
        const double amp2Born, 
        const double muR2
    ) const
    {
        double logMuS = std::log(muR2 / born.sHat);

        return m_config.alphaS(muR2) / 2.0 / PI * m_config.C_F 
            * (PI * PI - 8 - 3 * logMuS - logMuS * logMuS) * amp2Born;
    }

} // namespace powheg_dy
