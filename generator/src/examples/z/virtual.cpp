#include "drell_yan.h"
#include "couplings.h"

#include "powheg_dy/alpha_s.h"

namespace powheg_dy
{
    double DrellYanProcess::virtualContribution(const BornPhSpPt& born, double muR2) const
    {
        double logMuS = std::log(muR2 / born.sHat);

        return alphaS(m_config, muR2) / 2.0 / PI * m_config.C_F 
            * (PI * PI - 8 - 3 * logMuS - logMuS * logMuS) * bornContribution(born);
    }

    double DrellYanProcess::virtAmp2(
        const BornPhSpPt& born, 
        const double amp2Born, 
        const double muR2
    ) const
    {
        double logMuS = std::log(muR2 / born.sHat);

        return alphaS(m_config, muR2) / 2.0 / PI * m_config.C_F 
            * (PI * PI - 8 - 3 * logMuS - logMuS * logMuS) * amp2Born;
    }

} // namespace powheg_dy
