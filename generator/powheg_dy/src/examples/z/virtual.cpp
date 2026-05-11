#include "drell_yan.h"
#include "couplings.h"

namespace powheg_dy
{
    double DrellYanProcess::virtualOverBorn(const BornPhSpPt& born, double muR2) const
    {
        // double s = 2.0 * dot(born.p1Bar, born.p2Bar);
        double logMuS = std::log(muR2 / born.sHat);

        return m_config.C_F * (PI * PI - 8 - 3 * logMuS - logMuS * logMuS);
    }

} // namespace powheg_dy
