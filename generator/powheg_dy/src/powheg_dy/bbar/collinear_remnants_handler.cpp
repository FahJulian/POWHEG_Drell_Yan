#include "collinear_remnants_handler.h"

#include "powheg_dy/process.h"

namespace powheg_dy
{
namespace 
{
} // namespace

    double CollinearRemnantsHandler::dSigmaCollinearRemnants(
        const BornPhSpPt& born,
        const CollinearRemnantChannel& channel,
        const double u
    ) const
    {
        const double dSigmaBornPartonic = 1.0 / (2.0 * born.sHat) * born.jacobian * born.amp2;
    }


    double CollinearRemnantsHandler::oneMinusZTimesKerelQQ(const double z)
    {
        // eq. (2.103) times (1 - z)
        return m_config.C_F * (1.0 + z * z);
    }

    double CollinearRemnantsHandler::oneMinusZTimesKerelQG(const double z)
    {
        // eq. (2.104) times (1 - z)
        return m_config.C_F * (1.0 + (1.0 - z) * (1.0 - z)) * (1.0 - z) / z;
    }

    double CollinearRemnantsHandler::oneMinusZTimesKerelGQ(const double z)
    {
        // eq. (2.105) times (1 - z)
        return m_config.T_F * (1.0 - 2.0 * z * (1.0 - z)) * (1.0 - z);
    }

    double CollinearRemnantsHandler::oneMinusZTimesKerelGG(const double z)
    {
        // eq. (2.106) times (1 - z)
        return 2.0 * m_config.C_A * (z + (1.0 - z) * (1.0 - z) / z + z * (1.0 - z) * (1.0 - z));
    }

} // namespace powheg_dy
