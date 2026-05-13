#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/bbar/singular_region.h"
#include "powheg_dy/phase_space/born_phase_space.h"
#include "powheg_dy/phase_space/real_phase_space.h"

namespace powheg_dy
{
    enum class CollinearRemnantChannel
    {
        QQ_LEG1,        // q -> q g on leg 1
        QQ_LEG2,        // qbar -> qbar g on leg 2
        GLUON_LEG1,     // g -> q qbar, underlying Born quark on leg 1
        GLUON_LEG2      // g -> q qbar, underlying Born antiquark on leg 2
    };

    std::optional<CollinearRemnantChannel> remnantTypeFromRegion(
        const BornChannel& born,
        const RealChannel& real,
        const SingularRegion& region
    );

} // namespace powheg_dy
