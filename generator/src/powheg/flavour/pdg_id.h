#pragma once

#include "powheg/base.h"

namespace powheg
{
    enum PDG_ID : int
    {
        D = 1, DBAR = -1,
        U = 2, UBAR = -2,
        S = 3, SBAR = -3,
        C = 4, CBAR = -4,
        B = 5, BBAR = -5,
        T = 6, TBAR = -6,
        E_MINUS = 11, E_PLUS = -11,
        NU_E = 12, NU_E_BAR = -12,
        MU_MINUS = 13, MU_PLUS = -13,
        NU_MU = 14, NU_MU_BAR = -14,
        TAU_MINUS = 15, TAU_PLUS = -15,
        NU_TAU = 16, NU_TAU_BAR = -16,
        G = 21,
        GAMMA = 22,
        Z = 23,
        W_MINUS = -24, W_PLUS = 24
    };

    inline bool isAntiparticle(PDG_ID id1, PDG_ID id2)
    {
        return id1 == -id2;
    }

    inline bool isQuark(PDG_ID id)
    {
        return id != 0 && std::abs(id) <= 6;
    }

    inline bool isGluon(PDG_ID id)
    {
        return std::abs(id) == 21;
    }

    inline bool isColored(PDG_ID id)
    {
        return isGluon(id) || isQuark(id);
    }

} // namespace powheg
