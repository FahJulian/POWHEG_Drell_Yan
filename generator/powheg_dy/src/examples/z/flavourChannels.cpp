#include "drell_yan.h"

namespace powheg_dy
{
namespace 
{
    static const std::vector<BornChannel> BORN_CHANNELS = {
        { -5,  5,  5 },
        { -4,  4,  4 },
        { -3,  3,  3 },
        { -2,  2,  2 },
        { -1,  1,  1 },
        {  1, -1,  1 },
        {  2, -2,  2 },
        {  3, -3,  3 },
        {  4, -4,  4 },
        {  5, -5,  5 },
    };

} // anonymous namespace

    std::vector<BornChannel> DrellYanProcess::bornChannels() const
    {
        return BORN_CHANNELS;
    }

} // namespace powheg_dy