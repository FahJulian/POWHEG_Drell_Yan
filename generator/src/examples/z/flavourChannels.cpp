#include "drell_yan.h"

namespace powheg_dy
{
    std::vector<BornChannel> DrellYanProcess::bornChannels() const
    {
        return {
            { -5,  5, 5, { 11, -11 } },
            { -4,  4, 4, { 11, -11 } },
            { -3,  3, 3, { 11, -11 } },
            { -2,  2, 2, { 11, -11 } },
            { -1,  1, 1, { 11, -11 } },
            {  1, -1, 1, { 11, -11 } },
            {  2, -2, 2, { 11, -11 } },
            {  3, -3, 3, { 11, -11 } },
            {  4, -4, 4, { 11, -11 } },
            {  5, -5, 5, { 11, -11 } },
        };
    }

    std::vector<RealChannel> DrellYanProcess::realChannels(const BornChannel& bornChannel) const
    {   
        const int outID1 = bornChannel.outIDs[0];
        const int outID2 = bornChannel.outIDs[1];
        
        return {
            { bornChannel.id1, bornChannel.id2, { outID1, outID2, 21 } },                  
            { bornChannel.id1, 21,              { outID1, outID2, bornChannel.id1 } },      
            { 21,              bornChannel.id2, { outID1, outID2, bornChannel.id2 } },
        };
    }

} // namespace powheg_dy