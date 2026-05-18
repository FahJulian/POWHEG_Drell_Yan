#include "born_virtual_handler.h"

#include "powheg/process.h"

namespace powheg
{
    double BornVirtualHandler::dSigmaBorn(
        const BBarCache& cache,
        size_t channelIdx
    ) const
    {
        const ChannelCache& channel = cache.channels[channelIdx];

        return 
            cache.born.jacobian 
            * channel.luminosityBorn 
            * cache.fluxFactorBorn 
            * channel.amp2Born;
    }

    double BornVirtualHandler::dSigmaVirtual(
        const BBarCache& cache,
        size_t channelIdx
    ) const
    {
        const ChannelCache& channel = cache.channels[channelIdx];

        const double dSigmaVirtualFinite = 
            cache.born.jacobian
            * channel.luminosityBorn 
            * cache.fluxFactorBorn
            * channel.amp2Virt;
        
        // TODO: Generalize
        const double otherTwoTerms = cache.alphaS / (2.0 * PI) * m_config.C_F
            * (-2.0 * std::log(cache.muF2 / cache.born.sHat) - 1.0 / 3.0 * PI * PI)
            * channel.luminosityBorn * cache.fluxFactorBorn * cache.born.jacobian * channel.amp2Born;
            
        return dSigmaVirtualFinite + otherTwoTerms;
    }

} // namespace powheg
