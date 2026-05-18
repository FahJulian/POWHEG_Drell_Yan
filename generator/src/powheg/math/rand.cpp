#include "rand.h"

namespace powheg
{
    double rand()
    {
        static std::random_device dev;
        static std::mt19937 gen(dev());
        static std::uniform_real_distribution<double> dist(0.0, 1.0);

        return dist(gen);
    }

    double rand(double min, double max)
    {
        return rand() * (max - min) + min;
    }
    
} // namespace powheg
