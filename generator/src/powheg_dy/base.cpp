#include "base.h"

namespace powheg_dy
{
    const Log Log::info = ("\033[32m[INFO]    ");
    const Log Log::warn = ("\033[33m[WARNING] ");
    const Log Log::err  = ("\033[31m[ERROR]   ");
    const std::string Log::endl = "\033[0m\n";

    std::chrono::system_clock::time_point Timer::initTime = std::chrono::system_clock::time_point();

    void Timer::init()
    {
        initTime = std::chrono::system_clock::now();
    }

    double Timer::getTime()
    {
        std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
        std::chrono::duration<double> duration = now - initTime;

        return duration.count();
    }

} // namespace powheg_dy
