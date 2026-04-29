#pragma once

#include "powheg_dy/base.h"

#include <string>
#include <sstream>

namespace powheg_dy
{
    class Process;

    class LesHouchesSerializer
    {
    public:
        LesHouchesSerializer(const Process& process)
            : m_process(process)
        {
        }

        void serialize(const std::string& filePath);

    private:
        const Process& m_process;
    };

} // namespace powheg_dy
