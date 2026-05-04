#pragma once
 
#include "powheg_dy/event.h"

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
        void _writeEvent(const Event& event, std::stringstream& content) const;

    private:
        const Process& m_process;
    };

} // namespace powheg_dy
