#pragma once
 
#include "powheg_dy/base.h"
#include "powheg_dy/event.h"
#include "powheg_dy/config.h"

namespace powheg_dy
{
    class Process; 

    class LesHouchesSerializer
    {
    public:
        LesHouchesSerializer(const Process& process, const Config& config)
            : m_process(process), m_config(config)
        {
        }

        void serialize(const std::string& filePath);

    private:
        void _writeEvent(const Event& event, std::stringstream& content) const;

    private:
        const Process& m_process;
        const Config& m_config;
    };

} // namespace powheg_dy
