#pragma once
 
#include "powheg_dy/base.h"
#include "powheg_dy/event.h"
#include "powheg_dy/config/config.h"

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
        void writeEvent(const Event& event, std::stringstream& content) const;
        void writeEventHeader(std::stringstream& content, int sign, int nParticles, double scalup) const;
        void writeEventBorn(const Event& event, std::stringstream& content) const;
        void writeEventqqbar(const Event& event, std::stringstream& content) const;
        void writeEventGluonLeg1(const Event& event, std::stringstream& content) const;
        void writeEventGluonLeg2(const Event& event, std::stringstream& content) const;

    private:
        const Process& m_process;
        const Config& m_config;
    };

} // namespace powheg_dy
