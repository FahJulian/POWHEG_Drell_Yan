#pragma once
 
#include "powheg_dy/base.h"
#include "powheg_dy/event.h"
#include "powheg_dy/config/config.h"

namespace powheg
{
    class BaseProcess; 

    class LesHouchesSerializer
    {
    public:
        LesHouchesSerializer(const BaseProcess& process, const Config& config)
            : m_process(process), m_config(config)
        {
        }

        void serialize(const std::string& filePath);

    private:
        void writeEvent(const Event& event, std::ostream& content) const;
        void writeEventHeader(std::ostream& content, int sign, int nParticles, double scalup) const;
        void writeEventBorn(const Event& event, std::ostream& content) const;
        void writeEventqqbar(const Event& event, std::ostream& content) const;
        void writeEventGluonLeg1(const Event& event, std::ostream& content) const;
        void writeEventGluonLeg2(const Event& event, std::ostream& content) const;

    private:
        const BaseProcess& m_process;
        const Config& m_config;
    };

} // namespace powheg_dy
