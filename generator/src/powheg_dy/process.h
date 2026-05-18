#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/event.h"
#include "powheg_dy/config/config.h"
#include "powheg_dy/config/config_parser.h"
#include "powheg_dy/bbar/bbar_integrator.h"
#include "powheg_dy/emission/emission_generator.h"
#include "powheg_dy/phase_space/real_phase_space.h"
#include "powheg_dy/phase_space/born_phase_space_point.h"

namespace powheg
{
    class BaseProcess
    {
    private:
        explicit BaseProcess(std::unique_ptr<Config> config)
            : m_config(std::move(config))
        {
        }

        template<typename ConfigType>
        ConfigType& getConfig() { return static_cast<ConfigType&>(*m_config); }
        
    protected:
        virtual void initConfig(ConfigParser& parser) const = 0;

    public:
        virtual ~BaseProcess();

        virtual double bornAmp2(const BornPhSpPt& born, const BornChannel& bornChannel) const = 0;
        virtual double virtAmp2(const BornPhSpPt& born, const BornChannel& bornChannel, const double amp2Born, const double muR2) const = 0;
        virtual double realAmp2(const RealPhSpPt& real, const RealChannel& realChannel, const double alphaS) const = 0;

        virtual std::vector<BornChannel> bornChannels() const = 0;
        virtual std::vector<RealChannel> realChannels(const BornChannel& bornChannel) const = 0;

        virtual BornPhSpPt sampleBorn() const = 0;
        
        void init(const std::string& configPath);
        void run();
        void writeToFile(const std::string& filePath) const;

        double getSigma() const { return m_bbarIntegrator->getTotalCrossSection(); }
        auto getEvents() const { return m_events; }

    private:
        void clear();
        void initBaseConfig(ConfigParser& parser);
        void analyse();
        void generateEvents();

    private:
        bool m_initialized = false;
        std::vector<Event> m_events;

        std::shared_ptr<ISRRealPhaseSpace> m_realPhSp;
        std::shared_ptr<BBarIntegrator> m_bbarIntegrator;
        std::shared_ptr<EmissionGenerator> m_emissionGenerator;

    protected:
        std::unique_ptr<Config> m_config;

        template<typename ConfigType>
        friend class Process;
    };  

    template<typename ConfigType>
    class Process : public BaseProcess
    {
    public:
        Process()
            : BaseProcess(std::make_unique<ConfigType>()),
            m_config(getConfig<ConfigType>())
        {
        }
    
    protected:
        ConfigType& m_config;
    };

} // namespace powheg_dy
