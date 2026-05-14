#include "process.h"

#include "powheg_dy/math/rand.h"
#include "powheg_dy/util/file.h"
#include "powheg_dy/input_parser.h"
#include "powheg_dy/emission/emission.h"
#include "powheg_dy/born_event/born_event_generator.h"
#include "powheg_dy/les_houches/les_houches_serializer.h"

#include <iostream>
#include <sstream>

namespace powheg_dy
{
    Process::~Process()
    {
        LHAPDF::setVerbosity(0);
        Log::info("Using LHAPDF. Please make sure to cite the paper:");
        Log::info("Eur.Phys.J. C75 (2015) 3, 132  (http://arxiv.org/abs/1412.7420)");
    }
    
    void Process::init(const std::string& configPath)
    {
        Timer::init();

        Log::info("Starting intitialization");

        try
        {
            m_config = InputParser(configPath).parse();
        }
        catch(const std::exception& e)
        {
            Log::err("Initialization failed.");
            return;
        }
        
        m_config.initLHAPDF();
        m_config.setDependentParams();
        m_config.extractLambdaFromPdf();

        m_bornPhSp = std::make_shared<BornPhaseSpace>(m_config);
        m_realPhSp = std::make_shared<FKSRealPhaseSpace>(m_config);
        m_bornGenerator = std::make_shared<BornEventGenerator>(*this, m_config, m_bornPhSp);
        m_emissionGenerator = std::make_shared<EmissionGenerator>(*this, m_config, m_realPhSp);
        m_bbarIntegrator = std::make_shared<BBarIntegrator>(*this, m_config, m_bornPhSp, m_realPhSp);

        Log::info << "Intitialization complete" << std::endl << std::endl;
        
        m_initialized = true;
    }

    void Process::run()
    {
        if (!m_initialized)
        {
            Log::err("Not intitialized, aborting run...");
            return;
        }

        try
        {
            clear();
            generateEvents();
            analyse();
        }
        catch(const std::runtime_error& e)
        {   
            Log::err << "Aborting..." << std::endl;

            if (std::string(e.what()) != "")
                Log::err << "Error message: " << e.what() << std::endl;
        }
    }

    void Process::writeToFile(const std::string& filePath) const
    {
        if (!m_initialized)
        {
            Log::err("Not initialized, can't write LHE file.");
            return;
        }

        Log::info << "Generating LHE file at " << filePath << std::endl;
        LesHouchesSerializer(*this, m_config).serialize(filePath);
        Log::info << "Done generating LHE file" << std::endl << std::endl;
    }

    void Process::clear()
    {
        m_events.clear();
        m_events.reserve(m_config.N_ACCEPTED_EVENTS);
        m_bbarIntegrator->clear();
    }

    void Process::generateEvents()  
    {
        m_bbarIntegrator->determineMaxWeight();

        double timeStart = Timer::getTime();

        Log::info("Starting event generation");
        while (m_events.size() < m_config.N_ACCEPTED_EVENTS)
        {   
            BornPhSpPt born = m_bbarIntegrator->sampleAccordingtoBTilde();

            if (m_config.NO_EMISSIONS)
            {
                Emission emission = Emission().reject();
                RealPhSpPt real = m_realPhSp->reconstruct(born, emission.rad);

                m_events.push_back({ born, real, emission });
            }
            else
            {
                Emission highestPtEm = Emission().reject();

                // Highest bid procedure, for Drell Yan only ISR radiation 
                for (int region = 1; region < 2; region++)
                {
                    Emission emission = m_emissionGenerator->generateEmission(born, region);
                    if (emission.kt2 > highestPtEm.kt2)
                        highestPtEm = emission;
                }
                
                RealPhSpPt real = m_realPhSp->reconstruct(born, highestPtEm.rad);
                m_events.push_back({ born, real, highestPtEm });
            }

            if (m_events.size() % 1000 == 0)
                Log::info << m_events.size() << " Events generated." << std::endl;
        }

        Log::info << "Duration: " << Timer::getTime() - timeStart << "s" << std::endl;
        Log::info << "Event generation done." << std::endl << std::endl;
    }

    void Process::analyse()
    {   
        if (!m_config.NO_EMISSIONS)
        {
            double rejected = 0.0;
            for (const auto& event : m_events)
            {
                if (event.emission.rejected)
                    rejected++;
            }
        
            Log::info << "No emission probability: " << rejected / static_cast<double>(m_events.size()) << std::endl;
        }

        Log::info << "Acceptance ratio: " << m_bbarIntegrator->getAcceptanceRatio() << std::endl;
        Log::info << "Total cross section: " << m_bbarIntegrator->getTotalCrossSection() << " pb." << std::endl << std::endl;
    }

} // namespace powheg_dy
