#include "process.h"

#include "powheg_dy/emission/emission.h"
#include "powheg_dy/les_houches/les_houches_serializer.h"

namespace powheg_dy
{
    Process::~Process()
    {
        LHAPDF::setVerbosity(0);
        Log::info("Using LHAPDF. Please make sure to cite the paper:");
        Log::info("Eur.Phys.J. C75 (2015) 3, 132  (http://arxiv.org/abs/1412.7420)");
    }

    void Process::initBaseConfig(ConfigParser& parser)
    {
        powheg_assert(parser.extract("SQRT_S", m_config->SQRT_S), "No value for SQRT_S found.");
        powheg_assert(parser.extract("M_MIN", m_config->M_MIN), "No value for M_MIN found.");
        powheg_assert(parser.extract("M_MAX", m_config->M_MAX), "No value for M_MAX found.");
        powheg_assert(parser.extract("PT_SQ_CUTOFF", m_config->PT_SQ_CUTOFF), "No value for PT_SQ_CUTOFF found.");
        powheg_assert(parser.extract("N_ACCEPTED_EVENTS", m_config->N_ACCEPTED_EVENTS), "No value for N_ACCEPTED_EVENTS found.");
        powheg_assert(parser.extract("PDF_NAME", m_config->PDF_NAME), "No value for PDF_NAME found.");

        bool bornVetoWeight = parser.extract("BORN_VETO_WEIGHT", m_config->BORN_VETO_WEIGHT);
        bool nTrialEvents   = parser.extract("N_TRIAL_EVENTS", m_config->N_TRIAL_EVENTS);

        if (!bornVetoWeight && !nTrialEvents)
            throw std::runtime_error("No value for both N_TRIAL_EVENTS and BORN_VETO_WEIGHT found.");
        else if (bornVetoWeight && nTrialEvents)
            Log::info("Both N_TRIAL_EVENTS and BORN_VETO_WEIGHT found. Ignoring N_TRIAL_EVENTS");

        parser.extract("NO_EMISSIONS", m_config->NO_EMISSIONS);
        parser.extract("BTILDE_BORNONLY", m_config->BTILDE_BORNONLY);
        
        parser.extract("ALPHA_S_FROM_PDF", m_config->ALPHA_S_FROM_PDF);

        parser.extract("N_C", m_config->N_C);
        parser.extract("C_F", m_config->C_F);
        parser.extract("C_A", m_config->C_A);
        parser.extract("T_F", m_config->T_F);

        parser.extract("CMW_CHARM_TRSHLD_SQ", m_config->CMW_CHARM_TRSHLD_SQ);
        parser.extract("CMW_BOTTOM_TRSHLD_SQ", m_config->CMW_BOTTOM_TRSHLD_SQ);

        m_config->S = m_config->SQRT_S * m_config->SQRT_S;
        m_config->CMW_CHARM_TRSHLD_SQ = m_config->CMW_CHARM_TRSHLD * m_config->CMW_CHARM_TRSHLD;
        m_config->CMW_BOTTOM_TRSHLD_SQ = m_config->CMW_BOTTOM_TRSHLD * m_config->CMW_BOTTOM_TRSHLD;
    }
    
    void Process::init(const std::string& configPath)
    {
        Timer::init();

        Log::info("Starting intitialization");

        try
        {
            ConfigParser parser(configPath);
            parser.parse();

            initBaseConfig(parser);
            m_config->initPdf();
            m_config->extractLambdaFromPdf();

            initConfig(parser);

            parser.close();
        }
        catch(const std::exception& e)
        {
            if (std::string(e.what()) != "")
                Log::err << e.what() << Log::endl;

            Log::err("Initialization failed.");
            return;
        }

        m_realPhSp = std::make_shared<ISRRealPhaseSpace>(*m_config);
        m_emissionGenerator = std::make_shared<EmissionGenerator>(*this, *m_config, m_realPhSp);
        m_bbarIntegrator = std::make_shared<BBarIntegrator>(*this, *m_config, m_realPhSp);

        Log::info << "Intitialization complete" << Log::endl << Log::endl;
        
        m_initialized = true;
    }

    void Process::run()
    {
        if (!m_initialized)
        {
            Log::err("Not initialized, can't run.");
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
            if (std::string(e.what()) != "")
                Log::err << e.what() << Log::endl;

            Log::err << "Aborting." << Log::endl;
        }
    }

    void Process::writeToFile(const std::string& filePath) const
    {
        if (!m_initialized)
        {
            Log::err("Not initialized, can't write LHE file.");
            return;
        }

        Log::info << "Generating LHE file at " << filePath << Log::endl;
        LesHouchesSerializer(*this, *m_config).serialize(filePath);
        Log::info << "Done generating LHE file" << Log::endl << Log::endl;
    }

    void Process::clear()
    {
        m_events.clear();
        m_events.reserve(m_config->N_ACCEPTED_EVENTS);
        m_bbarIntegrator->clear();
    }

    void Process::generateEvents()  
    {
        m_bbarIntegrator->determineMaxWeight();

        double timeStart = Timer::getTime();

        Log::info("Starting event generation");
        while (m_events.size() < m_config->N_ACCEPTED_EVENTS)
        {   
            BornEvent bornEvent = m_bbarIntegrator->sampleAccordingtoBTilde();

            if (m_config->NO_EMISSIONS)
            {
                Emission emission = Emission().reject();
                RealPhSpPt real = m_realPhSp->reconstruct(bornEvent.phaseSpace, emission.rad);

                m_events.push_back({ 
                    bornEvent.phaseSpace, 
                    bornEvent.channel, 
                    bornEvent.weightSign,
                    real, 
                    emission 
                });
            }
            else
            {
                Emission highestPtEm = Emission().reject();

                // Highest bid procedure, for Drell Yan only ISR radiation 
                for (int region = 1; region < 2; region++)
                {
                    Emission emission = m_emissionGenerator->generateEmission(bornEvent.phaseSpace, bornEvent.channel, bornEvent.amp2Born, region);
                    if (emission.kt2 > highestPtEm.kt2)
                        highestPtEm = emission;
                }
                
                RealPhSpPt real = m_realPhSp->reconstruct(bornEvent.phaseSpace, highestPtEm.rad);

                m_events.push_back({ 
                    bornEvent.phaseSpace, 
                    bornEvent.channel, 
                    bornEvent.weightSign,
                    real, 
                    highestPtEm
                });
            }

            if (m_events.size() % 1000 == 0)
                Log::info << m_events.size() << " Events generated." << Log::endl;
        }

        Log::info << "Duration: " << Timer::getTime() - timeStart << "s" << Log::endl;
        Log::info << "Event generation done." << Log::endl << Log::endl;
    }

    void Process::analyse()
    {   
        if (!m_config->NO_EMISSIONS)
        {
            double rejected = 0.0;
            for (const auto& event : m_events)
            {
                if (event.emission.rejected)
                    rejected++;
            }
        
            Log::info << "No emission probability: " << rejected / static_cast<double>(m_events.size()) << Log::endl;
        }

        Log::info << "Acceptance ratio: " << m_bbarIntegrator->getAcceptanceRatio() << Log::endl;
        Log::info << "Total cross section: " << m_bbarIntegrator->getTotalCrossSection() << " pb." << Log::endl << Log::endl;
    }

} // namespace powheg_dy
