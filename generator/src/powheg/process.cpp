#include "process.h"

#include "powheg/emission/emission.h"
#include "powheg/les_houches/les_houches_serializer.h"

namespace powheg
{
    BaseProcess::~BaseProcess()
    {
        LHAPDF::setVerbosity(0);
        Log::info("Using LHAPDF. Please make sure to cite the paper:");
        Log::info("Eur.Phys.J. C75 (2015) 3, 132  (http://arxiv.org/abs/1412.7420)");
    }

    void BaseProcess::initBaseConfig(ConfigParser& parser)
    {
        powheg_assert(parser.extract("SQRT_S", m_config->SQRT_S), "No value for SQRT_S found.");
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
        const bool useBorn = parser.extract("BTILDE_USE_BORN", m_config->BTILDE_USE_BORN);
        const bool useVirt = parser.extract("BTILDE_USE_VIRT", m_config->BTILDE_USE_VIRT);
        const bool useReal = parser.extract("BTILDE_USE_REAL", m_config->BTILDE_USE_REAL);
        const bool useColl = parser.extract("BTILDE_USE_COLL", m_config->BTILDE_USE_COLL);

        if (parser.extract("BTILDE_BORNONLY", m_config->BTILDE_BORNONLY); m_config->BTILDE_BORNONLY)
        {
            if (useBorn && !m_config->BTILDE_USE_BORN) Log::warn("BTILDE_BORNONLY overwrite BTILDE_USE_BORN");
            if (useVirt && m_config->BTILDE_USE_VIRT) Log::warn("BTILDE_BORNONLY overwrite BTILDE_USE_VIRT");
            if (useReal && m_config->BTILDE_USE_REAL) Log::warn("BTILDE_BORNONLY overwrite BTILDE_USE_REAL");
            if (useColl && m_config->BTILDE_USE_COLL) Log::warn("BTILDE_BORNONLY overwrite BTILDE_USE_COLL");

            m_config->BTILDE_USE_BORN = true;
            m_config->BTILDE_USE_VIRT = false;
            m_config->BTILDE_USE_REAL = false;
            m_config->BTILDE_USE_COLL = false;
        }
        
        parser.extract("ALPHA_S_FROM_PDF", m_config->ALPHA_S_FROM_PDF);

        parser.extract("N_C", m_config->N_C);
        parser.extract("C_F", m_config->C_F);
        parser.extract("C_A", m_config->C_A);
        parser.extract("T_F", m_config->T_F);

        parser.extract("CMW_CHARM_TRSHLD_SQ", m_config->CMW_CHARM_TRSHLD_SQ);
        parser.extract("CMW_BOTTOM_TRSHLD_SQ", m_config->CMW_BOTTOM_TRSHLD_SQ);

        m_config->S = m_config->SQRT_S * m_config->SQRT_S;
        m_config->CMW_CHARM_TRSHLD_SQ  = m_config->CMW_CHARM_TRSHLD  * m_config->CMW_CHARM_TRSHLD;
        m_config->CMW_BOTTOM_TRSHLD_SQ = m_config->CMW_BOTTOM_TRSHLD * m_config->CMW_BOTTOM_TRSHLD;
    }
    
    void BaseProcess::init(const std::string& configPath)
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
        m_bbarIntegrator = std::make_shared<BBarIntegrator>(*this, *m_config, *m_realPhSp);

        m_channels = findChannels(bornChannels(), realChannels());
        
        Log::info << "Intitialization complete" << Log::endl << Log::endl;
        
        m_initialized = true;
    }

    void BaseProcess::run()
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

    void BaseProcess::writeToFile(const std::string& filePath) const
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

    void BaseProcess::clear()
    {
        m_events.clear();
        m_events.reserve(m_config->N_ACCEPTED_EVENTS);
        m_bbarIntegrator->clear();
    }

    void BaseProcess::generateEvents()  
    {
        m_bbarIntegrator->determineMaxWeight();

        double timeStart = Timer::getTime();

        Log::info("Starting event generation");
        while (m_events.size() < m_config->N_ACCEPTED_EVENTS)
        {   
            BornEvent bornEvent = m_bbarIntegrator->sampleAccordingtoBTilde();

            Emission emission = m_config->NO_EMISSIONS ? Emission().reject() : 
                m_emissionGenerator->generateEmission(bornEvent.phaseSpace, bornEvent.channel);
            RealPhSpPt real = m_realPhSp->reconstruct(bornEvent.phaseSpace, emission.rad);

            m_events.push_back({ 
                bornEvent.phaseSpace, 
                bornEvent.channel, 
                bornEvent.weightSign,
                real, 
                emission
            });

            if (m_events.size() % 1000 == 0)
                Log::info << m_events.size() << " Events generated." << Log::endl;
        }

        Log::info << "Duration: " << Timer::getTime() - timeStart << "s" << Log::endl;
        Log::info << "Event generation done." << Log::endl << Log::endl;
    }

    void BaseProcess::analyse()
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

} // namespace powheg
