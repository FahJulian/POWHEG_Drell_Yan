#include "config_parser.h"

#include "powheg_dy/util/file.h"

namespace powheg_dy
{
namespace
{   
    bool extractString(std::map<std::string, std::string>& params, const std::string& name, std::string& value)
    {
        if (auto pos = params.find(name); pos != params.end())
        {
            value = pos->second;  
            params.erase(pos);

            Log::info << "Found parameter " << name << " = " << value << std::endl;
            
            return true;
        }
        
        return false;
    }

    bool extractDouble(std::map<std::string, std::string>& params, const std::string& name, double& value)
    {
        if (auto pos = params.find(name); pos != params.end())
        {
            value = std::stod(pos->second);  
            params.erase(pos);

            Log::info << "Found parameter " << name << " = " << value << std::endl;
            
            return true;
        }
        
        return false;
    }

    template<typename T>
    bool extractInt(std::map<std::string, std::string>& params, const std::string& name, T& value)
    {
        if (auto pos = params.find(name); pos != params.end())
        {
            value = static_cast<T>(std::stol(pos->second));  
            params.erase(pos);

            Log::info << "Found parameter " << name << " = " << value << std::endl;

            return true;
        }
        
        return false;
    }

    bool extractBool(std::map<std::string, std::string>& params, const std::string& name, bool& value)
    {
        if (auto pos = params.find(name); pos != params.end())
        {
            if (pos->second == "1" || pos->second == "true")
            {
                value = true;
                params.erase(pos);
            }
            else if (pos->second == "0" || pos->second == "false")
            {
                value = false;
                params.erase(pos);
            }
            else
            {
                params.erase(pos);
                return false;
            }

            Log::info << "Found parameter " << name << " = " << value << std::endl;

            return true;
        }
        
        return false;
    }
}

    Config ConfigParser::parse() const
    {
        Log::info << "Reading config from file " << m_filePath << "." << std::endl;

        std::map<std::string, std::string> params = readParams();
        Config config = extractParams(params);

        for (auto pair : params)
            Log::warn << "Unknown parameter: " << pair.first << std::endl;

        Log::info << "Done reading config." << std::endl;

        return config;
    }

    Config ConfigParser::extractParams(std::map<std::string, std::string>& params) const
    {
        Config config = { };

        powheg_assert(extractDouble(params, "SQRT_S", config.SQRT_S), "No value for SQRT_S found.");
        powheg_assert(extractDouble(params, "M_MIN", config.M_MIN), "No value for M_MIN found.");
        powheg_assert(extractDouble(params, "M_MAX", config.M_MAX), "No value for M_MAX found.");
        powheg_assert(extractDouble(params, "PT_SQ_CUTOFF", config.PT_SQ_CUTOFF), "No value for PT_SQ_CUTOFF found.");
        powheg_assert(extractInt<size_t>(params, "N_ACCEPTED_EVENTS", config.N_ACCEPTED_EVENTS), "No value for N_ACCEPTED_EVENTS found.");
        powheg_assert(extractString(params, "PDF_NAME", config.PDF_NAME), "No value for PDF_NAME found.");

        bool bornVetoWeight = extractDouble(params, "BORN_VETO_WEIGHT", config.BORN_VETO_WEIGHT);
        bool nTrialEvents   = extractInt<size_t>(params, "N_TRIAL_EVENTS", config.N_TRIAL_EVENTS);

        if (!bornVetoWeight && !nTrialEvents)
            throw std::runtime_error("No value for both N_TRIAL_EVENTS and BORN_VETO_WEIGHT found.");
        else if (bornVetoWeight && nTrialEvents)
            Log::info("Both N_TRIAL_EVENTS and BORN_VETO_WEIGHT found. Ignoring N_TRIAL_EVENTS");

        extractBool(params, "NO_EMISSIONS", config.NO_EMISSIONS);
        extractBool(params, "BTILDE_BORNONLY", config.BTILDE_BORNONLY);
        
        extractBool(params, "ALPHA_S_FROM_PDF", config.ALPHA_S_FROM_PDF);

        return config;
    }

    std::map<std::string, std::string> ConfigParser::readParams() const
    {
        std::map<std::string, std::string> params;
        
        File file(m_filePath);

        std::string line;
        size_t lineNr = 0;
        while (file.getLine(line))
        {
            lineNr++;

            if (size_t pos = line.find('#'); pos != std::string::npos)
                line.erase(pos);

            std::istringstream iss(line);

            std::string name;
            std::string value;
            std::string eqSign;
            std::string extra;

            if (!(iss >> name))
                continue;

            if (!(iss >> eqSign) || eqSign != "=")
            {
                Log::warn << "Line " << lineNr << ": Invalid Syntax." << std::endl;
                continue;
            }

            if (!(iss >> value))
            {
                Log::warn << "Line " << lineNr << ": Missing value for " << name << "." << std::endl;
                continue;
            }

            if (iss >> extra)
            {
                Log::warn << "Line " << lineNr << ": Too many values for " << name << "." << std::endl;
                continue;
            }

            params.emplace(name, value);
        }

        return params;
    }

} // namespace powheg_dy
