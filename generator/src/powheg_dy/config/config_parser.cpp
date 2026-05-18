#include "config_parser.h"

namespace powheg
{
    void ConfigParser::close()
    {
        for (const auto& pair : m_params)
            Log::warn << "Unknown parameter: " << pair.first << Log::endl;
    }

    void ConfigParser::parse()
    {
        Log::info << "Reading config from file " << m_filePath << "." << Log::endl;

        std::ifstream file(m_filePath);

        if (!file)
            throw std::runtime_error("Could not open config file " + m_filePath);

        std::string line;
        size_t lineNr = 0;
        while (std::getline(file, line))
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
                Log::warn << "Line " << lineNr << ": Invalid Syntax." << Log::endl;
                continue;
            }

            if (!(iss >> value))
            {
                Log::warn << "Line " << lineNr << ": Missing value for " << name << "." << Log::endl;
                continue;
            }

            if (iss >> extra)
            {
                Log::warn << "Line " << lineNr << ": Too many values for " << name << "." << Log::endl;
                continue;
            }

            m_params.emplace(name, value);
        }

        Log::info << "Done reading config." << Log::endl;
    }

} // namespace powheg_dy
