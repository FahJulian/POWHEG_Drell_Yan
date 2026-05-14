#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/config/config.h"

namespace powheg_dy
{
    class ConfigParser
    {
    private:
        ConfigParser(const std::string& filePath)
            : m_filePath(filePath)
        {
        }

        void parse();
        void close();

    public:
        template<typename T>
        bool extract(const std::string& name, T& value)
        {
            if (auto pos = m_params.find(name); pos != m_params.end())
            {
                std::istringstream(pos->second) >> value;
                m_params.erase(pos);

                Log::info << "Found parameter " << name << " = " << value << Log::endl;

                return true;
            }
            
            return false;
        }

    private:
        std::string m_filePath;
        std::unordered_map<std::string, std::string> m_params;
        
        friend class Process;
    }; 

} // namespace powheg_dy
