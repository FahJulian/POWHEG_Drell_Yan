#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/config/config.h"

namespace powheg_dy
{
    class ConfigParser
    {
    public:
        ConfigParser(const std::string& filePath)
            : m_filePath(filePath)
        {
        }

        Config parse() const;

    private:
        std::map<std::string, std::string> readParams() const;
        Config extractParams(std::map<std::string, std::string>& params) const;

    private:
        std::string m_filePath;
    }; 

} // namespace powheg_dy
