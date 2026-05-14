#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/config.h"

namespace powheg_dy
{
    class InputParser
    {
    public:
        InputParser(const std::string& filePath)
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
