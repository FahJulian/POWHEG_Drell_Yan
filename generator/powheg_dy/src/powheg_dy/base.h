#pragma once

#include <array>
#include <tuple>
#include <math.h>
#include <memory>
#include <string>
#include <vector>
#include <random>
#include <complex>
#include <fstream>
#include <sstream>
#include <iostream>

#include <LHAPDF/LHAPDF.h>

namespace powheg_dy
{
    namespace Log
    {
        class Logger
        {
        public:
            Logger(const char* prefix)
                : m_prefix(prefix)
            {
            }

            template<typename T>
            std::ostream& operator<<(const T& t) const
            {
                return std::cout << m_prefix << t;
            }

            template<typename T>
            void operator()(const T& t) const 
            {
                std::cout << m_prefix << t << std::endl;
            }

        private:
            std::string m_prefix = "";
        };

        static const Logger info("\033[32m[INFO]    ");
        static const Logger warn("\033[33m[WARNING] ");
        static const Logger err ("\033[31m[ERROR]   ");

    } // namespace Log

    static constexpr double PI = 3.14159265358979323846;

} // namespace powheg_dy
