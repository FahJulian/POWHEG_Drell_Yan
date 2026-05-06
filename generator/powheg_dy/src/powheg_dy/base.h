#pragma once

#include <array>
#include <math.h>
#include <memory>
#include <string>
#include <vector>
#include <complex>
#include <sstream>
#include <iostream>

namespace powheg_dy
{
    namespace Log
    {
        class logger
        {
        public:
            logger(const char* prefix)
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

        // std::cout << "\033[31mThis is red\033[0m\n";
        // std::cout << "\033[32mThis is green\033[0m\n";
        // std::cout << "\033[33mThis is yellow\033[0m\n";
        // std::cout << "\033[34mThis is blue\033[0m\n";
        static const logger info("\033[32m[INFO]    ");
        static const logger warn("\033[33m[WARNING] ");
        static const logger err ("\033[31m[ERROR]   ");

    } // namespace Log

    static constexpr double PI = 3.14159265358979323846;

} // namespace powheg_dy
