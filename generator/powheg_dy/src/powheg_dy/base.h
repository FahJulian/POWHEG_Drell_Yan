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
#include <charconv>
#include <stdexcept>

#include <LHAPDF/LHAPDF.h>

#define POWHEG_DEBUG

#undef assert

#ifdef POWHEG_DEBUG
# define assert2(expr, message) \
    if (!(expr)) \
    { \
        Log::err << message << std::endl; \
        throw std::runtime_error(""); \
    } 
# define assert1(expr) assert2(expr, #expr) 
#else 
#define assert2(expr, message)
#define assert1(expr)
#endif

#define GET_MACRO(_1, _2, name, ...) name
#define assert(...) GET_MACRO(__VA_ARGS__, assert2, assert1)(__VA_ARGS__) 

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
