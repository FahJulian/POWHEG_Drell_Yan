#pragma once

#include <array>
#include <tuple>
#include <chrono>
#include <math.h>
#include <memory>
#include <string>
#include <vector>
#include <random>
#include <complex>
#include <fstream>
#include <sstream>
#include <optional>
#include <iostream>
#include <charconv>
#include <stdexcept>
#include <functional>

#include <LHAPDF/LHAPDF.h>

#undef powheg_assert

#ifdef POWHEG_DEBUG
# define powheg_assert2(expr, message) \
    if (!(expr)) \
    { \
        Log::err << message << Log::endl; \
        throw std::runtime_error(""); \
    } 
# define powheg_assert1(expr) powheg_assert2(expr, #expr) 
#else 
#define powheg_assert2(expr, message)
#define powheg_assert1(expr)
#endif

#define GET_MACRO(_1, _2, name, ...) name
#define powheg_assert(...) GET_MACRO(__VA_ARGS__, powheg_assert2, powheg_assert1)(__VA_ARGS__) 

namespace powheg
{
    class Log
    {
    private:
        Log(const char* prefix)
            : m_prefix(prefix)
        {
        }

    public:
        template<typename T>
        std::ostream& operator<<(const T& t) const
        {
            return std::cout << m_prefix << t;
        }

        template<typename T>
        void operator()(const T& t) const 
        {
            std::cout << m_prefix << t << Log::endl;
        }

    private:
        std::string m_prefix = "";

    public:
        static const Log info;
        static const Log warn;
        static const Log err; 
        static const std::string endl;
    };

    class Timer
    {
    public:
        static void init();
        static double getTime();

    private:
        static std::chrono::system_clock::time_point initTime;
    };

    static constexpr double PI = 3.14159265358979323846;

} // namespace powheg_dy
