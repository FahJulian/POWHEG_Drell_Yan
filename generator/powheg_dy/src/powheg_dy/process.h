#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/event.h"

#include <memory>
#include <vector>
#include <string>

#include <LHAPDF/LHAPDF.h>

namespace powheg_dy
{
    class Process
    {
    public:
        Process() = delete; 

        Process(double sqrtS, double mMin, double mMax)
            : m_sqrtS(sqrtS), m_mMin(mMin), m_mMax(mMax)
        {
        }

        ~Process() = default;

        void init(const std::string& pdfDataLocation, const std::string& pdfSet);
        void run();
        void writeToFile(const std::string& filePath) const;

        double getMMin() const { return m_mMin; }
        double getMMax() const { return m_mMax; }
        double getSqrtS() const { return m_sqrtS; }
        double getSigma() const { return m_totalCrossSection; }
        auto getEvents() const { return m_events; }
        const std::unique_ptr<LHAPDF::PDF>& getPdfs() const { return m_pdfs; }

    private:
        void _clear();
        void _determineMaxWeight();
        void _computeTotalCrossSection();

    private:
        const double m_sqrtS;   
        const double m_mMin, m_mMax;
        int m_nEventTrials = 0;
        double m_maxDSigma = 0.0;
        double m_totalCrossSection = 0.0;   // pb

        std::unique_ptr<LHAPDF::PDF> m_pdfs;
        std::vector<Event> m_events;
    };  

} // namespace powheg_dy
