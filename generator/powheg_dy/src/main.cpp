#include "./powheg_dy/process.h"
#include "./powheg_dy/math/rand.h"
#include "./powheg_dy/phase_space/real_phase_space.h"

const std::string exportFilePath = "/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/data/events/events.lhe";
const std::string pdfSet = "NNPDF40MC_nlo_as_01180";

class DYProcess : public powheg_dy::Process
{   
public:
    inline int nBornLegs() const override { return 4; }
    inline int nRealLegs() const override { return 5; }
};

int main()
{
    DYProcess process;

    process.init(pdfSet);
    process.run();
    process.writeToFile(exportFilePath);
}
