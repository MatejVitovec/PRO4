#include "MeanPressureOutlet.hpp"

double MeanPressureOutlet::calculateCorrectionConstant(const std::vector<int>& ownerIndexList, const std::vector<Face>& faces, const Field<Compressible>& w) const
{
    double pressureSum = 0.0;
    double pressureSumMachLessOne = 0.0;

    for (auto & faceIndex : boundary.facesIndex)
    {
        pressureSum += w[ownerIndexList[faceIndex]].pressure();

        if(w[ownerIndexList[faceIndex]].normalVelocity(faces[faceIndex].normalVector)/w[ownerIndexList[faceIndex]].soundSpeed() < 1.0)
        {
            pressureSumMachLessOne += w[ownerIndexList[faceIndex]].pressure();
        }        
    }

    return (pressure*w.size() - pressureSum)/pressureSumMachLessOne + 1.0;
}

void MeanPressureOutlet::apply(const std::vector<int>& ownerIndexList, const std::vector<Face>& faces, const Field<Compressible>& w, Field<Compressible>& wr, const Thermo * const thermoModel) const
{
    double pressureCorrection = calculateCorrectionConstant(ownerIndexList, faces, w);

    for (auto & faceIndex : boundary.facesIndex)
    {
        if (w[ownerIndexList[faceIndex]].normalVelocity(faces[faceIndex].normalVector)/w[ownerIndexList[faceIndex]].soundSpeed() < 1.0)
        {
            wr[faceIndex] = thermoModel->primitiveToConservative(Vars<5>({w[ownerIndexList[faceIndex]].density(),
                                                                          w[ownerIndexList[faceIndex]].velocityU(),
                                                                          w[ownerIndexList[faceIndex]].velocityV(),
                                                                          w[ownerIndexList[faceIndex]].velocityW(),
                                                                          w[ownerIndexList[faceIndex]].pressure()*pressureCorrection}));
        }
    }
}