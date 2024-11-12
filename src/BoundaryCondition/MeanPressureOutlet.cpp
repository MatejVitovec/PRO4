#include "MeanPressureOutlet.hpp"
#include <iostream>

double MeanPressureOutlet::calculateCorrectionConstant(const Mesh& mesh, const Field<Compressible>& w) const
{
    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();

    double pressureSum = 0.0;
    double pressureSumMachLessOne = 0.0;

    for (auto & faceIndex : boundary.facesIndex)
    {
        pressureSum += w[ownerIndexList[faceIndex]].pressure();

        if(w[ownerIndexList[faceIndex]].normalVelocity(faceList[faceIndex].normalVector)/w[ownerIndexList[faceIndex]].soundSpeed() < 1.0)
        {
            pressureSumMachLessOne += w[ownerIndexList[faceIndex]].pressure();
        }        
    }

    return (pressure*boundary.facesIndex.size() - pressureSum)/pressureSumMachLessOne + 1.0;
}

Compressible MeanPressureOutlet::calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const
{
    return Compressible();
}

std::vector<Compressible> MeanPressureOutlet::calc(const Field<Compressible>& w, const Mesh& mesh, const Thermo * const thermoModel) const
{
    double pressureCorrection = calculateCorrectionConstant(mesh, w);

    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();

    std::vector<Compressible> out(boundary.facesIndex.size());

    for (int i = 0; i < boundary.facesIndex.size(); i++)
    {
        if (w[ownerIndexList[boundary.facesIndex[i]]].normalVelocity(faceList[boundary.facesIndex[i]].normalVector)/w[ownerIndexList[boundary.facesIndex[i]]].soundSpeed() < 1.0)
        {
            if (pressureCorrection > 0.5 && pressureCorrection < 1.5)
            {
                out[i] = thermoModel->primitiveToConservative(Vars<5>({w[ownerIndexList[boundary.facesIndex[i]]].density(),
                                                                       w[ownerIndexList[boundary.facesIndex[i]]].velocityU(),
                                                                       w[ownerIndexList[boundary.facesIndex[i]]].velocityV(),
                                                                       w[ownerIndexList[boundary.facesIndex[i]]].velocityW(),
                                                                       w[ownerIndexList[boundary.facesIndex[i]]].pressure()*pressureCorrection}));
            }
            else
            {
                out[i] = thermoModel->primitiveToConservative(Vars<5>({w[ownerIndexList[boundary.facesIndex[i]]].density(),
                                                                       w[ownerIndexList[boundary.facesIndex[i]]].velocityU(),
                                                                       w[ownerIndexList[boundary.facesIndex[i]]].velocityV(),
                                                                       w[ownerIndexList[boundary.facesIndex[i]]].velocityW(),
                                                                       pressure}));
            }
        }
        else
        {
            out[i] = w[ownerIndexList[boundary.facesIndex[i]]];
        }
    }
    
    return out;
}