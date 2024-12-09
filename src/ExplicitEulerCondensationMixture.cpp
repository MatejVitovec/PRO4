#include "ExplicitEulerCondensationMixture.hpp"
#include <iostream>
#include "outputCFD.hpp"

#include "VolField.hpp"


void ExplicitEulerCondensationMixture::solve()
{
    init();

    mixture.updateVaporThermo(vaporThermoField, liquidThermoField, w);
    thermo->updateThermoInternal(vaporThermoField);
    mixture.updateMixtureThermo(mixtureThermoField, vaporThermoField, liquidThermoField, w);

    iter = 0;
    bool exitLoop = false;
    Vars<9> resNorm;

    auto startAll = std::chrono::high_resolution_clock::now();

    while (iter < maxIter && !exitLoop)
    {
        iter++;

        //boundField();

        updateTimeStep();

        calcBoundaryConditionFields();

        interpolateToFaces();

        calculateFluxes();

        Field<Vars<9>> res = calculateResidual();

        w += res*timeSteps;

        //UPDATE T liquid TODO

        liquidThermo.updateThermoFromT(liquidThermoField);
        mixture.updateVaporThermo(vaporThermoField, liquidThermoField, w);
        thermo->updateThermo(vaporThermoField);
        mixture.updateMixtureThermo(mixtureThermoField, vaporThermoField, liquidThermoField, w);


        if(iter % 100 == 0)
        {
            outputCFD::saveValue(savePath + "/time.txt", std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startAll).count());

            resNorm = res.norm();
            outputCFD::saveResidual(savePath + "/residuals.txt", resNorm[0]);
            std::cout << "iter: " << iter << " density res: " << resNorm[0] << std::endl;
            if(resNorm[0] < targetError) exitLoop = true;

        }

        if(iter % saveEveryIter == 0)
        {
            outputCFD::outputVTK(savePath + "/results/results." + std::to_string(iter) + ".vtk", mesh, w, mixtureThermoField);
        }
    }

    outputCFD::outputVTK(savePath + "/results/results." + std::to_string(iter) + ".vtk", mesh, w, mixtureThermoField);
    std::cout << "iter: " << iter << std::endl;
}