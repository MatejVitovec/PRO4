#include "HeunScheme.hpp"
#include <iostream>
#include "outputCFD.hpp"

void HeunScheme::solve()
{
    init();

    Field<Compressible> wOld = Field<Compressible>(w.size());

    w = thermo->updateField(w);

    iter = 0;

    bool exitLoop = false;

    Field<Compressible> wn = Field<Compressible>(w.size());
    Vars<5> resNorm;

    while (iter < maxIter && !exitLoop)
    {
        iter++;
        wOld = w;

        //boundField();

        updateTimeStep();

        //applyBoundaryConditions();
        calcBoundaryConditionFields();

        //calculateWlWr();
        interpolateToFaces();

        //reconstruct();

        calculateFluxes();

        Field<Vars<5>> k1 = calculateResidual();
        
        wn = w + (k1*timeSteps);

        w = thermo->updateField(wn);



        //applyBoundaryConditions();
        calcBoundaryConditionFields();

        //calculateWlWr();
        interpolateToFaces();

        //reconstruct();

        calculateFluxes();

        Field<Vars<5>> k2 = calculateResidual();

        wn = wOld + (k1 + k2)*(timeSteps/2.0);

        w = thermo->updateField(wn);

        if(iter % 100 == 0)
        {
            resNorm = (k1 + k2).norm();
            outputCFD::saveResidual(savePath + "/residuals.txt", resNorm);
            std::cout << "iter: " << iter << " density res: " << resNorm[0] << std::endl;
            if(resNorm[0] < targetError) exitLoop = true;
        }
        
        if(iter % saveEveryIter == 0)
        {
            outputCFD::outputVTK(savePath + "/results/results." + std::to_string(iter) + ".vtk", mesh, w);
        }
    }

    outputCFD::outputVTK("../results/results." + std::to_string(iter) + ".vtk", mesh, w);
    std::cout << "iter: " << iter << std::endl;

    std::cout << "time: " << time << std::endl;
}