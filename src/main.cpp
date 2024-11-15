#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <chrono>
#include <fenv.h>

#include "outputCFD.hpp"
#include "CaseSetter.hpp"


int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    std::string savePath;
    if (argc >1 )
    {
        savePath = argv[1];
        savePath = "../results/" + savePath;
    }
    else
    { 
        //savePath = "../resultOK/lowTempAir/idealGasHLLE";
        savePath = "../results/SE1050/1ord/lowTempAir/IdealGas2";
        //savePath = "../results/SE1050/2ord/lowTempAir/idealGas";
    }
    
    CaseSetter setter = CaseSetter();
    setter.loadSettingFile(savePath + "/setup.txt");

    std::unique_ptr<FVMScheme> mySolver = setter.createAndSetSolver();
    mySolver->setSavePath(savePath);

    outputCFD::outputVTK(savePath + "/results/results.0.vtk", mySolver->getMesh(), mySolver->getResults(),  mySolver->getResultsThermo());

    auto stop1 = std::chrono::high_resolution_clock::now();

    mySolver->solve();    
    auto stop2 = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - stop1).count() << " ms\n";

    outputCFD::saveFieldOnBoundary(savePath + "/pressure.txt", "wall", mySolver->getMesh(), mySolver->getResults(), mySolver->getResultsThermo());

    /*Field<Compressible> w = outputCFD::loadCompressibleFieldFromVTK(savePath + "/results/results.665000.vtk");
    w = mySolver->getThermoRef()->updateField(w);
    outputCFD::saveFieldOnBoundary(savePath + "/pressure.txt", "wall", mySolver->getMesh(), w);*/
    //outputCFD::outputVTK(savePath + "/testResult.vtk",  mySolver->getMesh(), w);
    //outputCFD::outputVTKPeriodicBoundary(savePath + "/periodicResult.vtk", mySolver->getMesh(), w, Vars<3>({0.0, 0.0551168, 0.0}));

    outputCFD::outputVTKPeriodicBoundary(savePath + "/periodicResult.vtk", mySolver->getMesh(), mySolver->getResults(), mySolver->getResultsThermo(), Vars<3>({0.0, 0.0551168, 0.0}));

    return 0;
}