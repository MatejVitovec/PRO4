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

/*double saturatedDensity(double T)
{
    double theta = 1.0-(T/647.096);
    std::array<double, 6> b = {1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -674694.450};
    return (1 + b[0]*std::pow(theta, 1.0/3.0) + b[1]*std::pow(theta, 2.0/3.0) + b[2]*std::pow(theta, 5.0/3.0) + b[3]*std::pow(theta, 16.0/3.0) + b[4]*std::pow(theta, 43.0/3.0) + b[5]*std::pow(theta, 110.0/3.0))/322.0 ;
}

double saturatedPressure(double T)
{
    double theta = 1.0-(T/647.096);
    std::array<double, 6> a = {-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502};
    return std::exp((647.096/T)*(a[0]*theta + a[1]*std::pow(theta, 1.5) + a[2]*std::pow(theta, 3) + a[3]*std::pow(theta, 3.5) + a[4]*std::pow(theta, 4) + a[5]*std::pow(theta, 7.5)))*22064000;
}*/


int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    //////////////////////////////////
    /*std::shared_ptr<Helmholtz> thermo = std::make_shared<Iapws95>();


    std::vector<double> T = std::vector<double>(601);
    std::vector<double> p = std::vector<double>(601);
    std::vector<std::vector<double>> rho = std::vector<std::vector<double>>(601);
    for (size_t i = 0; i < T.size(); i++)
    {
        T[i] = 373 + i*0.37833333;
    }

    double first = log(10);
    double last = log(10000000);
    double dx = std::abs(first - last)/600;

    for (size_t i = 0; i < p.size(); i++)
    {
        p[i] = exp(first + i*dx);
    }

    for (size_t j = 0; j < T.size(); j++)
    {
        rho[j] = std::vector<double>(601);
        rho[j][0] = thermo->rhoFromTP(T[j], p[0], p[0]/(461.51805*T[j]));        
        for (size_t i = 1; i < p.size(); i++)
        {
            if (saturatedPressure(T[j]) < p[i])
            {
                break;
            }
            
            rho[j][i] = thermo->rhoFromTP(T[j], p[i], rho[j][i-1]);
        }
        std::cout << j << "temp ok" << std::endl;
    }

    std::ofstream f;
    f.open("../compFactor.txt", std::ios_base::app);

    for (size_t j = 0; j < T.size(); j++)
    {
        for (size_t i = 0; i < p.size(); i++)
        {
            f << T[j] << " " << p[i] << " " << p[i]/(rho[j][i]*T[j]*461.51805) << " " << rho[j][i] << std::endl;
        }
    }
    f.close();
    
    return 0;*/

    //////////////////////////////////

    /*std::shared_ptr<Iapws95Thermo> thermo = std::make_shared<Iapws95Thermo>();
    for (int i = 0; i < 100; i++)
    {
        std::cout << "M2is: " << thermo->calcM2is(5000000.0, 573.15,  2196600 + i*1) << " p2: " <<  2196600 + i*1 << std::endl;
    }*/

    std::string savePath;
    if (argc >1 )
    {
        savePath = argv[1];
        savePath = "../results/" + savePath;
    }
    else
    {
        //savePath = "../resultOK/lowTempAir/idealGasHLLE";
        savePath = "../results/test2";
        //savePath = "../results/SE1050/2ord/lowTempAir/idealGas";
    }
    
    CaseSetter setter = CaseSetter();
    setter.loadSettingFile(savePath + "/setup.txt");

    std::unique_ptr<FVMScheme> mySolver = setter.createAndSetSolver();
    mySolver->setSavePath(savePath);

    outputCFD::outputVTK(savePath + "/results/results.0.vtk", mySolver->getMesh(), mySolver->getResults());

    auto stop1 = std::chrono::high_resolution_clock::now();

    mySolver->solve();    
    auto stop2 = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - stop1).count() << " ms\n";

    outputCFD::saveFieldOnBoundary(savePath + "/pressure.txt", "wall", mySolver->getMesh(), mySolver->getResults());

    /*Field<Compressible> w = outputCFD::loadCompressibleFieldFromVTK(savePath + "/results/results.665000.vtk");
    w = mySolver->getThermoRef()->updateField(w);
    outputCFD::saveFieldOnBoundary(savePath + "/pressure.txt", "wall", mySolver->getMesh(), w);*/
    //outputCFD::outputVTK(savePath + "/testResult.vtk",  mySolver->getMesh(), w);
    //outputCFD::outputVTKPeriodicBoundary(savePath + "/periodicResult.vtk", mySolver->getMesh(), w, Vars<3>({0.0, 0.0551168, 0.0}));

    outputCFD::outputVTKPeriodicBoundary(savePath + "/periodicResult.vtk", mySolver->getMesh(), mySolver->getResults(), Vars<3>({0.0, 0.0551168, 0.0}));

    return 0;
}