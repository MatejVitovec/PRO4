#ifndef INTERPOLATIONEOS
#define INTERPOLATIONEOS


#include "ThermoEoS.hpp"
#include "../Interpolation/Interpolation.hpp"
#include "../StateEquations/NonLinearSolver/NewtonMethod.hpp"

//TODO predalat na libovolny EOS (ThermoEoS), zatim by melo fungovat pro basic helmholtz

template <typename INTERPOLATION, typename EOS>
class InterpolationEoS : public ThermoEoS
{
    public:
        enum InterpolationType{BILINEAR, BIQUADRATIC};

        InterpolationEoS();

        Vars<3> updateThermo(const Compressible& data, const ThermoVar& thermoData) const;
        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;

        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn, ThermoVar thermoIn) const;
        std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const;


    private:
        EOS eos;

        INTERPOLATION pressureInterpolationFromRhoE;
        INTERPOLATION soundSpeedInterpolationFromRhoE;
        INTERPOLATION temperatureInterpolationFromRhoE;
        INTERPOLATION entropyInterpolationFromRhoE;

        INTERPOLATION energyInterpolationFromRhoPAux;

        NewtonMethod newtonSolver;
};

template <typename INTERPOLATION, typename EOS>
InterpolationEoS<INTERPOLATION, EOS>::InterpolationEoS() : Thermo(), Iapws95(), newtonSolver(NewtonMethod())
{
    std::vector<int> densityGridSize({200});
    std::vector<double> densityBoundary({1/1188.87, 1/0.008});
    Interpolation::Transformation densityTransformation(Interpolation::LOGINV);

    std::vector<int> internalEnergyGridSize({50, 100, 100});
    std::vector<double> internalEnergyBoundary({1800000.0, 2005000.0, 2650000.0, 4085000.27});
    Interpolation::Transformation internalEnergyTransformation(Interpolation::NONE);

    std::vector<int> pressureGridSize({800});
    std::vector<double> pressureBoundary({7000.0, 200000.0});
    //std::vector<double> pressureBoundary({400000.0, 6000000.0});
    Interpolation::Transformation pressureTransformation(Interpolation::LOGINV);

    pressureInterpolationFromRhoE = INTERPOLATION(
        densityGridSize, internalEnergyGridSize, densityBoundary, internalEnergyBoundary, densityTransformation, internalEnergyTransformation,
        [=](double density, double energy) { return eos.p(density, eos.tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });

    soundSpeedInterpolationFromRhoE = INTERPOLATION(
        densityGridSize, internalEnergyGridSize, densityBoundary, internalEnergyBoundary, densityTransformation, internalEnergyTransformation,
        [=](double density, double energy) { return eos.a(density, eos.tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });

    temperatureInterpolationFromRhoE = INTERPOLATION(
        densityGridSize, internalEnergyGridSize, densityBoundary, internalEnergyBoundary, densityTransformation, internalEnergyTransformation,
        [=](double density, double energy) { return eos.tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805); });

    energyInterpolationFromRhoPAux = INTERPOLATION(
        densityGridSize, pressureGridSize, densityBoundary, pressureBoundary, densityTransformation, pressureTransformation,
        [=](double density, double pressure) { return eos.e(density, eos.tFromRhoP(density, pressure, pressure/(density*461.51805))); });

    entropyInterpolationFromRhoE = INTERPOLATION(
        densityGridSize, internalEnergyGridSize, densityBoundary, internalEnergyBoundary, densityTransformation, internalEnergyTransformation,
        [=](double density, double energy) { return eos.s(density, eos.tFromRhoE(density, energy, (energy*(1.32 - 1.0))/461.51805)); });
}

template <typename INTERPOLATION, typename EOS>
Vars<3> InterpolationEoS<INTERPOLATION, EOS>::updateThermo(const Compressible& data, const ThermoVar& thermoData) const
{
    double pressure = pressureInterpolationFromRhoE.calcFastFind(data.density(), data.internalEnergy());
    double soundSpeed = soundSpeedInterpolationFromRhoE.calcFastFind(data.density(), data.internalEnergy());
    double temperature = temperatureInterpolationFromRhoE.calcFastFind(data.density(), data.internalEnergy());

    return Vars<3>({temperature, pressure, soundSpeed});
}

template <typename INTERPOLATION, typename EOS>
Compressible InterpolationEoS<INTERPOLATION, EOS>::primitiveToConservative(const Vars<5>& primitive) const
{
    double density = primitive[0];
    double pressure = primitive[4];

    double energyAux = energyInterpolationFromRhoPAux.calcFastFind(density, pressure);
    
    double energy = pressureInterpolationFromRhoE.calcInverseY(density, pressure, energyAux);
    double soundSpeed = soundSpeedInterpolationFromRhoE.calcFastFind(density, energy);

    return Compressible({density,
                         density*primitive[1],
                         density*primitive[2],
                         density*primitive[3],
                         density*(energy + 0.5*(primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3]))});
}

template <typename INTERPOLATION, typename EOS>
Compressible InterpolationEoS<INTERPOLATION, EOS>::stagnationState(double TTot, double pTot) const
{
    // FULL IAPWS95
    double rhoTot = eos.rhoFromTP(TTot, pTot, pTot/(eos.R()*TTot));

    return Compressible({rhoTot,
                         0.0,
                         0.0,
                         0.0,
                         rhoTot*(e(rhoTot, TTot))});
}

template <typename INTERPOLATION, typename EOS>
Compressible InterpolationEoS<INTERPOLATION, EOS>::isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn, ThermoVar thermoIn) const
{
    double pIn = std::min(thermoIn.pressure(), pTot);

    double guessRho = stateIn.density();
    double guess_e = stateIn.internalEnergy();

    std::pair<double, double> result = newtonSolver.solve([=](double val1, double val2) { return entropyInterpolationFromRhoE.calcFastFind(val1, val2)/sTot - 1; },
                                                          [=](double val1, double val2) { return entropyInterpolationFromRhoE.calcFastFindDiffX(val1, val2)/sTot; },
                                                          [=](double val1, double val2) { return entropyInterpolationFromRhoE.calcFastFindDiffY(val1, val2)/sTot; },
                                                          [=](double val1, double val2) { return pressureInterpolationFromRhoE.calcFastFind(val1, val2)/pIn - 1; },
                                                          [=](double val1, double val2) { return pressureInterpolationFromRhoE.calcFastFindDiffX(val1, val2)/pIn; },
                                                          [=](double val1, double val2) { return pressureInterpolationFromRhoE.calcFastFindDiffY(val1, val2)/pIn; },
                                                          guessRho,
                                                          guess_e);

    double rho = result.first;
    double e = result.second;

    double absU2 = std::max(2.0*hTot - 2.0*(e + pressureInterpolationFromRhoE.calcFastFind(rho, e)/rho), 0.0);
    double absU = std::sqrt(absU2);

    return Compressible({rho,
                         rho*absU*velocityDirection[0],
                         rho*absU*velocityDirection[1],
                         rho*absU*velocityDirection[2],
                         rho*(0.5*absU2 + e)});
}

template <typename INTERPOLATION, typename EOS>
std::array<double, 3> InterpolationEoS<INTERPOLATION, EOS>::initPressureTemperatureInlet(double pTot, double TTot) const
{
    double rhoTot = eos.rhoFromTP(TTot, pTot, pTot/(eos.R()*TTot));

    return std::array<double, 3>({rhoTot, eos.s(rhoTot, TTot), eos.h(rhoTot, TTot)});
}

#endif // INTERPOLATIONEOS