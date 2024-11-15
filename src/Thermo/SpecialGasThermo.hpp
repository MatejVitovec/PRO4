#ifndef SPECIALGASTHERMO
#define SPECIALGASTHERMO

#include "Thermo.hpp"
#include "StateEquations/SpecialGas.hpp"

class SpecialGasThermo : public Thermo, SpecialGas
{
    public:

        SpecialGasThermo() : Thermo(), SpecialGas() {}

        Vars<3> updateThermo(const Compressible& data, const ThermoVar& thermoData) const;
        
        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;
        
        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn, ThermoVar thermoIn) const;
        std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const;
};

#endif // SPECIALGASTHERMO