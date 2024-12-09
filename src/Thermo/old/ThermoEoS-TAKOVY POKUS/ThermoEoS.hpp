#ifndef THERMOEOS_HPP
#define THERMOEOS_HPP

#include "../Compressible.hpp"
#include "../ThermoVar.hpp"
#include "../Field.hpp"
#include "../VolField.hpp"
#include "../Mesh/Mesh.hpp"
#include "../ComponentThermoVar.hpp"

#include <chrono>

class ThermoEoS
{
    public:
    
        ThermoEoS() {}

        virtual ~ThermoEoS() {}

        virtual Vars<3> updateThermo(const Compressible& data, const ThermoVar& thermoData) const = 0;
        virtual void updateThermo(ComponentThermoVar& theromVar);

        virtual Compressible primitiveToConservative(const Vars<5>& primitive) const = 0;
        virtual Compressible stagnationState(double TTot, double pTot) const = 0;

        virtual std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const = 0;
        virtual Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn, ThermoVar thermoeosIn) const = 0;
        
};

#endif // THERMOEOS_HPP