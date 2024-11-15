#ifndef THERMO_HPP
#define THERMO_HPP

#include "../Compressible.hpp"
#include "../ThermoVar.hpp"
#include "../Field.hpp"
#include "../VolField.hpp"
#include "../Mesh/Mesh.hpp"

#include <chrono>

class Thermo
{
    public:
    
        Thermo() {}

        virtual ~Thermo() {}

        void updateThermo(VolField<ThermoVar>& thermoFiled, const VolField<Compressible>& w) const;
        void updateThermoInternal(VolField<ThermoVar>& thermoFiled, const VolField<Compressible>& w) const;
        void updateThermo(Field<ThermoVar>& thermoFiled, const Field<Compressible>& w) const;

        virtual Vars<3> updateThermo(const Compressible& data, const ThermoVar& thermoData) const = 0;

        virtual Compressible primitiveToConservative(const Vars<5>& primitive) const = 0;
        virtual Compressible stagnationState(double TTot, double pTot) const = 0;

        virtual std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const = 0;

        virtual Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn, ThermoVar thermoIn) const = 0;
        
};

#endif // THERMO_HPP