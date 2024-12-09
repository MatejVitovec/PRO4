#ifndef IAPWSLIQUIDTHERMO_HPP
#define IAPWSLIQUIDTHERMO_HPP

#include "Thermo.hpp"
#include "Saturation/IapwsSaturation.hpp"

class IapwsLiquidThermo : public IapwsSaturation
{
    public:

        IapwsLiquidThermo() {}

        Vars<3> updateThermo(const Compressible& data, const ThermoVar& thermoData) const;
        
        void updateThermo(ComponentThermoVar& thermoData) const;
        void updateThermoFromT(ComponentThermoVar& thermoData) const;

        void updateThermoFromT(Field<ComponentThermoVar>& thermoField) const;

    private:

};

#endif // IAPWSLIQUIDTHERMO_HPP