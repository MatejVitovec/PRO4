#include "IapwsLiquidThermo.hpp"


void IapwsLiquidThermo::updateThermoFromT(ComponentThermoVar& thermoData) const
{
    double temperature= thermoData.temperature();
    if(temperature < 270) { return; }
    
    double rho = saturatedLiquidDensity(temperature);
    double e = saturatedLiquidInternalEnergy(temperature);
    thermoData.updateDensityAndInternalEnergy(rho, e);
}

void IapwsLiquidThermo::updateThermoFromT(Field<ComponentThermoVar>& thermoField) const
{
    #pragma omp parallel for
    for (size_t i = 0; i < thermoField.size(); i++)
    {
        updateThermoFromT(thermoField[i]);
    }
}


