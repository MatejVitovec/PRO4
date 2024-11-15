#include <cmath>
#include <omp.h>

#include "Thermo.hpp"

void Thermo::updateThermo(VolField<ThermoVar>& thermoFiled, const VolField<Compressible>& w) const
{
    #pragma omp parallel for
    for (size_t i = 0; i < thermoFiled.size(); i++)
    {
        thermoFiled[i] = updateThermo(w[i], thermoFiled[i]);
    }

    #pragma omp parallel for
    for (size_t j = 0; j < thermoFiled.boundarySize(); j++)
    {
        for (size_t i = 0; i < thermoFiled.boundary(j).size(); i++)
        {
            thermoFiled.boundary(j)[i] = updateThermo(w.boundary(j)[i], thermoFiled.boundary(j)[i]);
        }        
    }
}

void Thermo::updateThermoInternal(VolField<ThermoVar>& thermoFiled, const VolField<Compressible>& w) const
{
    #pragma omp parallel for
    for (size_t i = 0; i < thermoFiled.size(); i++)
    {
        thermoFiled[i] = updateThermo(w[i], thermoFiled[i]);
    }
}


void Thermo::updateThermo(Field<ThermoVar>& thermoFiled, const Field<Compressible>& w) const
{
    #pragma omp parallel for
    for (size_t i = 0; i < thermoFiled.size(); i++)
    {
        thermoFiled[i] = updateThermo(w[i], thermoFiled[i]);
    }
}
