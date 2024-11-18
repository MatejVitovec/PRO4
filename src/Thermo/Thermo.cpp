#include <cmath>
#include <omp.h>

#include "Thermo.hpp"

void Thermo::updateThermo(VolField<ThermoVar>& thermoField, const VolField<Compressible>& w) const
{
    #pragma omp parallel for
    for (size_t i = 0; i < thermoField.size(); i++)
    {
        thermoField[i] = updateThermo(w[i], thermoField[i]);
    }

    #pragma omp parallel for
    for (size_t j = 0; j < thermoField.boundarySize(); j++)
    {
        for (size_t i = 0; i < thermoField.boundary(j).size(); i++)
        {
            thermoField.boundary(j)[i] = updateThermo(w.boundary(j)[i], thermoField.boundary(j)[i]);
        }        
    }
}

void Thermo::updateThermoInternal(VolField<ThermoVar>& thermoField, const VolField<Compressible>& w) const
{
    #pragma omp parallel for
    for (size_t i = 0; i < thermoField.size(); i++)
    {
        thermoField[i] = updateThermo(w[i], thermoField[i]);
    }
}


void Thermo::updateThermo(Field<ThermoVar>& thermoField, const Field<Compressible>& w) const
{
    #pragma omp parallel for
    for (size_t i = 0; i < thermoField.size(); i++)
    {
        thermoField[i] = updateThermo(w[i], thermoField[i]);
    }
}
