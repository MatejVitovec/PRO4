#include <cmath>
#include "ComponentThermoVar.hpp"


ComponentThermoVar::ComponentThermoVar(const Vars<5>& in)
{
    ComponentThermoVar(in.getData());
}

ComponentThermoVar::ComponentThermoVar(const std::array<double, 5>& in) : ThermoVar(), compressibleData()
{
    std::copy(in.begin(), in.begin() + 3, data.begin());
    std::copy(in.begin() + 3, in.begin() + 5, compressibleData.getData().begin());
}


double ComponentThermoVar::density() const
{
    return compressibleData[0];
}

double ComponentThermoVar::internalEnergy() const
{
    return compressibleData[1];
}

void ComponentThermoVar::updateDensityAndInternalEnergy(double rho, double e)
{
    compressibleData[0] = rho;
    compressibleData[1] = e;
}