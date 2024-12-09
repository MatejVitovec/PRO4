#include "IsentropicInlet.hpp"

Compressible IsentropicInlet::calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const
{
    return thermoModel->isentropicInlet(totalPressure, totalTemperature, totalDensity, totalEntropy, totalEnthalpy, velocityDirection, w, thermoVar);
}

CompressibleMixture IsentropicInlet::calculateState(const CompressibleMixture& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const
{
    return thermoModel->isentropicInlet(totalPressure, totalTemperature, totalDensity, totalEntropy, totalEnthalpy, velocityDirection, w, thermoVar);
}