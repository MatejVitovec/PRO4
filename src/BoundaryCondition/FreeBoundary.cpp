#include "FreeBoundary.hpp"

Compressible FreeBoundary::calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const
{
    return w;
}

CompressibleMixture FreeBoundary::calculateState(const CompressibleMixture& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const
{
    return w;
}