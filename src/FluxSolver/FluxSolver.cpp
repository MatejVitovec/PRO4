#include <cmath>

#include "FluxSolver.hpp"

Field<Vars<5>> FluxSolver::calculateFluxes(const Field<Compressible>& wl, const Field<Compressible>& wr, const Field<ThermoVar>& thermoFieldL, const Field<ThermoVar>& thermoFieldR, const std::vector<Face>& faceList) const
{
    Field<Vars<5>> out(wl.size());

    #pragma omp parallel for
    for (int i = 0; i < wl.size(); i++)
    {
        out[i] = claculateFlux(wl[i], wr[i], thermoFieldL[i], thermoFieldR[i], faceList[i].normalVector)*faceList[i].area;
    }
    
    return out;
}


Field<Vars<9>> FluxSolver::calculateFluxes(const Field<CompressibleMixture>& wl, const Field<CompressibleMixture>& wr, const Field<ThermoVar>& thermoFieldL, const Field<ThermoVar>& thermoFieldR, const std::vector<Face>& faceList) const
{
    Field<Vars<9>> out(wl.size());

    #pragma omp parallel for
    for (int i = 0; i < wl.size(); i++)
    {
        out[i] = claculateFlux(wl[i], wr[i], thermoFieldL[i], thermoFieldR[i], faceList[i].normalVector)*faceList[i].area;
    }
    
    return out;
}
