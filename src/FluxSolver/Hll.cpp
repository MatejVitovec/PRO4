#include "Hll.hpp"

Vars<5> Hll::claculateFlux(const Compressible& wl, const Compressible& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const
{
    enum {sl, ss, sr};
    Vars<3> wSpeed = waveSpeedsEstimate(wl, wr, thermoL, thermoR, normalVector);

    if (0 <= wSpeed[sl])
    {
        //FL
        return wl.flux(thermoL, normalVector);
    }
    else if(0 < wSpeed[sr])
    {
        //F*
        return ((wSpeed[sr]*wl.flux(thermoL, normalVector) - wSpeed[sl]*wr.flux(thermoR, normalVector) + wSpeed[sr]*wSpeed[sl]*(wr - wl)) / (wSpeed[sr] - wSpeed[sl]));
    }
    else
    {
        //FR
        return wr.flux(thermoR, normalVector);
    }
}

Vars<3> Hll::waveSpeedsEstimate(const Compressible& wl, const Compressible& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const
{
    double ul = wl.normalVelocity(normalVector);
    double ur = wr.normalVelocity(normalVector);
    double al = thermoL.soundSpeed();
    double ar = thermoR.soundSpeed();
    double rhol = wl.density();
    double rhor = wr.density();

    double sl = std::min(ul - al, ur - ar);
    double sr = std::max(ul + al, ur + ar);
    double ss = (thermoR.pressure() - thermoL.pressure() + rhol*ul*(sl - ul) - rhor*ur*(sr - ur))/(rhol*sl - rhol*ul - rhor*sr + rhor*ur);

    return Vars<3>({sl, ss, sr});
}