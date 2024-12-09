#ifndef IDEALGASEOS_HPP
#define IDEALGASEOS_HPP

#include "ThermoEoS.hpp"

class IdealGasEoS : public ThermoEoS
{
    public:

        IdealGasEoS() : ThermoEoS(), gamma(1.3326), R(461.51805), energyShift(1995917.326) {}
        IdealGasEoS(double gamma_, double R_) : ThermoEoS(), gamma(gamma_), R(R_), energyShift(1995917.326) {}
        IdealGasEoS(double gamma_, double R_, double energyShift_) : ThermoEoS(), gamma(gamma_), R(R_), energyShift(energyShift_) {}

        Vars<3> updateThermo(const Compressible& data, const ThermoVar& thermoData) const;

        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;

        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn, ThermoVar thermoIn) const;
        std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const;

    private:
        double gamma;
        double R;
        double energyShift;

        double pressure(const Compressible& data) const;
        double soundSpeed(const Compressible& data) const;
        double temperature(const Compressible& data) const;

};

#endif // IDEALGASEOS_HPP