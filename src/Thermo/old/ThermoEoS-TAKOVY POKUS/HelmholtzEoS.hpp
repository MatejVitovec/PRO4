#ifndef HELMHOLTZEOS
#define HELMHOLTZEOS

#include "ThermoEoS.hpp"

template <typename EOS>
class HelmholtzEoS : public ThermoEoS
{
    public:

        HelmholtzEoS() : ThermoEoS() {}

        virtual ~HelmholtzEoS() {}

        Vars<3> updateThermo(const Compressible& data, const ThermoVar& thermoData) const;
        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;
        
        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn, ThermoVar thermoIn) const;
        std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const;

        double calcM2is(double p0, double T0, double p2) const;

        private:
            EOS helmholtz;
        
};

template <typename EOS>
Vars<3> HelmholtzEoS<EOS>::updateThermo(const Compressible& data, const ThermoVar& thermoData) const
{
    double rho = data.density();

    double T = EOS.tFromRhoE(rho, data.internalEnergy(), thermoData.temperature());

    return Vars<3>({T, EOS.p(rho, T), EOS.a(rho, T)});
}

template <typename EOS>
Compressible HelmholtzEoS<EOS>::primitiveToConservative(const Vars<5>& primitive) const
{
    double rho = primitive[0];
    double p = primitive[4];
    double T = EOS.tFromRhoP(rho, p, p/(EOS.R()*rho)); //odhad T pomoci idealniho plynu

    return Compressible({rho,
                         rho*primitive[1],
                         rho*primitive[2],
                         rho*primitive[3],
                         rho*(EOS.e(rho, T) + 0.5*(primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3]))});
}

template <typename EOS>
Compressible HelmholtzEoS<EOS>::stagnationState(double TTot, double pTot) const
{
    double rhoTot = EOS.rhoFromTP(TTot, pTot, pTot/(EOS.R()*TTot));

    return Compressible({rhoTot,
                         0.0,
                         0.0,
                         0.0,
                         rhoTot*(EOS.e(rhoTot, TTot))});
}

template <typename EOS>
double HelmholtzEoS<EOS>::calcM2is(double p0, double T0, double p2) const
{
    double rho0 = EOS.rhoFromTP(T0, p0, p0/(EOS.R()*T0));
    double s0 = EOS.s(rho0, T0);

    std::pair<double, double> state2 = EOS.RhoTFromSP(s0, p2, rho0, T0*1.2);
    double h0 = EOS.h(rho0, T0);
    double h2 = EOS.h(state2.first, state2.second);

    return std::sqrt((2*(h0 - h2))/std::pow(EOS.a(state2.first, state2.second), 2));
}

template <typename EOS>
Compressible HelmholtzEoS<EOS>::isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn, ThermoVar thermoIn) const
{
    thermoIn = updateThermo(stateIn, thermoIn);

    double pIn = std::min(thermoIn.pressure(), pTot);

    double guessRho = stateIn.density();
    double guessT = thermoIn.temperature();

    std::pair<double, double> result = EOS.RhoTFromSP(sTot, pIn, guessRho, guessT);

    double rho = result.first;
    double T = result.second;

    double absU2 = std::max(2.0*hTot - 2.0*EOS.h(rho, T), 0.0);
    double absU = std::sqrt(absU2);

    return Compressible({rho,
                         rho*absU*velocityDirection[0],
                         rho*absU*velocityDirection[1],
                         rho*absU*velocityDirection[2],
                         rho*(0.5*absU2 + EOS.e(rho, T))});
}

template <typename EOS>
std::array<double, 3> HelmholtzEoS<EOS>::initPressureTemperatureInlet(double pTot, double TTot) const
{
    double rhoTot = EOS.rhoFromTP(TTot, pTot, pTot/(EOS.R()*TTot));

    return std::array<double, 3>({rhoTot, EOS.s(rhoTot, TTot), EOS.h(rhoTot, TTot)});
}

#endif // HELMHOLTZEOS