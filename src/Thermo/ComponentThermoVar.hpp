#ifndef ComponentThermoVar_HPP
#define ComponentThermoVar_HPP

#include "../ThermoVar.hpp"

class ComponentThermoVar : public ThermoVar
{
    public:

        ComponentThermoVar() : ThermoVar() {}
        ComponentThermoVar(const Vars<3>& in) : ThermoVar(in) {}
        ComponentThermoVar(const std::array<double, 3>& in) : ThermoVar(in) {}

        ComponentThermoVar(const Vars<5>& in);
        ComponentThermoVar(const std::array<double, 5>& in);

        virtual ~ComponentThermoVar() {}

        double density() const;
        double internalEnergy() const;

        void updateDensityAndInternalEnergy(double rho, double e);


    private:
        Vars<2> compressibleData;
};




#endif // ComponentThermoVar_HPP