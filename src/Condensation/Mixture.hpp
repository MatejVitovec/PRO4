#ifndef MIXTURE_HPP
#define MIXTURE_HPP

#include "../VolField.hpp"

#include "CompressibleMixture.hpp"
#include "../Thermo/ComponentThermoVar.hpp"

class Mixture
{
    public:
    
        Mixture() {}

        virtual ~Mixture() {}

        void updateVaporThermo(VolField<ComponentThermoVar>& vaporThermoField, const VolField<ComponentThermoVar>& liquidThermoField, const VolField<CompressibleMixture>& w) const;
        void updateVaporThermoInternal(VolField<ComponentThermoVar>& vaporThermoField, const VolField<ComponentThermoVar>& liquidThermoField, const VolField<CompressibleMixture>& w) const;
        void updateVaporThermo(Field<ComponentThermoVar>& vaporThermoField, const Field<ComponentThermoVar>& liquidThermoField, const Field<CompressibleMixture>& w) const;

        void updateMixtureThermo(VolField<ThermoVar>& mixtureThermoField, const VolField<ComponentThermoVar>& vaporThermoField, const VolField<ComponentThermoVar>& liquidThermoField, const VolField<CompressibleMixture>& w) const;
        void updateVaporThermoInternal(VolField<ThermoVar>& mixtureThermoField, const VolField<ComponentThermoVar>& vaporThermoField, const VolField<ComponentThermoVar>& liquidThermoField, const VolField<CompressibleMixture>& w) const;
        void updateVaporThermo(VolField<ThermoVar>& mixtureThermoField, const VolField<ComponentThermoVar>& vaporThermoField, const VolField<ComponentThermoVar>& liquidThermoField, const VolField<CompressibleMixture>& w) const;

    private:
        //LiguidEoS liquidEoS
        //Saturation 
        
};

#endif // MIXTURE_HPP