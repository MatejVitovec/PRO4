#include <omp.h>

#include "Mixture.hpp"

void Mixture::updateVaporThermo(VolField<ComponentThermoVar>& vaporThermoField, const VolField<ComponentThermoVar>& liquidThermoField, const VolField<CompressibleMixture>& w) const
{
    #pragma omp parallel for
    for (size_t i = 0; i < vaporThermoField.size(); i++)
    {
        double rhoM = w[i].density();
        double eM = w[i].internalEnergy();
        double alpha = w[i].alpha();

        double rhoL = liquidThermoField[i].density();
        double eL = liquidThermoField[i].internalEnergy();
        
        double rhoV = (rhoM - alpha*rhoL);
        double eV = (rhoM*eM - alpha*rhoL*eL)/(rhoM - rhoL*alpha);

        vaporThermoField[i].updateDensityAndInternalEnergy(rhoV, eV);
    }

    #pragma omp parallel for
    for (size_t j = 0; j < vaporThermoField.boundarySize(); j++)
    {
        for (size_t i = 0; i < vaporThermoField.boundary(j).size(); i++)
        {
            double rhoM = w.boundary(j)[i].density();
            double eM = w.boundary(j)[i].internalEnergy();
            double alpha = w.boundary(j)[i].alpha();

            double rhoL = liquidThermoField.boundary(j)[i].density();
            double eL = liquidThermoField.boundary(j)[i].internalEnergy();
            

            double rhoV = (rhoM - alpha*rhoL);
            double eV = (rhoM*eM - alpha*rhoL*eL)/(rhoM - rhoL*alpha);

            vaporThermoField.boundary(j)[i].updateDensityAndInternalEnergy(rhoV, eV);
        }        
    }
}

void Mixture::updateVaporThermoInternal(VolField<ComponentThermoVar>& vaporThermoField, const VolField<ComponentThermoVar>& liquidThermoField, const VolField<CompressibleMixture>& w) const
{
    #pragma omp parallel for
    for (size_t i = 0; i < vaporThermoField.size(); i++)
    {
        double rhoM = w[i].density();
        double eM = w[i].internalEnergy();
        double alpha = w[i].alpha();

        double rhoL = liquidThermoField[i].density();
        double eL = liquidThermoField[i].internalEnergy();
        
        double rhoV = (rhoM - alpha*rhoL);
        double eV = (rhoM*eM - alpha*rhoL*eL)/(rhoM - rhoL*alpha);

        vaporThermoField[i].updateDensityAndInternalEnergy(rhoV, eV);
    }
}

void Mixture::updateVaporThermo(Field<ComponentThermoVar>& vaporThermoField, const Field<ComponentThermoVar>& liquidThermoField, const Field<CompressibleMixture>& w) const
{
    #pragma omp parallel for
    for (size_t i = 0; i < vaporThermoField.size(); i++)
    {
        double rhoM = w[i].density();
        double eM = w[i].internalEnergy();
        double alpha = w[i].alpha();

        double rhoL = liquidThermoField[i].density();
        double eL = liquidThermoField[i].internalEnergy();
        
        double rhoV = (rhoM - alpha*rhoL);
        double eV = (rhoM*eM - alpha*rhoL*eL)/(rhoM - rhoL*alpha);

        vaporThermoField[i].updateDensityAndInternalEnergy(rhoV, eV);
    }
}



void Mixture::updateMixtureThermo(VolField<ThermoVar>& mixtureThermoField, const VolField<ComponentThermoVar>& vaporThermoField, const VolField<ComponentThermoVar>& liquidThermoField, const VolField<CompressibleMixture>& w) const
{
    #pragma omp parallel for
    for (size_t i = 0; i < vaporThermoField.size(); i++)
    {
        double rhoM = w[i].density();
        double alpha = w[i].alpha();

        double rhoL = liquidThermoField[i].density();
        double aL = liquidThermoField[i].soundSpeed();
        
        double rhoV = vaporThermoField[i].density();
        double aV = vaporThermoField[i].soundSpeed();

        mixtureThermoField[i][ThermoVar::P] = vaporThermoField[i].pressure();
        mixtureThermoField[i][ThermoVar::A] = (1.0/rhoM)/((alpha/(rhoL*aL*aL)) + (1.0 - rhoL*alpha)/(rhoV*rhoV*aV*aV));
    }

    #pragma omp parallel for
    for (size_t j = 0; j < vaporThermoField.boundarySize(); j++)
    {
        for (size_t i = 0; i < vaporThermoField.boundary(j).size(); i++)
        {
        double rhoM = w.boundary(j)[i].density();
        double alpha = w.boundary(j)[i].alpha();

        double rhoL = liquidThermoField.boundary(j)[i].density();
        double aL = liquidThermoField.boundary(j)[i].soundSpeed();
        
        double rhoV = vaporThermoField.boundary(j)[i].density();
        double aV = vaporThermoField.boundary(j)[i].soundSpeed();

        mixtureThermoField.boundary(j)[i][ThermoVar::P] = vaporThermoField.boundary(j)[i].pressure();
        mixtureThermoField.boundary(j)[i][ThermoVar::A] = (1.0/rhoM)/((alpha/(rhoL*aL*aL)) + (1.0 - rhoL*alpha)/(rhoV*rhoV*aV*aV));
        }        
    }
}

void Mixture::updateVaporThermoInternal(VolField<ThermoVar>& mixtureThermoField, const VolField<ComponentThermoVar>& vaporThermoField, const VolField<ComponentThermoVar>& liquidThermoField, const VolField<CompressibleMixture>& w) const
{
    #pragma omp parallel for
    for (size_t i = 0; i < vaporThermoField.size(); i++)
    {
        double rhoM = w[i].density();
        double alpha = w[i].alpha();

        double rhoL = liquidThermoField[i].density();
        double aL = liquidThermoField[i].soundSpeed();
        
        double rhoV = vaporThermoField[i].density();
        double aV = vaporThermoField[i].soundSpeed();

        mixtureThermoField[i][ThermoVar::P] = vaporThermoField[i].pressure();
        mixtureThermoField[i][ThermoVar::A] = (1.0/rhoM)/((alpha/(rhoL*aL*aL)) + (1.0 - rhoL*alpha)/(rhoV*rhoV*aV*aV));
    }
}

void Mixture::updateVaporThermo(VolField<ThermoVar>& mixtureThermoField, const VolField<ComponentThermoVar>& vaporThermoField, const VolField<ComponentThermoVar>& liquidThermoField, const VolField<CompressibleMixture>& w) const
{
    #pragma omp parallel for
    for (size_t i = 0; i < vaporThermoField.size(); i++)
    {
        double rhoM = w[i].density();
        double alpha = w[i].alpha();

        double rhoL = liquidThermoField[i].density();
        double aL = liquidThermoField[i].soundSpeed();
        
        double rhoV = vaporThermoField[i].density();
        double aV = vaporThermoField[i].soundSpeed();

        mixtureThermoField[i][ThermoVar::P] = vaporThermoField[i].pressure();
        mixtureThermoField[i][ThermoVar::A] = (1.0/rhoM)/((alpha/(rhoL*aL*aL)) + (1.0 - rhoL*alpha)/(rhoV*rhoV*aV*aV));
    }
}
