#include <cmath>
#include <iostream>
#include <fstream>

#include "Helmholtz.hpp"

#include <cmath>

#include "IapwsSaturation.hpp"

double IapwsSaturation::lnPsigmaPerCritP(double T) const
{
    double theta = calcTheta(T);
    return calcTau(T)*(satPCoeffs[0]*theta
                     + satPCoeffs[1]*std::pow(theta, 1.5)
                     + satPCoeffs[2]*std::pow(theta, 3)
                     + satPCoeffs[3]*std::pow(theta, 3.5)
                     + satPCoeffs[4]*std::pow(theta, 4)
                     + satPCoeffs[5]*std::pow(theta, 7.5));
}

double IapwsSaturation::vaporPressure(double T) const
{    
    return critP*std::exp(lnPsigmaPerCritP(T));
}

double IapwsSaturation::vaporPressureDiffT(double T) const
{
    double theta = calcTheta(T);
    return (-vaporPressure(T)/T)*(lnPsigmaPerCritP(T) +     satPCoeffs[0]
                                                   + 1.5*satPCoeffs[1]*std::pow(theta, 0.5)
                                                   +   3*satPCoeffs[2]*std::pow(theta, 2)
                                                   + 3.5*satPCoeffs[3]*std::pow(theta, 2.5)
                                                   +   4*satPCoeffs[4]*std::pow(theta, 3)
                                                   + 7.5*satPCoeffs[5]*std::pow(theta, 6.5));
}

double IapwsSaturation::saturatedVaporDensity(double T) const
{
    double theta = calcTheta(T);
    return critRho*(1.0 + satLiqRhoCoeffs[0]*std::pow(theta, 1.0/3.0)
                        + satLiqRhoCoeffs[1]*std::pow(theta, 2.0/3.0)
                        + satLiqRhoCoeffs[2]*std::pow(theta, 5.0/3.0)
                        + satLiqRhoCoeffs[3]*std::pow(theta, 16.0/3.0)
                        + satLiqRhoCoeffs[4]*std::pow(theta, 43.0/3.0)
                        + satLiqRhoCoeffs[5]*std::pow(theta, 110.0/3.0));
}

double IapwsSaturation::saturatedLiquidDensity(double T) const
{
    double theta = calcTheta(T);
    return critRho*std::exp(satVapRhoCoeffs[0]*std::pow(theta, 2.0/6.0)
                          + satVapRhoCoeffs[1]*std::pow(theta, 4.0/6.0)
                          + satVapRhoCoeffs[2]*std::pow(theta, 8.0/6.0)
                          + satVapRhoCoeffs[3]*std::pow(theta, 18.0/6.0)
                          + satVapRhoCoeffs[4]*std::pow(theta, 37.0/6.0)
                          + satVapRhoCoeffs[5]*std::pow(theta, 71.0/6.0));
}

double IapwsSaturation::alpha(double T) const
{
    double theta = calcTheta2(T);
    return alpha0*(dAlpha + alphaCoeffs[0]*std::pow(theta, -19)
                          + alphaCoeffs[1]*theta
                          + alphaCoeffs[2]*std::pow(theta, 4.5)
                          + alphaCoeffs[3]*std::pow(theta, 5)
                          + alphaCoeffs[4]*std::pow(theta, 54.5));
}

double IapwsSaturation::saturatedVaporEnthalpy(double T) const
{
    return alpha(T)/alpha0 + std::pow(10, 3)*(T/saturatedVaporDensity(T))*vaporPressureDiffT(T);
}

double IapwsSaturation::saturatedLiquidEnthalpy(double T) const
{
    return alpha(T)/alpha0 + std::pow(10, 3)*(T/saturatedLiquidDensity(T))*vaporPressureDiffT(T);
}

double IapwsSaturation::equalibriumEnthalpy(double T, double y) const
{
    return (y*(saturatedVaporEnthalpy(T) - saturatedLiquidEnthalpy(T)) + saturatedLiquidEnthalpy(T))/1000.0;
}

double IapwsSaturation::saturatedVaporInternalEnergy(double T) const
{
    return alpha(T)/alpha0 + std::pow(10, 3)*(T*vaporPressureDiffT(T) - vaporPressure(T))/saturatedVaporDensity(T);
}

double IapwsSaturation::saturatedLiquidInternalEnergy(double T) const
{
    return alpha(T)/alpha0 + std::pow(10, 3)*(T*vaporPressureDiffT(T) - vaporPressure(T))/saturatedLiquidDensity(T);
}