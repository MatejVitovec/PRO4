#ifndef IAPWSSATURATION
#define IAPWSSATURATION

#include <array>
#include <vector>

#include "NonLinearSolver/NewtonMethod.hpp"

class IapwsSaturation
{
    public:

        IapwsSaturation();

        virtual ~IapwsSaturation() {}

        double vaporPressure(double T) const;
        double vaporPressureDiffT(double T) const;

        double saturatedVaporDensity(double T) const;
        double saturatedLiquidDensity(double T) const;

        double saturatedVaporEnthalpy(double T) const;
        double saturatedLiquidEnthalpy(double T) const;

        double saturatedVaporInternalEnergy(double T) const;
        double saturatedLiquidInternalEnergy(double T) const;

        double equalibriumEnthalpy(double T, double y) const;

    private:
        double critT = 647.096;
        double critRho = 322.0;
        double specGasConst = 461.51805;
        double critP = 22064000.0;

        double alpha0 = 1000.0;
        double dAlpha = -1135.905627715;
        double dPsi = 2319.5246;

        std::array<double, 6> satPCoeffs = {-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, -1.80122502};
        std::array<double, 6> satLiqRhoCoeffs = {1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -674694.450};
        std::array<double, 6> satVapRhoCoeffs = {-2.03150240, -2.68302940, -5.38626492, -17.2991605, -44.7586581, -63.9201063};
        std::array<double, 5> alphaCoeffs = {-0.0000000565134998, 2690.66631, 127.287297, -135.003439, 0.981825814};

        inline double calcTheta(double T) const { return 1.0 - T/critT; }
        inline double calcTheta2(double T) const { return T/critT; }
        inline double calcDelta(double rho) const { return rho/critRho; }
        inline double calcTau(double T) const { return critT/T; }
        inline double psi0() { return alpha0/critT; }

        double lnPsigmaPerCritP(double T) const;
        double alpha(double T) const;

};

#endif // IAPWSSATURATION