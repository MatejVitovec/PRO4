#include <cmath>
#include "Compressible.hpp"

double Compressible::density() const
{
    return data[RHO];
}

Vars<3> Compressible::velocity() const
{
    return Vars<3>({data[RHO_U] / data[RHO], data[RHO_V] / data[RHO], data[RHO_W] / data[RHO]});
}

double Compressible::absVelocity() const
{
    return sqrt((data[RHO_U]*data[RHO_U] + data[RHO_V]*data[RHO_V] + data[RHO_W]*data[RHO_W]) / (data[RHO]*data[RHO]));
}

double Compressible::absVelocity2() const
{
    return (data[RHO_U]*data[RHO_U] + data[RHO_V]*data[RHO_V] + data[RHO_W]*data[RHO_W]) / (data[RHO]*data[RHO]);
}

double Compressible::normalVelocity(const Vars<3>& normalVector) const
{
    return dot(this->velocity(), normalVector);
}

double Compressible::velocityU() const
{
    return data[RHO_U] / data[RHO];
}

double Compressible::velocityV() const
{
    return data[RHO_V] / data[RHO];
}

double Compressible::velocityW() const
{
    return data[RHO_W] / data[RHO];
}

double Compressible::totalEnergy() const
{
    return data[RHO_E] / data[RHO];
}

double Compressible::internalEnergy() const
{
    return data[RHO_E]/data[RHO] - 0.5*this->absVelocity2();
}

Vars<5> Compressible::flux(const Vars<3>& thermoData, const Vars<3>& normalVector) const
{
    Vars<3> velocity = this->velocity();

    double normalVelocity = dot(velocity, normalVector);

    return Compressible({data[RHO] * normalVelocity,
                         data[RHO] * velocity[0] * normalVelocity + thermoData[1] * normalVector[0],
                         data[RHO] * velocity[1] * normalVelocity + thermoData[1] * normalVector[1],
                         data[RHO] * velocity[2] * normalVelocity + thermoData[1] * normalVector[2],
                         (data[RHO_E] + thermoData[1]) * normalVelocity});
}
