#include <cmath>
#include "CompressibleMixture.hpp"


CompressibleMixture::CompressibleMixture(const Compressible& in) : Vars<9>()
{
    data[0] = in[0];
    data[1] = in[1];
    data[2] = in[2];
    data[3] = in[3];
    data[4] = in[4];
}


void CompressibleMixture::updateCompressiblePart(const Compressible& w)
{
    data[0] = w[0];
    data[1] = w[1];
    data[2] = w[2];
    data[3] = w[3];
    data[4] = w[4];
}

Compressible CompressibleMixture::getCompressible() const
{
    return Compressible({data[RHO], data[RHO_U], data[RHO_V], data[RHO_W], data[RHO_E]});
}        

double CompressibleMixture::alpha() const
{
    return (4.0/3.0)*M_PI*(data[RHO_Q3]);
}


double CompressibleMixture::density() const
{
    return data[RHO];
}

Vars<3> CompressibleMixture::velocity() const
{
    return Vars<3>({data[RHO_U] / data[RHO], data[RHO_V] / data[RHO], data[RHO_W] / data[RHO]});
}

double CompressibleMixture::absVelocity() const
{
    return sqrt((data[RHO_U]*data[RHO_U] + data[RHO_V]*data[RHO_V] + data[RHO_W]*data[RHO_W]) / (data[RHO]*data[RHO]));
}

double CompressibleMixture::absVelocity2() const
{
    return (data[RHO_U]*data[RHO_U] + data[RHO_V]*data[RHO_V] + data[RHO_W]*data[RHO_W]) / (data[RHO]*data[RHO]);
}

double CompressibleMixture::normalVelocity(const Vars<3>& normalVector) const
{
    return dot(this->velocity(), normalVector);
}

double CompressibleMixture::velocityU() const
{
    return data[RHO_U] / data[RHO];
}

double CompressibleMixture::velocityV() const
{
    return data[RHO_V] / data[RHO];
}

double CompressibleMixture::velocityW() const
{
    return data[RHO_W] / data[RHO];
}

double CompressibleMixture::totalEnergy() const
{
    return data[RHO_E] / data[RHO];
}

double CompressibleMixture::internalEnergy() const
{
    return data[RHO_E]/data[RHO] - 0.5*this->absVelocity2();
}

double CompressibleMixture::Q0() const
{
    return data[RHO_Q0]/data[RHO];
}

double CompressibleMixture::Q1() const
{
    return data[RHO_Q1]/data[RHO];
}

double CompressibleMixture::Q2() const
{
    return data[RHO_Q2]/data[RHO];
}

double CompressibleMixture::Q3() const
{
    return data[RHO_Q3]/data[RHO];
}

Vars<9> CompressibleMixture::flux(const Vars<3>& thermoData, const Vars<3>& normalVector) const
{
    Vars<3> velocity = this->velocity();

    double normalVelocity = dot(velocity, normalVector);

    return CompressibleMixture({data[RHO] * normalVelocity,
                         data[RHO] * velocity[0] * normalVelocity + thermoData[1] * normalVector[0],
                         data[RHO] * velocity[1] * normalVelocity + thermoData[1] * normalVector[1],
                         data[RHO] * velocity[2] * normalVelocity + thermoData[1] * normalVector[2],
                         (data[RHO_E] + thermoData[1]) * normalVelocity,
                         data[RHO_Q0] * normalVelocity,
                         data[RHO_Q1] * normalVelocity,
                         data[RHO_Q2] * normalVelocity,
                         data[RHO_Q3] * normalVelocity,});
}
