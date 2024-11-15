#include <cmath>
#include "ThermoVar.hpp"

/*void ThermoVar::setThermoVar(Vars<3> thermoProp)
{
    thermoVar = thermoProp;
}

double ThermoVar::density() const
{
    return data[RHO];
}

Vars<3> ThermoVar::velocity() const
{
    return Vars<3>({data[RHO_U] / data[RHO], data[RHO_V] / data[RHO], data[RHO_W] / data[RHO]});
}

double ThermoVar::absVelocity() const
{
    return sqrt((data[RHO_U]*data[RHO_U] + data[RHO_V]*data[RHO_V] + data[RHO_W]*data[RHO_W]) / (data[RHO]*data[RHO]));
}

double ThermoVar::absVelocity2() const
{
    return (data[RHO_U]*data[RHO_U] + data[RHO_V]*data[RHO_V] + data[RHO_W]*data[RHO_W]) / (data[RHO]*data[RHO]);
}

double ThermoVar::normalVelocity(const Vars<3>& normalVector) const
{
    return dot(this->velocity(), normalVector);
}

double ThermoVar::velocityU() const
{
    return data[RHO_U] / data[RHO];
}

double ThermoVar::velocityV() const
{
    return data[RHO_V] / data[RHO];
}

double ThermoVar::velocityW() const
{
    return data[RHO_W] / data[RHO];
}

double ThermoVar::totalEnergy() const
{
    return data[RHO_E] / data[RHO];
}*/

double ThermoVar::temperature() const
{
    return data[T];
}

double ThermoVar::pressure() const
{
    return data[P];
}

double ThermoVar::soundSpeed() const
{
    return data[A];
}

/*double ThermoVar::internalEnergy() const
{
    return data[RHO_E]/data[RHO] - 0.5*this->absVelocity2();
}

double ThermoVar::machNumber() const
{
    return absVelocity()/soundSpeed();
}

Vars<3> ThermoVar::thermo() const
{
    return thermoVar;
}

Vars<5> ThermoVar::flux(const Vars<3>& normalVector) const
{
    Vars<3> velocity = this->velocity();

    double normalVelocity = dot(velocity, normalVector);

    return ThermoVar({data[RHO] * normalVelocity,
                         data[RHO] * velocity[0] * normalVelocity + thermoVar[P] * normalVector[0],
                         data[RHO] * velocity[1] * normalVelocity + thermoVar[P] * normalVector[1],
                         data[RHO] * velocity[2] * normalVelocity + thermoVar[P] * normalVector[2],
                         (data[RHO_E]+ thermoVar[P]) * normalVelocity});  //tady byla entalpie - muze byt blbe
}

Vars<5> ThermoVar::primitive() const
{
    return Vars<5>({data[RHO],
                    data[RHO_U] / data[RHO],
                    data[RHO_V] / data[RHO],
                    data[RHO_W] / data[RHO],
                    thermoVar[P]});
}

void ThermoVar::operator+=(const ThermoVar& v)
{
    data[0] += v[0];
    data[1] += v[1];
    data[2] += v[2];
    data[3] += v[3];
    data[4] += v[4];
}

void ThermoVar::operator-=(const ThermoVar& v)
{
    data[0] -= v[0];
    data[1] -= v[1];
    data[2] -= v[2];
    data[3] -= v[3];
    data[4] -= v[4];
}

void ThermoVar::operator+=(const Vars<5>& v)
{
    data[0] += v[0];
    data[1] += v[1];
    data[2] += v[2];
    data[3] += v[3];
    data[4] += v[4];
}

void ThermoVar::operator-=(const Vars<5>& v)
{
    data[0] -= v[0];
    data[1] -= v[1];
    data[2] -= v[2];
    data[3] -= v[3];
    data[4] -= v[4];
}


//////////////Non member operators///////////////////

// u == v
bool operator== (const ThermoVar& u, const ThermoVar& v)
{
    if(u[0] == v[0] && u[1] == v[1] && u[2] == v[2] && u[3] == v[3] && u[4] == v[4]) return true;
    return false;
}

// u + v
ThermoVar operator+ (const ThermoVar& u, const ThermoVar& v)
{
    ThermoVar out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] += v[i];
    }
    
    return out;
}

// u - v
ThermoVar operator- (const ThermoVar& u, const ThermoVar& v)
{
    ThermoVar out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] -= v[i];
    }
    
    return out;
}

// w * u
ThermoVar operator* (const ThermoVar& u, const ThermoVar& v)
{
    ThermoVar out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] *= v[i];
    }
    
    return out;
}

ThermoVar operator/ (const ThermoVar& u, const ThermoVar& v)
{
    ThermoVar out = u;
    for (int i = 0; i < 5; i++)
    {
        if (out[i] != 0)
        {
            out[i] /= v[i];
        }
    }
    
    return out;
}

// a * u
ThermoVar operator* (double a, const ThermoVar& u)
{
    ThermoVar out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] *= a;
    }
    
    return out;
}

// u * a
ThermoVar operator* (const ThermoVar& u, double a)
{
    ThermoVar out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] *= a;
    }
    
    return out;
}

// u / a
ThermoVar operator/ (const ThermoVar& u, double a)
{
    ThermoVar out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] /= a;
    }
    
    return out;
}

// ThermoVar, Vars<5>

// u + v
ThermoVar operator+ (const ThermoVar& u, const Vars<5>& v)
{
    ThermoVar out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] += v[i];
    }
    
    return out;
}

// u - v
ThermoVar operator- (const ThermoVar& u, const Vars<5>& v)
{
    ThermoVar out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] -= v[i];
    }
    
    return out;
}

// w * u
ThermoVar operator* (const ThermoVar& u, const Vars<5>& v)
{
    ThermoVar out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] *= v[i];
    }
    
    return out;
}


//////////////Non member function///////////////////

ThermoVar sqrt(const ThermoVar& u)
{
    ThermoVar out = u;
    for (int i = 0; i < 5; i++)
    {
        out[i] = std::sqrt(u[i]);
    }
    
    return out;
}*/