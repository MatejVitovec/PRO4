#include "Vars.hpp"

Vars<3> cross(const Vars<3>& u, const Vars<3>& v)
{
    return Vars<3>({u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]});
}

Vars<3> angleAngleToUnit(double xy, double xz)
{
    return unit(Vars<3>({1.0, std::atan(xy*(M_PI/180.0)), std::atan(xz*(M_PI/180.0))}));
}