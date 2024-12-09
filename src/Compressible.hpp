#ifndef COMPRESSIBLE_HPP
#define COMPRESSIBLE_HPP

#include <vector>
#include <memory>

#include "Vars.hpp"

class Compressible : public Vars<5>
{
    public:
        using Vars<5>::operator+=;
        using Vars<5>::operator-=;

        enum {RHO, RHO_U, RHO_V, RHO_W, RHO_E};

        Compressible() : Vars<5>() {}
        Compressible(const Vars<5>& varsIn) : Vars<5>(varsIn) {}
        Compressible(const std::array<double, 5>& in) : Vars<5>(in) {}
        //Compressible(const std::array<double, 5>& inRef ) : 

        virtual ~Compressible() {}

        double density() const;
        Vars<3> velocity() const;
        double absVelocity() const;
        double absVelocity2() const;
        double normalVelocity(const Vars<3>& normalVector) const;
        double velocityU() const;
        double velocityV() const;
        double velocityW() const;
        double totalEnergy() const;
        double internalEnergy() const;

        Vars<5> flux(const Vars<3>& thermoData, const Vars<3>& normalVector) const;

};

#endif // COMPRESSIBLE_HPP