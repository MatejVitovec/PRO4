#ifndef COMPRESSIBLEMIXTURE_HPP
#define COMPRESSIBLEMIXTURE_HPP

#include <vector>
#include <memory>

#include "../Compressible.hpp"

//TODO [] operator

class CompressibleMixture : public Vars<9>
{
    public:
        using Vars<9>::operator+=;
        using Vars<9>::operator-=;

        enum {RHO, RHO_U, RHO_V, RHO_W, RHO_E, RHO_Q0, RHO_Q1, RHO_Q2, RHO_Q3};

        CompressibleMixture() : Vars<9>() {}
        CompressibleMixture(const Vars<9>& varsIn) : Vars<9>(varsIn) {}
        CompressibleMixture(const std::array<double, 9>& in) : Vars<9>(in) {}
        CompressibleMixture(const Compressible& in);
        virtual ~CompressibleMixture() {}

        operator Compressible() const //Conversion operator
        {
            return Compressible({data[RHO], data[RHO_U], data[RHO_V], data[RHO_W], data[RHO_E]});
        }

        void updateCompressiblePart(const Compressible& w);
        Compressible getCompressible() const;

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

        double Q0() const;
        double Q1() const;
        double Q2() const;
        double Q3() const;

        double alpha() const;

        Vars<9> flux(const Vars<3>& thermoData, const Vars<3>& normalVector) const;

};

#endif // COMPRESSIBLEMIXTURE_HPP