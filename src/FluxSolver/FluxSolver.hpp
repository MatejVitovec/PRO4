#ifndef FLUXSOLVER_HPP
#define FLUXSOLVER_HPP

#include "../Mesh/Mesh.hpp"
#include "../Field.hpp"
#include "../Compressible.hpp"
#include "../ThermoVar.hpp"

#include "../Condensation/CompressibleMixture.hpp"

class FluxSolver
{
    public:

        FluxSolver() {}

        virtual ~FluxSolver() {}

        Field<Vars<5>> calculateFluxes(const Field<Compressible>& wl, const Field<Compressible>& wr, const Field<ThermoVar>& thermoFieldL, const Field<ThermoVar>& thermoFieldR, const std::vector<Face>& faceList) const;
        virtual Vars<5> claculateFlux(const Compressible& wl, const Compressible& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const = 0;

        Field<Vars<9>> calculateFluxes(const Field<CompressibleMixture>& wl, const Field<CompressibleMixture>& wr, const Field<ThermoVar>& thermoFieldL, const Field<ThermoVar>& thermoFieldR, const std::vector<Face>& faceList) const;
        virtual Vars<9> claculateFlux(const CompressibleMixture& wl, const CompressibleMixture& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const = 0;

    protected:

};

#endif // FLUXSOLVER_HPP