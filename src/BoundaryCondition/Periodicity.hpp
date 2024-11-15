#ifndef PERIODICITY_HPP
#define PERIODICITY_HPP

#include "BoundaryCondition.hpp"

class Periodicity : public BoundaryCondition
{
    public:

        Periodicity(Boundary meshBoundary, Vars<3> faceMidpointShift_, std::string associatedBoundaryName_, const Mesh& mesh);

        void init(const Mesh& mesh);

        std::vector<int> getPeriodicityFacesIndex() const;
        std::vector<int> getPeriodicityFacesOwnersIndexes() const;
        Vars<3> getFaceShift() const;

        Compressible calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const;

        std::vector<Compressible> calc(const VolField<Compressible>& w, const VolField<ThermoVar>& thermoField, const Mesh& mesh, const Thermo * const thermoModel) const;

        void correct(const Field<Compressible>& w, Field<Compressible>& wl, Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const;

    private:
        std::vector<int> periodicityFacesIndex;
        std::vector<int> periodicityFacesOwnersIndexes;
        std::string associatedBoundaryName; //mozna nepotrebuju - potom zjistit, pripadne odstranit
        Vars<3> faceMidpointShift;
};

#endif // PERIODICITY_HPP