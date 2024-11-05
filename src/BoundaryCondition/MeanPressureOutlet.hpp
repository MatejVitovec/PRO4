#ifndef MEANPRESSUREOUTLET_HPP
#define MEANPRESSUREOUTLET_HPP

#include "BoundaryCondition.hpp"

class MeanPressureOutlet : public BoundaryCondition
{
    public:

        MeanPressureOutlet(Boundary meshBoundary, double pressure_) : BoundaryCondition(meshBoundary, MEANPRESSUREOUTLET), pressure(pressure_) {}

        Compressible calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const;
        void apply(const std::vector<int>& ownerIndexList, const std::vector<Face>& faces, const Field<Compressible>& w, Field<Compressible>& wr, const Thermo * const thermoModel) const;
        

    private:
        double pressure;

        double calculateCorrectionConstant(const std::vector<int>& ownerIndexList, const std::vector<Face>& faces, const Field<Compressible>& w) const;
};

#endif // MEANPRESSUREOUTLET