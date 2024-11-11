#ifndef MEANPRESSUREOUTLET_HPP
#define MEANPRESSUREOUTLET_HPP

#include "BoundaryCondition.hpp"

class MeanPressureOutlet : public BoundaryCondition
{
    public:

        MeanPressureOutlet(Boundary meshBoundary, double pressure_) : BoundaryCondition(meshBoundary, MEANPRESSUREOUTLET), pressure(pressure_) {}

        Compressible calculateState(const Compressible& w, const Face& f, const Thermo * const thermoModel) const;
        std::vector<Compressible> calc(const Field<Compressible>& w, const Mesh& mesh, const Thermo * const thermoModel) const;
        

    private:
        double pressure;

        double calculateCorrectionConstant(const Mesh& mesh, const Field<Compressible>& w) const;
};



#endif // MEANPRESSUREOUTLET