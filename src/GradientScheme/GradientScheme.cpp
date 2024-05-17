#include <cmath>

#include "GradientScheme.hpp"

void GradientScheme::init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{

}

Field<Mat<5,3>> GradientScheme::calculateGradient(const Field<Compressible>& wl, const Field<Compressible>& wr, const Mesh& mesh) const
{
    Field<Mat<5,3>> grad(mesh.getCellsSize());

    return grad;
}

Field<Mat<5,3>> GradientScheme::calculateGradient(const Field<Compressible>& w, const std::vector<std::vector<Compressible>>& boundaryFields, const Mesh& mesh) const
{
    Field<Mat<5,3>> grad(w.size());

    return grad;
}