#include <cmath>

#include "GradientScheme.hpp"

void GradientScheme::init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{

}

Field<Mat<5,3>> GradientScheme::calculateGradient(const VolField<Compressible>& w, const Mesh& mesh) const
{
    Field<Mat<5,3>> grad(w.size());

    return grad;
}

Field<Mat<9,3>> GradientScheme::calculateGradient(const VolField<CompressibleMixture>& w, const Mesh& mesh) const
{
    Field<Mat<9,3>> grad(w.size());

    return grad;
}