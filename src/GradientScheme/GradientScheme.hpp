#ifndef GRADIENTSCHEME_HPP
#define GRADIENTSCHEME_HPP

#include "../Mesh/Mesh.hpp"
#include "../BoundaryCondition/BoundaryCondition.hpp"
#include "../BoundaryCondition/Periodicity.hpp"
#include "../Field.hpp"
#include "../VolField.hpp"
#include "../Compressible.hpp"
#include "../Condensation/CompressibleMixture.hpp"

#include "../Mat.hpp"

class GradientScheme
{
    public:

        GradientScheme() {}

        virtual ~GradientScheme() {}

        virtual void init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList);

        virtual Field<Mat<5,3>> calculateGradient(const VolField<Compressible>& w, const Mesh& mesh) const;
        virtual Field<Mat<9,3>> calculateGradient(const VolField<CompressibleMixture>& w, const Mesh& mesh) const;
  
    protected:


};

#endif // GRADIENTSCHEME_HPP