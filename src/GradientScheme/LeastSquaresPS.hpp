#ifndef LEASTSQUARESPS_HPP
#define LEASTSQUARESPS_HPP

#include "GradientScheme.hpp"

class LeastSquaresPS : public GradientScheme
{
    public:

        LeastSquaresPS() {}

        virtual ~LeastSquaresPS() {}

        virtual void init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList);

        Field<Mat<5,3>> calculateGradient(const Field<Compressible>& wl, const Field<Compressible>& wr, const Mesh& mesh) const;

        virtual Field<Mat<5,3>> calculateGradient(const Field<Compressible>& w, const std::vector<std::vector<Compressible>>& boundaryFields, const Mesh& mesh) const; //TODO template Field<Vars<5>>

    protected:
        Field<Mat<3,3>> MInv;

        std::vector<std::vector<int>> cellsStencil; //CellField
        std::vector<std::vector<std::pair<int, int>>> boundaryStencil; //CellField
        std::vector<std::vector<Vars<3>>> cellToCellDelta;
        std::vector<std::vector<Vars<3>>> cellToCellDeltaBoundary;

        void calculateInverseM(Field<Mat<3,3>> M); //TODO presun do Mat.hpp
        
        void createNodeStencil(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList);
        void calculateCellToCellDelta(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList);

};

#endif // LEASTSQUARESPS_HPP