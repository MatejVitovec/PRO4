#include "LeastSquare.hpp"
#include <iostream>

void LeastSquare::init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{
    calculateCellToCellDelta(mesh, boundaryConditionList);

    Field<Mat<3,3>> M(mesh.getCellsSize());

    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    for (int i = 0; i < owners.size(); i++)
    {
        Mat<3,3> auxM = outerProd(cellToCellDelta[i], cellToCellDelta[i]);
        M[owners[i]] += auxM;
        if (neighbours[i] >= 0)
        {
            M[neighbours[i]] += auxM;
        }
    }

    calculateInverseM(M);
}

void LeastSquare::calculateCellToCellDelta(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    cellToCellDelta = Field<Vars<3>>(faces.size());

    for (int i = 0; i < faces.size(); i++)
    {
        if (neighbours[i] != -1)
        {
            cellToCellDelta[i] = cells[neighbours[i]].center - cells[owners[i]].center;
        }
    }

    for (auto & boundaryCondition : boundaryConditionList)
    {
        if (boundaryCondition->getType() == BoundaryCondition::PERIODICITY)
        {
            std::vector<int> boundaryFaces = boundaryCondition->getBoundary().facesIndex;
            std::vector<int> associatedBoundaryFaces = static_cast<Periodicity*>(boundaryCondition.get())->getPeriodicityFacesIndex();

            Vars<3> faceMidpointShift = static_cast<Periodicity*>(boundaryCondition.get())->getFaceShift();

            for (int i = 0; i < boundaryFaces.size(); i++)
            {
                cellToCellDelta[boundaryFaces[i]] = cells[owners[associatedBoundaryFaces[i]]].center - cells[owners[boundaryFaces[i]]].center - faceMidpointShift;
            }
        }
        else
        {
            std::vector<int> boundaryFaces = boundaryCondition->getBoundary().facesIndex;

            for (int i = 0; i < boundaryFaces.size(); i++)
            {
                //TODO - mozna to bude fungovat
                cellToCellDelta[boundaryFaces[i]] = 2*(faces[boundaryFaces[i]].midpoint - cells[owners[boundaryFaces[i]]].center);
            }
        }
    }
}



void LeastSquare::calculateInverseM(Field<Mat<3,3>> M)
{
    MInv = Field<Mat<3,3>>(M.size());

    for (int i = 0; i < MInv.size(); i++)
    {
        MInv[i] = inv(M[i]);
    }
}


Field<Mat<5,3>> LeastSquare::calculateGradient(const Field<Compressible>& wl, const Field<Compressible>& wr, const Mesh& mesh) const
{
    Field<Mat<5,3>> grad(wl.size()); //TODO mozna ma byt cellList.size() misto wl.size()
    
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    Field<Mat<3,5>> b = Field<Mat<3,5>>(MInv.size());

    for (int i = 0; i < owners.size(); i++)
    {
        int neighbour = neighbours[i];

        Mat<3,5> rhs = outerProd(cellToCellDelta[i], wr[i] - wl[i]);
        b[owners[i]] += rhs;
        if (neighbour >= 0)
        {
            b[neighbour] -= rhs;
        }
    }

    for (int i = 0; i < cells.size(); i++)
    {
        grad[i] = transpose(dot(MInv[i], b[i]));        
    }

    return grad;
}