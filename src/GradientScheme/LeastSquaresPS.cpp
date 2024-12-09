#include "LeastSquaresPS.hpp"
#include <iostream>
#include <algorithm>



void LeastSquaresPS::createNodeStencil(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{
    const std::vector<Vars<3>>& nodeList = mesh.getNodeList();
    const std::vector<Cell>& cellList = mesh.getCellList();
    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& neighbors = mesh.getNeighborIndexList();
    const std::vector<Boundary>& boundarylist = mesh.getBoundaryList();

    cellsStencil = std::vector<std::vector<int>>(cellList.size()); //presun do konstruktoru (init)

    ////internal cells
    std::vector<std::vector<int>> nodeCells = std::vector<std::vector<int>>(nodeList.size());
    
    for (int cellId = 0; cellId < cellList.size(); cellId++)
    {
        const std::vector<int>& cellNodes = cellList[cellId].nodesIndex; //TODO encapsulate in mesh classes
        for (int nodeId = 0; nodeId < cellNodes.size(); nodeId++)
        {
            nodeCells[cellNodes[nodeId]].push_back(cellId);
        }        
    }

    for (int i = 0; i < cellList.size(); i++)
    {
        std::vector<int> cellStencil;
        for (int j = 0; j < cellList[i].nodesIndex.size(); j++)
        {            
            const std::vector<int>& aux = nodeCells[cellList[i].nodesIndex[j]];
            cellStencil.insert(cellStencil.end(), aux.begin(), aux.end());
        }

        std::sort(cellStencil.begin(), cellStencil.end());
        cellStencil.erase(std::unique(cellStencil.begin(), cellStencil.end()), cellStencil.end());
        cellStencil.erase(std::remove(cellStencil.begin(), cellStencil.end(), i), cellStencil.end());
        cellStencil.shrink_to_fit();

        cellsStencil[i] = cellStencil;
    }


    ////boundary cells

    boundaryStencil = std::vector<std::vector<std::pair<int, int>>>(cellList.size());

    for (int boundaryId = 0; boundaryId < boundarylist.size(); boundaryId++)
    {
        for (int i = 0; i < boundarylist[boundaryId].facesIndex.size(); i++)
        {
            for (int j = 0; j < faceList[boundarylist[boundaryId].facesIndex[i]].nodesIndex.size(); j++)
            {
                for (int k = 0; k < nodeCells[faceList[boundarylist[boundaryId].facesIndex[i]].nodesIndex[j]].size(); k++)
                {
                    //cellBoundaryFaces[nodeCells[faceList[boundary.facesIndex[i]].nodesIndex[j]][k]].push_back(boundary.facesIndex[i]);
                    boundaryStencil[nodeCells[faceList[boundarylist[boundaryId].facesIndex[i]].nodesIndex[j]][k]].push_back(
                        std::pair<int, int>({boundaryId, /*boundarylist[boundaryId].facesIndex[i]*/i}));
                }
            }
        }
    }

    for (int i = 0; i < boundaryStencil.size(); i++)
    {
        std::sort(boundaryStencil[i].begin(), boundaryStencil[i].end());
        boundaryStencil[i].erase(std::unique(boundaryStencil[i].begin(), boundaryStencil[i].end()), boundaryStencil[i].end());
        boundaryStencil[i].shrink_to_fit();
    }

    boundaryStencil.shrink_to_fit();
}

void LeastSquaresPS::calculateCellToCellDelta(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{
    const std::vector<Vars<3>>& nodeList = mesh.getNodeList();
    const std::vector<Cell>& cellList = mesh.getCellList();
    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();
    const std::vector<int>& neighbors = mesh.getNeighborIndexList();
    const std::vector<Boundary>& boundaryList = mesh.getBoundaryList();

    cellToCellDelta = std::vector<std::vector<Vars<3>>>(cellsStencil.size());
    cellToCellDeltaBoundary = std::vector<std::vector<Vars<3>>>(cellsStencil.size());

    for (int cellId = 0; cellId < cellsStencil.size(); cellId++)
    {
        for (int i = 0; i < cellsStencil[cellId].size(); i++)
        {
            cellToCellDelta[cellId].push_back(cellList[cellsStencil[cellId][i]].center - cellList[cellId].center);
        }

        for (int i = 0; i < boundaryStencil[cellId].size(); i++)
        {
            if (boundaryConditionList[boundaryStencil[cellId][i].first]->getType() == BoundaryCondition::PERIODICITY)
            {
                const std::vector<int>& associatedBoundaryFaces = static_cast<Periodicity*>(boundaryConditionList[boundaryStencil[cellId][i].first].get())->getPeriodicityFacesIndex();
                const std::vector<int>& associatedBoundaryCells = static_cast<Periodicity*>(boundaryConditionList[boundaryStencil[cellId][i].first].get())->getPeriodicityFacesOwnersIndexes();
                Vars<3> faceMidpointShift = static_cast<Periodicity*>(boundaryConditionList[boundaryStencil[cellId][i].first].get())->getFaceShift();

                cellToCellDeltaBoundary[cellId].push_back(cellList[owners[associatedBoundaryFaces[boundaryStencil[cellId][i].second]]].center - cellList[cellId].center - faceMidpointShift);
            }
            else
            {
                cellToCellDeltaBoundary[cellId].push_back(2*(faceList[boundaryList[boundaryStencil[cellId][i].first].facesIndex[boundaryStencil[cellId][i].second]].midpoint - cellList[cellId].center));
            }
        }
    }

    for (int i = 0; i < cellToCellDelta.size(); i++)
    {
        for (int j = 0; j < cellToCellDelta[i].size(); j++)
        {
            cellToCellDelta[i][j] = zeroSmallNumbers(cellToCellDelta[i][j]);
        }
    }
    for (int i = 0; i < cellToCellDeltaBoundary.size(); i++)
    {
        for (int j = 0; j < cellToCellDeltaBoundary[i].size(); j++)
        {
            cellToCellDeltaBoundary[i][j] = zeroSmallNumbers(cellToCellDeltaBoundary[i][j]);
        }
    }
}

void LeastSquaresPS::init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{
    createNodeStencil(mesh, boundaryConditionList);
    calculateCellToCellDelta(mesh, boundaryConditionList);

    Field<Mat<3,3>> M(mesh.getCellsSize());

    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    for (int cellId = 0; cellId < cellsStencil.size(); cellId++)
    {
        for (int i = 0; i < cellsStencil[cellId].size(); i++)
        {
            M[cellId] += outerProd(cellToCellDelta[cellId][i], cellToCellDelta[cellId][i]);
        }
    }

    calculateInverseM(M);
}


void LeastSquaresPS::calculateInverseM(Field<Mat<3,3>> M)
{
    MInv = Field<Mat<3,3>>(M.size());

    for (int i = 0; i < MInv.size(); i++)
    {
        MInv[i] = invSingularCheck(M[i]);
    }
}

Field<Mat<5,3>> LeastSquaresPS::calculateGradient(const VolField<Compressible>& w, const Mesh& mesh) const
{
    Field<Mat<5,3>> grad(w.size());
    
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    Field<Mat<3,5>> b = Field<Mat<3,5>>(MInv.size());

    for (int cellId = 0; cellId < cellsStencil.size(); cellId++)
    {
        for (int i = 0; i < cellsStencil[cellId].size(); i++)
        {
            b[cellId] += outerProd(cellToCellDelta[cellId][i], w[cellsStencil[cellId][i]] - w[cellId]);
        }
    }

    for (int i = 0; i < cells.size(); i++)
    {
        grad[i] = transpose(dot(MInv[i], b[i]));        
    }

    return grad;
}

Field<Mat<9,3>> LeastSquaresPS::calculateGradient(const VolField<CompressibleMixture>& w, const Mesh& mesh) const
{
    Field<Mat<9,3>> grad(w.size());
    
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& neighbours = mesh.getNeighborIndexList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();

    Field<Mat<3,9>> b = Field<Mat<3,9>>(MInv.size());

    for (int cellId = 0; cellId < cellsStencil.size(); cellId++)
    {
        for (int i = 0; i < cellsStencil[cellId].size(); i++)
        {
            b[cellId] += outerProd(cellToCellDelta[cellId][i], w[cellsStencil[cellId][i]] - w[cellId]);
        }
    }

    for (int i = 0; i < cells.size(); i++)
    {
        grad[i] = transpose(dot(MInv[i], b[i]));        
    }

    return grad;
}