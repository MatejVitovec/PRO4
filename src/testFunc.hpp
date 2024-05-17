#include <iostream>
#include <algorithm>

#include "Mesh/Mesh.hpp"
#include "BoundaryCondition/BoundaryCondition.hpp"
#include "BoundaryCondition/Periodicity.hpp"
#include "Field.hpp"


std::vector<std::vector<int>> createNodeStencil(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList)
{
    const std::vector<Vars<3>>& nodeList = mesh.getNodeList();
    const std::vector<Cell>& cellList = mesh.getCellList();
    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& neighbors = mesh.getNeighborIndexList();
    const std::vector<Boundary>& boundaries = mesh.getBoundaryList();

    std::vector<std::vector<int>> cellsStencil(cellList.size());

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
    std::cout << "start" <<std::endl;

    std::vector<std::vector<std::pair<int, int>>> boundaryStencil(cellList.size()); //presun do class

    for (int boundaryId = 0; boundaryId < boundaries.size(); boundaryId++)
    {
        //std::vector<std::vector<int>> cellBoundaryFaces = std::vector<std::vector<int>>();

        for (int i = 0; i < boundaries[boundaryId].facesIndex.size(); i++)
        {
            for (int j = 0; j < faceList[boundaries[boundaryId].facesIndex[i]].nodesIndex.size(); j++)
            {
                for (int k = 0; k < nodeCells[faceList[boundaries[boundaryId].facesIndex[i]].nodesIndex[j]].size(); k++)
                {
                    //cellBoundaryFaces[nodeCells[faceList[boundary.facesIndex[i]].nodesIndex[j]][k]].push_back(boundary.facesIndex[i]);
                    boundaryStencil[nodeCells[faceList[boundaries[boundaryId].facesIndex[i]].nodesIndex[j]][k]].push_back(std::pair<int, int>({boundaryId, boundaries[boundaryId].facesIndex[i]}));
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
    

    //return boundaryStencil;
    return cellsStencil;
}





