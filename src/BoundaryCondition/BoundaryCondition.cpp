#include "BoundaryCondition.hpp"


BoundaryCondition::BoundaryConditionType BoundaryCondition::getType() const
{
    return type;
}

Boundary BoundaryCondition::getBoundary() const
{
    return boundary;
}

std::vector<Compressible> BoundaryCondition::calc(const VolField<Compressible>& w, const VolField<ThermoVar>& thermoField, const Mesh& mesh, const Thermo * const thermoModel) const
{
    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();

    std::vector<Compressible> out(boundary.facesIndex.size());

    for (int i = 0; i < boundary.facesIndex.size(); i++)
    {
       out[i] = calculateState(w[ownerIndexList[boundary.facesIndex[i]]], thermoField[ownerIndexList[boundary.facesIndex[i]]], faceList[boundary.facesIndex[i]], thermoModel);
    }
    
    return out;
}

void BoundaryCondition::correct(const Field<Compressible>& w, Field<Compressible>& wl, Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{

}

void BoundaryCondition::updateMeshBoundary(const Mesh& mesh)
{
    const std::vector<Boundary>& boundaryList = mesh.getBoundaryList();
    for (int i = 0; i < boundaryList.size(); i++)
    {
        if(boundary.boundaryConditionName == boundaryList[i].boundaryConditionName)
        {
            boundary = boundaryList[i];
            break;
        }
    }
}


/////////////+

std::vector<CompressibleMixture> BoundaryCondition::calc(const VolField<CompressibleMixture>& w, const VolField<ThermoVar>& thermoField, const Mesh& mesh, const Thermo * const thermoModel) const
{
    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();

    std::vector<CompressibleMixture> out(boundary.facesIndex.size());

    for (int i = 0; i < boundary.facesIndex.size(); i++)
    {
       out[i] = calculateState(w[ownerIndexList[boundary.facesIndex[i]]], thermoField[ownerIndexList[boundary.facesIndex[i]]], faceList[boundary.facesIndex[i]], thermoModel);
    }
    
    return out;
}

void BoundaryCondition::correct(const Field<CompressibleMixture>& w, Field<CompressibleMixture>& wl, Field<CompressibleMixture>& wr, const Field<Mat<9,3>>& grad, const Field<Vars<9>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{

}

void BoundaryCondition::updateMeshBoundary(const Mesh& mesh)
{
    const std::vector<Boundary>& boundaryList = mesh.getBoundaryList();
    for (int i = 0; i < boundaryList.size(); i++)
    {
        if(boundary.boundaryConditionName == boundaryList[i].boundaryConditionName)
        {
            boundary = boundaryList[i];
            break;
        }
    }
}