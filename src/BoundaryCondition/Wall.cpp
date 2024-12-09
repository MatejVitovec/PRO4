#include "Wall.hpp"

Compressible Wall::calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const
{
    Compressible out = w;

    Vars<3> normalVector = f.normalVector;
    Vars<3> ghostVelocity = w.velocity() - 2*w.normalVelocity(normalVector)*normalVector;

    out[Compressible::RHO_U] = w.density()*ghostVelocity[0];
    out[Compressible::RHO_V] = w.density()*ghostVelocity[1];
    out[Compressible::RHO_W] = w.density()*ghostVelocity[2];

    return out;
}

void Wall::correct(const Field<Compressible>& w, Field<Compressible>& wl, Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    for (auto & faceIndex : boundary.facesIndex)
    {
        Vars<5> wlDiff = dot(grad[ownerIndexList[faceIndex]], faces[faceIndex].midpoint - cells[ownerIndexList[faceIndex]].center);

        wl[faceIndex] = w[ownerIndexList[faceIndex]] + phi[ownerIndexList[faceIndex]]*wlDiff;
        wr[faceIndex] = calculateState(wl[faceIndex], ThermoVar(), faces[faceIndex], thermoModel);
    }
}


////////


CompressibleMixture Wall::calculateState(const CompressibleMixture& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const
{
    CompressibleMixture out = w;

    Vars<3> normalVector = f.normalVector;
    Vars<3> ghostVelocity = w.velocity() - 2*w.normalVelocity(normalVector)*normalVector;

    out[CompressibleMixture::RHO_U] = w.density()*ghostVelocity[0];
    out[CompressibleMixture::RHO_V] = w.density()*ghostVelocity[1];
    out[CompressibleMixture::RHO_W] = w.density()*ghostVelocity[2];

    return out;
}

void Wall::correct(const Field<CompressibleMixture>& w, Field<CompressibleMixture>& wl, Field<CompressibleMixture>& wr, const Field<Mat<9,3>>& grad, const Field<Vars<9>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    for (auto & faceIndex : boundary.facesIndex)
    {
        Vars<9> wlDiff = dot(grad[ownerIndexList[faceIndex]], faces[faceIndex].midpoint - cells[ownerIndexList[faceIndex]].center);

        wl[faceIndex] = w[ownerIndexList[faceIndex]] + phi[ownerIndexList[faceIndex]]*wlDiff;
        wr[faceIndex] = calculateState(wl[faceIndex], ThermoVar(), faces[faceIndex], thermoModel);
    }
}