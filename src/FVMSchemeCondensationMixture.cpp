#include <set>
#include <iostream>
#include <algorithm>
#include "FVMSchemeCondensationMixture.hpp"

#include "BoundaryCondition/Periodicity.hpp"

#include "outputCFD.hpp"

void FVMSchemeCondensationMixture::setReconstructionGradient(std::unique_ptr<GradientScheme> gradScheme_)
{
    gradientScheme = std::move(gradScheme_);
}
void FVMSchemeCondensationMixture::setReconstructionLimiter(std::unique_ptr<Limiter> limiter_)
{
    limiter = std::move(limiter_);
}

void FVMSchemeCondensationMixture::setCfl(double cfl_)
{
    cfl = cfl_;
}

void FVMSchemeCondensationMixture::setMaxIter(int maxIter_)
{
    maxIter = maxIter_;
}

void FVMSchemeCondensationMixture::setSaveEveryIter(int saveEveryIter_)
{
    saveEveryIter = saveEveryIter_;
}

void FVMSchemeCondensationMixture::setTargetError(double targetError_)
{
    targetError = targetError_;
}

void FVMSchemeCondensationMixture::setLocalTimeStep(bool localTimeStep_)
{
    localTimeStep = localTimeStep_;
}

void FVMSchemeCondensationMixture::setReconstructionSettings(bool reconstruction_)
{
    reconstruction = reconstruction_;
}

double FVMSchemeCondensationMixture::getCfl() const
{
    return cfl;
}

int FVMSchemeCondensationMixture::getMaxIter() const
{
    return maxIter;
}

double FVMSchemeCondensationMixture::getTargetError() const
{
    return targetError;
}

bool FVMSchemeCondensationMixture::getTimeStepsettings() const
{
    return localTimeStep;
}

bool FVMSchemeCondensationMixture::getReconstructionSettings() const
{
    return reconstruction;
}

const Mesh& FVMSchemeCondensationMixture::getMesh() const
{
    return mesh;
}

const Thermo* FVMSchemeCondensationMixture::getThermoRef()
{
    return thermo.get();
}


void FVMSchemeCondensationMixture::setInitialConditions(Compressible initialCondition)
{
    //w = Field<Compressible>(mesh.getCellsSize());
    for (int i = 0; i < mesh.getCellsSize(); i++)
    {
        w[i] = initialCondition;
    }

    //thermo->updateThermoInternal(thermoField, w);
}

void FVMSchemeCondensationMixture::setInitialConditionsPrimitive(Vars<5> initialCondition)
{
    Compressible CompressibleIC = thermo->primitiveToConservative(initialCondition);

    //w = Field<Compressible>(mesh.getCellsSize());
    //w = VolField<Compressible>(mesh);
    for (int i = 0; i < mesh.getCellsSize(); i++)
    {
        w[i] = CompressibleIC;
    }  
}

void FVMSchemeCondensationMixture::init()
{
    applyFreeBoundaryCondition(); // presenout do CaseSetter, pote mouno zredukovat tuto init() funkci

    gradientScheme->init(mesh, boundaryConditionList);

    w.getBoundaryData() = std::vector<std::vector<CompressibleMixture>>(boundaryConditionList.size());
    mixtureThermoField.getBoundaryData() = std::vector<std::vector<ThermoVar>>(boundaryConditionList.size());

    wl = Field<CompressibleMixture>(mesh.getFacesSize());
    wr = Field<CompressibleMixture>(mesh.getFacesSize());
    mixtureThermoFieldL = Field<ThermoVar>(mesh.getFacesSize());
    mixtureThermoFieldR = Field<ThermoVar>(mesh.getFacesSize());

    grad = Field<Mat<9,3>>(mesh.getCellsSize());
    phi = Field<Vars<9>>(mesh.getCellsSize());

    for (int boundaryId = 0; boundaryId < boundaryConditionList.size(); boundaryId++)
    {
        int size = boundaryConditionList[boundaryId]->getBoundary().facesIndex.size();
        w.boundary(boundaryId) = std::vector<CompressibleMixture>(size);
        mixtureThermoField.boundary(boundaryId) = std::vector<ThermoVar>(size);
    }    

    timeSteps = Field<double>(mesh.getCellsSize());

    fixGradient = 1000000;
    fixedGradStep = 10000000;
}

void FVMSchemeCondensationMixture::applyFreeBoundaryCondition()
{
    for(auto& boundaryCondition : boundaryConditionList)
    {
        if(boundaryCondition->getType() == BoundaryCondition::FREEBOUNDARY)
        {
            mesh.deleteBoundary(boundaryCondition->getBoundary().boundaryConditionName);
        }
    }

    boundaryConditionList.erase(std::remove_if(boundaryConditionList.begin(), boundaryConditionList.end(),
        [](const std::shared_ptr<BoundaryCondition> & boundary) { return (boundary->getType() == BoundaryCondition::FREEBOUNDARY); }),
        boundaryConditionList.end());

    for(auto& boundaryCondition : boundaryConditionList)
    {
        boundaryCondition->updateMeshBoundary(mesh);
    }

    for(auto& boundaryCondition : boundaryConditionList)
    {
        if(boundaryCondition->getType() == BoundaryCondition::PERIODICITY)
        {
            static_cast<Periodicity*>(boundaryCondition.get())->init(mesh);
        }
    }
}

void FVMSchemeCondensationMixture::setBoundaryConditions(std::vector<std::shared_ptr<BoundaryCondition>> boundaryConditions)
{
    boundaryConditionList = std::move(boundaryConditions);
}

void FVMSchemeCondensationMixture::calculateWlWr()
{
    //Without reconstruction

    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    for (int i = 0; i < mesh.getFacesSize(); i++)
    {        
        wl[i] = w[ownerIndexList[i]];

        int neighbour = neighborIndexList[i];
        if(neighbour >= 0)
        {
            wr[i] = w[neighbour];
        }
    }    
}

void FVMSchemeCondensationMixture::interpolateToFaces()
{
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    reconstruction = false; //TODO

    if(reconstruction)
    {
        for (int boundaryId = 0; boundaryId < w.getBoundaryData().size(); boundaryId++)
        {
            const std::vector<int>& boundaryFacesIndexList = boundaryConditionList[boundaryId]->getBoundary().facesIndex;
            for (int i = 0; i < boundaryFacesIndexList.size(); i++)
            {
                wr[boundaryFacesIndexList[i]] = w.boundary(boundaryId)[i];
            }
        }

        if (iter < fixGradient || iter % fixedGradStep == 0)
        {
            grad = gradientScheme->calculateGradient(w, mesh);
            phi = limiter->calculateLimiter(w, grad, mesh);
        }

        const std::vector<Cell>& cells = mesh.getCellList();
        const std::vector<Face>& faces = mesh.getFaceList();
        const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
        const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

        //EOS INTERVAL PRESERVING - todo mozna smazu
        //for (int i = 0; i < cells.size(); i++)
        //{
        //    if (w[i].density() < 0.11 || w[i].internalEnergy() < 2200000.0)
        //    {
        //        phi[i] = Vars<5>(0.0);
        //        //std::cout << "Err, EOS range; rho: " << w[i].density() << " e: " << w[i].internalEnergy() << std::endl;
        //    }
        //}

        for (int i = 0; i < faces.size(); i++)
        {
            int neighbour = neighborIndexList[i];
            if(neighbour >= 0)
            {
                Vars<9> wlDiff = dot(grad[ownerIndexList[i]], faces[i].midpoint - cells[ownerIndexList[i]].center);
                Vars<9> wrDiff = dot(grad[neighborIndexList[i]], faces[i].midpoint - cells[neighborIndexList[i]].center);

                wl[i] = w[ownerIndexList[i]] + phi[ownerIndexList[i]]*wlDiff;
                wr[i] = w[neighborIndexList[i]] + phi[neighborIndexList[i]]*wrDiff;
            }
            else
            {
                wl[i] = w[ownerIndexList[i]];
            }
        }

        for (auto & boundaryCondition : boundaryConditionList)
        {
            boundaryCondition->correct(w, wl, wr, grad, phi, mesh, thermo.get());
        }

        // TODO thermo->updateThermo(thermoFieldL, wl);
        // TODO thermo->updateThermo(thermoFieldR, wr);

    }
    else
    {
        for (int i = 0; i < mesh.getFacesSize(); i++)
        {        
            wl[i] = w[ownerIndexList[i]];
            mixtureThermoFieldL[i] = mixtureThermoField[ownerIndexList[i]];

            int neighbour = neighborIndexList[i];
            if(neighbour >= 0)
            {
                wr[i] = w[neighbour];
                mixtureThermoFieldR[i] = mixtureThermoField[neighbour];
            }
        }

        for (int boundaryId = 0; boundaryId < w.getBoundaryData().size(); boundaryId++)
        {
            const std::vector<int>& boundaryFacesIndexList = boundaryConditionList[boundaryId]->getBoundary().facesIndex;

            for (int i = 0; i < boundaryFacesIndexList.size(); i++)
            {
                wr[boundaryFacesIndexList[i]] = w.boundary(boundaryId)[i];
                mixtureThermoFieldR[boundaryFacesIndexList[i]] = mixtureThermoField.boundary(boundaryId)[i];
            }
        }
    }
}

void FVMSchemeCondensationMixture::calcBoundaryConditionFields()
{
    for (int boundaryConditionId = 0; boundaryConditionId < boundaryConditionList.size(); boundaryConditionId++)
    {
        w.boundary(boundaryConditionId) = boundaryConditionList[boundaryConditionId]->calc(w, mixtureThermoField, mesh, thermo.get());
    }
}

void FVMSchemeCondensationMixture::boundField()
{
    constexpr double minDensity = 0.1;
    constexpr double maxMach = 2.2;
    constexpr double minInternalEnergy = 2100000.0;
    constexpr double maxInternalEnergy = 2800000.0;
    

    const std::vector<Cell>& cells = mesh.getCellList();
    for (int i = 0; i < cells.size(); i++)
    {
        if (w[i].density() < minDensity)
        {
            //std::cout << "rho: " << w[i].density() << " U: " << w[i].absVelocity() << " energy: " << w[i].internalEnergy() << std::endl;
            w[i][Compressible::RHO_U] = w[i].velocityU()*minDensity;
            w[i][Compressible::RHO_V] = w[i].velocityV()*minDensity;
            w[i][Compressible::RHO_W] = w[i].velocityW()*minDensity;
            w[i][Compressible::RHO_E] = w[i].totalEnergy()*minDensity;
            w[i][Compressible::RHO] = minDensity;
        }

        /*if (w[i].machNumber() > maxMach)
        {
            double velocityCoeff = (maxMach*w[i].soundSpeed())/(w[i].absVelocity());
            w[i][Compressible::RHO_U] = w[i][Compressible::RHO_U]*velocityCoeff;
            w[i][Compressible::RHO_V] = w[i][Compressible::RHO_V]*velocityCoeff;
            w[i][Compressible::RHO_W] = w[i][Compressible::RHO_W]*velocityCoeff;
            //std::cout << "M i: " << std::endl;
        }*/

        if (w[i].internalEnergy() < minInternalEnergy)
        {
            w[i][Compressible::RHO_E] = minInternalEnergy*w[i].density();
            //std::cout << "e min i: " << i << std::endl;
        }

        /*if (w[i].internalEnergy() > maxInternalEnergy)
        {
            w[i][Compressible::RHO_E] = maxInternalEnergy*w[i].density();
            //std::cout << "e min i: " << i << std::endl;
        }*/
    }
}

void FVMSchemeCondensationMixture::updateTimeStep()
{
    const std::vector<Cell>& cells = mesh.getCellList();

    for (int i = 0; i < w.size(); i++)
    {
        Vars<3> projectedArea = cells[i].projectedArea;
        timeSteps[i] = cfl*(cells[i].volume/sum(projectedArea*(abs(w[i].velocity()) + Vars<3>(mixtureThermoField[i].soundSpeed()))));
        if (timeSteps[i] < 0.0)
        {
            std::cout << "cell index: " << i << std::endl;
        }
    }

    if (!localTimeStep)
    {
        double timeStep = min(timeSteps);
        timeSteps = Field<double>(timeSteps.size(), timeStep);

        time += timeStep;
    }
}

void FVMSchemeCondensationMixture::calculateFluxes()
{
    fluxes = fluxSolver->calculateFluxes(wl, wr, mixtureThermoFieldL, mixtureThermoFieldR, mesh.getFaceList());
}

Field<Vars<9>> FVMSchemeCondensationMixture::calculateResidual()
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();
    const std::vector<int>& neighbors = mesh.getNeighborIndexList();

    Field<Vars<9>> res(w.size());

    for (int i = 0; i < faces.size(); i++)
    {
        int owner = owners[i];
        int neighbor = neighbors[i];
        
        res[owner] -= fluxes[i];
        if (neighbor >= 0)
        {
            res[neighbor] += fluxes[i];
        }
    }

    for (int i = 0; i < cells.size(); i++)
    {
        res[i] = res[i]/cells[i].volume;
    }
    
    return res;
}

Field<CompressibleMixture> FVMSchemeCondensationMixture::getResults() const
{
    return w;
}

Field<ThermoVar> FVMSchemeCondensationMixture::getResultsThermo() const
{
    return mixtureThermoField;
}
