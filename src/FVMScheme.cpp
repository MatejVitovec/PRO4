#include <set>
#include <iostream>
#include <algorithm>
#include "FVMScheme.hpp"

#include "BoundaryCondition/Periodicity.hpp"

#include "outputCFD.hpp"

void FVMScheme::setReconstructionGradient(std::unique_ptr<GradientScheme> gradScheme_)
{
    gradientScheme = std::move(gradScheme_);
}
void FVMScheme::setReconstructionLimiter(std::unique_ptr<Limiter> limiter_)
{
    limiter = std::move(limiter_);
}

void FVMScheme::setCfl(double cfl_)
{
    cfl = cfl_;
}

void FVMScheme::setMaxIter(int maxIter_)
{
    maxIter = maxIter_;
}

void FVMScheme::setSaveEveryIter(int saveEveryIter_)
{
    saveEveryIter = saveEveryIter_;
}

void FVMScheme::setTargetError(double targetError_)
{
    targetError = targetError_;
}

void FVMScheme::setLocalTimeStep(bool localTimeStep_)
{
    localTimeStep = localTimeStep_;
}

void FVMScheme::setReconstructionSettings(bool reconstruction_)
{
    reconstruction = reconstruction_;
}

double FVMScheme::getCfl() const
{
    return cfl;
}

int FVMScheme::getMaxIter() const
{
    return maxIter;
}

double FVMScheme::getTargetError() const
{
    return targetError;
}

bool FVMScheme::getTimeStepsettings() const
{
    return localTimeStep;
}

bool FVMScheme::getReconstructionSettings() const
{
    return reconstruction;
}

const Mesh& FVMScheme::getMesh() const
{
    return mesh;
}

const Thermo* FVMScheme::getThermoRef()
{
    return thermo.get();
}


void FVMScheme::setInitialConditions(Compressible initialCondition)
{
    //w = Field<Compressible>(mesh.getCellsSize());
    for (int i = 0; i < mesh.getCellsSize(); i++)
    {
        w[i] = initialCondition;
    }

    thermo->updateThermoInternal(thermoField, w);
}

void FVMScheme::setInitialConditionsPrimitive(Vars<5> initialCondition)
{
    Compressible CompressibleIC = thermo->primitiveToConservative(initialCondition);

    //w = Field<Compressible>(mesh.getCellsSize());
    //w = VolField<Compressible>(mesh);
    for (int i = 0; i < mesh.getCellsSize(); i++)
    {
        w[i] = CompressibleIC;
    }  
}

void FVMScheme::init()
{
    applyFreeBoundaryCondition(); // presenout do CaseSetter, pote mouno zredukovat tuto init() funkci

    gradientScheme->init(mesh, boundaryConditionList);

    w.getBoundaryData() = std::vector<std::vector<Compressible>>(boundaryConditionList.size());
    thermoField.getBoundaryData() = std::vector<std::vector<ThermoVar>>(boundaryConditionList.size());

    wl = Field<Compressible>(mesh.getFacesSize());
    wr = Field<Compressible>(mesh.getFacesSize());
    thermoFieldL = Field<ThermoVar>(mesh.getFacesSize());
    thermoFieldR = Field<ThermoVar>(mesh.getFacesSize());
    //boundaryFields = std::vector<std::vector<Compressible>>(boundaryConditionList.size());

    grad = Field<Mat<5,3>>(mesh.getCellsSize());
    phi = Field<Vars<5>>(mesh.getCellsSize());

    for (int boundaryId = 0; boundaryId < boundaryConditionList.size(); boundaryId++)
    {
        int size = boundaryConditionList[boundaryId]->getBoundary().facesIndex.size();
        //boundaryFields[boundaryId] = std::vector<Compressible>(size);
        w.boundary(boundaryId) = std::vector<Compressible>(size);
        thermoField.boundary(boundaryId) = std::vector<ThermoVar>(size);
    }    

    timeSteps = Field<double>(mesh.getCellsSize());

    fixGradient = 1000000;
    fixedGradStep = 10000000;
}

void FVMScheme::applyFreeBoundaryCondition()
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

void FVMScheme::setBoundaryConditions(std::vector<std::shared_ptr<BoundaryCondition>> boundaryConditions)
{
    boundaryConditionList = std::move(boundaryConditions);
}

void FVMScheme::calculateWlWr()
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

void FVMScheme::interpolateToFaces()
{
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    if(reconstruction)
    {
        /*for (int boundaryId = 0; boundaryId < boundaryFields.size(); boundaryId++)
        {
            const std::vector<int>& boundaryFacesIndexList = boundaryConditionList[boundaryId]->getBoundary().facesIndex;
            for (int i = 0; i < boundaryFacesIndexList.size(); i++)
            {
                wr[boundaryFacesIndexList[i]] = boundaryFields[boundaryId][i];
            }
        }*/
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
        /*for (int i = 0; i < cells.size(); i++)
        {
            if (w[i].density() < 0.11 || w[i].internalEnergy() < 2200000.0)
            {
                phi[i] = Vars<5>(0.0);
                //std::cout << "Err, EOS range; rho: " << w[i].density() << " e: " << w[i].internalEnergy() << std::endl;
            }
        }*/

        for (int i = 0; i < faces.size(); i++)
        {
            int neighbour = neighborIndexList[i];
            if(neighbour >= 0)
            {
                Vars<5> wlDiff = dot(grad[ownerIndexList[i]], faces[i].midpoint - cells[ownerIndexList[i]].center);
                Vars<5> wrDiff = dot(grad[neighborIndexList[i]], faces[i].midpoint - cells[neighborIndexList[i]].center);

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
        
        //wl = thermo->updateField(wl);
        //wr = thermo->updateField(wr);

        thermo->updateThermo(thermoFieldL, wl);
        thermo->updateThermo(thermoFieldR, wr);

    }
    else
    {
        for (int i = 0; i < mesh.getFacesSize(); i++)
        {        
            wl[i] = w[ownerIndexList[i]];
            thermoFieldL[i] = thermoField[ownerIndexList[i]];

            int neighbour = neighborIndexList[i];
            if(neighbour >= 0)
            {
                wr[i] = w[neighbour];
                thermoFieldR[i] = thermoField[neighbour];
            }
        }

        for (int boundaryId = 0; boundaryId < w.getBoundaryData().size(); boundaryId++)
        {
            const std::vector<int>& boundaryFacesIndexList = boundaryConditionList[boundaryId]->getBoundary().facesIndex;

            for (int i = 0; i < boundaryFacesIndexList.size(); i++)
            {
                wr[boundaryFacesIndexList[i]] = w.boundary(boundaryId)[i];
                thermoFieldR[boundaryFacesIndexList[i]] = thermoField.boundary(boundaryId)[i];
            }
        }
    }
}

void FVMScheme::calcBoundaryConditionFields()
{
    for (int boundaryConditionId = 0; boundaryConditionId < boundaryConditionList.size(); boundaryConditionId++)
    {
        //boundaryFields[boundaryConditionId] = boundaryConditionList[boundaryConditionId]->calc(w, mesh, thermo.get());
        w.boundary(boundaryConditionId) = boundaryConditionList[boundaryConditionId]->calc(w, thermoField, mesh, thermo.get());
    }
}

void FVMScheme::boundField()
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

void FVMScheme::updateTimeStep()
{
    const std::vector<Cell>& cells = mesh.getCellList();

    for (int i = 0; i < w.size(); i++)
    {
        Vars<3> projectedArea = cells[i].projectedArea;
        //timeSteps[i] = cfl*(cells[i].volume/sum(projectedArea*(abs(w[i].velocity()) + Vars<3>(w[i].soundSpeed()))));
        timeSteps[i] = cfl*(cells[i].volume/sum(projectedArea*(abs(w[i].velocity()) + Vars<3>(thermoField[i].soundSpeed()))));
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

void FVMScheme::calculateFluxes()
{
    fluxes = fluxSolver->calculateFluxes(wl, wr, thermoFieldL, thermoFieldR, mesh.getFaceList());
}

Field<Vars<5>> FVMScheme::calculateResidual()
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();
    const std::vector<int>& neighbors = mesh.getNeighborIndexList();

    Field<Vars<5>> res(w.size());

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

Field<Compressible> FVMScheme::getResults() const
{
    return w;
}

Field<ThermoVar> FVMScheme::getResultsThermo() const
{
    return thermoField;
}

/*std::vector<std::vector<Compressible>> FVMScheme::calcBoundaryConditionsToBoundaryFields()
{
    const std::vector<int>& owners = mesh.getOwnerIndexList();
    const std::vector<Face>& faceList = mesh.getFaceList();
    std::vector<std::vector<Compressible>> out(boundaryConditionList.size());

    for (int boundaryConditionId = 0; boundaryConditionId < boundaryConditionList.size(); boundaryConditionId++)
    {
        out[boundaryConditionId] = boundaryConditionList[boundaryConditionId]->calc(w, mesh, thermo.get());
    }

    return out;
}*/